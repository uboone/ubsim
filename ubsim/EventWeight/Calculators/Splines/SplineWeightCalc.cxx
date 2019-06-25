/// @class evwgh::SplineWeightCalc
///
/// @brief Calculator for reweighting from one set of GENIE total cross section
/// splines to another
///
/// Eventually, this should be replaced with a tool in GENIE's own Reweight
/// product
///
/// @author Steven Gardiner <gardiner@fnal.gov> (24 June 2019)

// Standard library includes
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

// GENIE includes
// TODO: add preprocessor macro to check GENIE version, use v2 headers if needed
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/XmlParserUtils.h"

// libxml2 includes (needed for the LoadSplinesFromXML() member function)
#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/Exception.h"
#include "larsim/EventWeight/Base/WeightCalcCreator.h"
#include "larsim/EventWeight/Base/WeightCalc.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nutools/EventGeneratorBase/GENIE/GENIE2ART.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"

namespace {
  constexpr double BOGUS = -99999.;
}

namespace evwgh {

  class SplineWeightCalc : public WeightCalc {

    public:

      SplineWeightCalc() {}

      // Returns true if the splines were loaded successfully, or false if
      // there was a problem. The spline_map will be cleared and filled
      // with new splines from the requested XML file.
      bool LoadSplinesFromXML(const std::string& filename,
        std::map<std::string, std::unique_ptr<genie::Spline> >& spline_map);

      // evwgh::WeightCalc interface
      void Configure(fhicl::ParameterSet const& p,
        CLHEP::HepRandomEngine& engine);

      std::vector<std::vector<double > > GetWeight(art::Event& e);

    private:

     void GetSplines();

     /// Whether to check the total cross section from the event record
     /// against the old splines (true) or not (false)
     bool fCheckOldXSec;

     // Fractional tolerance for complaining when checking the stored
     // total cross section against the old splines
     bool fCheckOldToleranceFrac;

     /// Label used for the GENIEGen producer module
     std::string fGenieModuleLabel;

     /// Name of the XML file containing the "old" splines (used only if
     /// checking is enabled)
     std::string fOldFileName;

     /// Name of the XML file containing the "new" splines (compared to
     /// the total cross section in the event record to compute event weights)
     std::string fNewFileName;

     /// GENIE v3 tune name to use when retrieving the "old" splines
     std::string fOldTune;

     /// GENIE v3 tune name to use when retrieving the "new" splines
     std::string fNewTune;

     /// Container for the "old" splines. Keys are "tune/interaction" strings,
     /// values are splines
     std::map<std::string, std::unique_ptr<genie::Spline> > fOldSplineMap;

     /// Container for the "new" splines. Keys are "tune/interaction" strings,
     /// values are splines
     std::map<std::string, std::unique_ptr<genie::Spline> > fNewSplineMap;

     DECLARE_WEIGHTCALC(SplineWeightCalc)
  };

  void SplineWeightCalc::Configure(fhicl::ParameterSet const& p,
    CLHEP::HepRandomEngine& engine)
  {
    // Global config
    fGenieModuleLabel= p.get<std::string>("genie_module_label");

    // Config for this weight calculator
    fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>( this->GetName() );
    fNewFileName = pset.get<std::string>("new_splines_file");

    fOldTune = pset.get<std::string>("old_tune", "$GENIE_XSEC_TUNE");
    fNewTune = pset.get<std::string>("new_tune", "$GENIE_XSEC_TUNE");

    // Convert from an environment variable name (starting with '$') to the
    // corresponding value if needed.
    fOldTune = evgb::ExpandEnvVar( fOldTune );
    fNewTune = evgb::ExpandEnvVar( fNewTune );

    fCheckOldXSec = pset.get<bool>("check_old_xsec", false);
    if ( fCheckOldXSec ) {
      fOldFileName = pset.get<std::string>("old_splines_file");
    }

    fCheckOldToleranceFrac = pset.get<double>("check_old_frac_tol", 0.01);

    GetSplines();
  }

  void SplineWeightCalc::GetSplines() {

    cet::search_path sp("FW_SEARCH_PATH");

    // Don't bother to load the old splines if we're not going to use them
    // for checking
    if ( fCheckOldXSec ) {
      std::string old_file_path = sp.find_file( fOldFileName );
      bool old_ok = this->LoadSplinesFromXML( old_file_path, fOldSplineMap );
      if ( !old_ok ) mf::LogWarning("evwgh::SplineWeightCalc")
        << "Problem encountered while loading \"old\" splines from " << old_file_path;
    }

    std::string new_file_path = sp.find_file( fNewFileName );
    bool new_ok = this->LoadSplinesFromXML( new_file_path, fNewSplineMap );
    if ( !new_ok ) {
      throw cet::exception("evwgh::SplineWeightCalc")
        << "Problem encountered while loading \"new\" splines from " << new_file_path;
    }
  }

  // Returns a vector of weights for each neutrino interaction in the event
  std::vector<std::vector<double> > SplineWeightCalc::GetWeight(art::Event & e)
  {
    // Get truth-level information created by GENIE from the event
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;
    e.getByLabel(fGenieModuleLabel, mcTruthHandle);

    const art::FindOneP<simb::GTruth> gTruths(mcTruthHandle, e, fGenieModuleLabel);
    assert( gTruths.isValid() );

    // Initialize the vector of event weights. Each MCTruth will have a single
    // associated weight (multisims don't make sense for this calculator)
    std::vector< std::vector<double> > weights( mcTruthHandle->size() );

    for ( size_t mc_idx = 0; mc_idx < mcTruthHandle->size(); ++mc_idx ) {

      // Reconsitute the full GENIE event record using the current
      // MCTruth object and its associated GTruth object
      const simb::MCTruth& mc_truth = mcTruthHandle->at( mc_idx );
      const simb::GTruth& g_truth = *gTruths.at( mc_idx );

      // Note that the caller takes ownership of the event record produced by
      // evgb::RetrieveGHEP(). We wrap it in a std::unique_ptr here so that it
      // will be auto-deleted.
      std::unique_ptr<genie::EventRecord> ev_rec(evgb::RetrieveGHEP( mc_truth, g_truth ));

      // Translate the interaction mode into a string that can be used to look
      // up matching cross section splines
      std::string interaction_str = ev_rec->Summary()->AsString();

      // Incident neutrino energy (lab frame)
      double E_probe = mc_truth.GetNeutrino().Nu().E();

      std::string old_spline_key = fOldTune + '/' + interaction_str;
      std::string new_spline_key = fNewTune + '/' + interaction_str;

      double weight = 1.;

      bool have_old_spline = fOldSplineMap.count( old_spline_key );
      bool have_new_spline = fNewSplineMap.count( new_spline_key );

      // Get the (old) total cross section from the event record
      double total_xsec = ev_rec->XSec();

      if ( fCheckOldXSec && have_old_spline ) {
        const genie::Spline& old_spline = *fOldSplineMap.at( old_spline_key );
        double old_total_xsec = old_spline.Evaluate( E_probe );

        // Do the check and complain if it fails
        if ( std::abs( total_xsec - old_total_xsec ) / total_xsec > fCheckOldToleranceFrac ) {
          mf::LogWarning("evwgh::SplineWeightCalc") << "For a " << old_spline_key << " event, the"
            << " stored total cross section " << total_xsec << " GeV^(-2) differs from the one"
            << " calculated using the \"old\" splines " << old_total_xsec << " GeV^(-2) by more than "
            << fCheckOldToleranceFrac*100. << "%";
        }
      }

      // Only bother to calculate a new event weight if we have a new spline to
      // use for the chosen new tune and interaction type
      if ( have_new_spline ) {
        const genie::Spline& new_spline = *fNewSplineMap.at( new_spline_key );

        double new_total_xsec = new_spline.Evaluate( E_probe );

        // Protect against a NaN weight in case the total cross section
        // in the event record is zero. This shouldn't ever happen, but
        // it might if, e.g., the GTruth object wasn't loaded properly.
        if ( total_xsec != 0. ) weight = new_total_xsec / total_xsec;
        else {
          mf::LogWarning("evwgh::SplineWeightCalc") << " Stored total cross section is zero"
            << " for a " << interaction_str << " event. Setting weight to " << BOGUS;
          weight = BOGUS;
        }
      }

      weights[ mc_idx ].push_back( weight );
    }

    return weights;
  }

  // This function is almost entirely copied from the file
  // src/Framework/Utils/XSecSplineList.cxx from GENIE R-3_00_04.
  //
  // To avoid code duplication, it would have been much better to use the
  // genie::XSecSplineList class to manage loading the splines. However, since
  // this class is a singleton in the current version of GENIE, there could be
  // problems using this weight calculator simultaneously with an event
  // generation job or another weight calculator, e.g., race conditions while
  // accessing/modifying the global spline list. Rather than risk dealing with
  // those, I decided to copy the one function I needed here. In the future, I
  // think GENIE should allow this function to be used by other classes that
  // need to manipulate spline files. -- gardiner
  //
  // The LoadSplinesFromXML() function below is a derivative work based on code that is
  // Copyright (c) 2003-2019, The GENIE Collaboration
  // For the full text of the license visit http://copyright.genie-mc.org
  // or see $GENIE/LICENSE
  bool SplineWeightCalc::LoadSplinesFromXML(const std::string& filename,
    std::map<std::string, std::unique_ptr<genie::Spline> >& spline_map)
  {
    spline_map.clear();

    mf::LogInfo("evwgh::SplineWeightCalc") << "Loading splines from: " << filename;

    const int kNodeTypeStartElement = 1;
    const int kNodeTypeEndElement   = 15;
    const int kKnotX                = 0;
    const int kKnotY                = 1;

    xmlTextReaderPtr reader;

    int ret = 0, val_type = -1, iknot = 0, nknots = 0;
    double * E = 0, * xsec = 0;
    std::string spline_name;
    std::string temp_tune ;

    reader = xmlNewTextReaderFilename(filename.c_str());
    if (reader != NULL) {
      ret = xmlTextReaderRead(reader);
      while (ret == 1) {
        xmlChar * name  = xmlTextReaderName     (reader);
        xmlChar * value = xmlTextReaderValue    (reader);
        int       type  = xmlTextReaderNodeType (reader);
        int       depth = xmlTextReaderDepth    (reader);

        if (depth == 0 && type == kNodeTypeStartElement) {
           if(xmlStrcmp(name, (const xmlChar *) "genie_xsec_spline_list")) {
               mf::LogError("evwgh::SplineWeightCalc")
                 << "\nXML doc. has invalid root element! [filename: " << filename << "]";
               return false;
           }

           xmlChar * xvrs   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"version");
           xmlChar * xinlog = xmlTextReaderGetAttribute(reader,(const xmlChar*)"uselog");
           std::string svrs      = genie::utils::str::TrimSpaces((const char *)xvrs);
           std::string sinlog    = genie::utils::str::TrimSpaces((const char *)xinlog);

           mf::LogInfo("evwgh::SplineWeightCalc")
             << "Input x-section spline XML file format version: " << svrs;

           xmlFree(xvrs);
           xmlFree(xinlog);
        }

        if( (!xmlStrcmp(name, (const xmlChar *) "genie_tune")) && type==kNodeTypeStartElement) {
           xmlChar * xtune = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
           temp_tune    = genie::utils::str::TrimSpaces((const char *)xtune);
           mf::LogInfo("evwgh::SplineWeightCalc")
             << "Loading x-section splines for GENIE tune: " << temp_tune;
           xmlFree(xtune);
        }

        if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeStartElement) {
           xmlChar * xname = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
           xmlChar * xnkn  = xmlTextReaderGetAttribute(reader,(const xmlChar*)"nknots");
           std::string sname    = genie::utils::str::TrimSpaces((const char *)xname);
           std::string snkn     = genie::utils::str::TrimSpaces((const char *)xnkn);

           spline_name = sname;
           mf::LogInfo("evwgh::SplineWeightCalc")
             << "Loading spline: " << spline_name;

           nknots = atoi( snkn.c_str() );
           iknot=0;
           E     = new double[nknots];
           xsec  = new double[nknots];

           xmlFree(xname);
           xmlFree(xnkn);
        }

        if ( (!xmlStrcmp(name, (const xmlChar *) "E"))    && type==kNodeTypeStartElement) { val_type = kKnotX; }
        if ( (!xmlStrcmp(name, (const xmlChar *) "xsec")) && type==kNodeTypeStartElement) { val_type = kKnotY; }

        if ( (!xmlStrcmp(name, (const xmlChar *) "#text")) && depth==5) {
            if      (val_type==kKnotX) E   [iknot] = atof((const char *)value);
            else if (val_type==kKnotY) xsec[iknot] = atof((const char *)value);
        }
        if ( (!xmlStrcmp(name, (const xmlChar *) "knot")) && type==kNodeTypeEndElement) {
           iknot++;
        }
        if ( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeEndElement) {

          // Done looping over knots. Before creating the spline and putting it in the
          // map, build the spline key using the tune name and the string representation
          // of the interaction mode. We don't need the components of the algorithm ID
          // for our weight calculator.
          std::string spline_key = temp_tune;

          // The spline name format is algorithm_name / configuration / interaction
          std::vector<std::string> split_spline_name = genie::utils::str::Split(spline_name, "/");
          if ( split_spline_name.size() != 3 ) {
            mf::LogWarning("evwgh::SplineWeightCalc") << "Invalid spline name " << spline_name
              << " encountered!";
            return false;
          }
          std::string interaction_str = split_spline_name.back();

          // The new spline key format is tune / interaction
          spline_key += '/' + interaction_str;

          // Put the spline in the map. Note that this
          // will not overwrite an existing spline if one with the same key is
          // already present.
          spline_map.emplace( spline_key, new genie::Spline(nknots, E, xsec) );
          delete [] E;
          delete [] xsec;
        }
        xmlFree(name);
        xmlFree(value);
        ret = xmlTextReaderRead(reader);
      }
      xmlFreeTextReader(reader);
      if (ret != 0) {
        mf::LogError("evwgh::SplineWeightCalc")
          << "\nXML file could not be parsed! [filename: " << filename << "]";
        return false;
      }
    } else {
        mf::LogError("evwgh::SplineWeightCalc")
          << "\nXML file could not be found! [filename: " << filename << "]";
    }

    return true;
  }

  REGISTER_WEIGHTCALC(SplineWeightCalc)
}
