/// @class evwgh::COHNormWeightCalc
///
/// @brief Scales the coherent pion production total cross section up or down.
///        Together with generous 1-sigma uncertainties (50%), we hope to use
///        this as a stopgap replacement for GENIE's MaCOHpi and R0COHpi knobs,
///        which do not work correctly for our chosen tune in v3.0.6.
///
///        Ideally this would be implemented as an extra weight calculator in
///        GENIE itself, but given that it took a lot of effort to get GENIE
///        v3.0.4 uBooNE patch 01 out the door and into LArSoft/uboonecode,
///        this is a much faster workaround.
///
/// @author Steven Gardiner <gardiner@fnal.gov> (6 January 2020)

// Standard library includes
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

// GENIE includes
// TODO: add preprocessor macro to check GENIE version, use v2 headers if needed
#include "Framework/Interaction/Interaction.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/BaryonResonance.h"
#include "Framework/ParticleData/PDGCodes.h"

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

namespace evwgh {

  class COHNormWeightCalc : public WeightCalc {

    public:

      COHNormWeightCalc() {}

      // evwgh::WeightCalc interface
      void Configure(fhicl::ParameterSet const& p,
        CLHEP::HepRandomEngine& engine);

      std::vector<std::vector<double > > GetWeight(art::Event& e);

    private:

     /// Label used for the GENIEGen producer module
     std::string fGenieModuleLabel;

     /// FHiCL-configurable one-sigma uncertainty (as a fractional error)
     double fFractionalOneSigma;

     /// Variations in each universe (stored in units of one-sigma)
     std::vector<double> fSigmas;

     DECLARE_WEIGHTCALC(COHNormWeightCalc)
  };

  void COHNormWeightCalc::Configure(fhicl::ParameterSet const& p,
    CLHEP::HepRandomEngine& engine)
  {
    // Global config
    fGenieModuleLabel= p.get<std::string>("genie_module_label");

    // Config for this weight calculator
    fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>( this->GetName() );

    // Get the user-defined one-sigma for this weight calculator from the FHiCL
    // configuration
    fFractionalOneSigma = pset.get<double>("fractional_one_sigma");

    // Get the reweighting mode (multisim or ±1-sigma currently supported
    std::string mode = pset.get< std::string >( "mode" );

    fSigmas.clear();

    if ( mode.find("multisim") != std::string::npos ) {
      int num_universes = pset.get<int>( "number_of_multisims" );
      assert( num_universes >= 0 );
      for ( int u = 0; u < num_universes; ++u ) {
        fSigmas.push_back( CLHEP::RandGaussQ::shoot(&engine, 0., 1.) );
      }
    }
    else {
      // Assume ±1-sigma mode, do minus one-sigma in the first universe of two
      fSigmas = std::vector<double>( {-1., 1.} );
    }
  }

  // TODO: reduce code duplication between this function and SplineWeightCalc::GetWeight()
  // Returns a vector of weights for each neutrino interaction in the event
  std::vector<std::vector<double> > COHNormWeightCalc::GetWeight(art::Event & e)
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

      // Reconstitute the full GENIE event record using the current
      // MCTruth object and its associated GTruth object
      const simb::MCTruth& mc_truth = mcTruthHandle->at( mc_idx );
      const simb::GTruth& g_truth = *gTruths.at( mc_idx );

      // Note that the caller takes ownership of the event record produced by
      // evgb::RetrieveGHEP(). We wrap it in a std::unique_ptr here so that it
      // will be auto-deleted.
      std::unique_ptr<genie::EventRecord> ev_rec(evgb::RetrieveGHEP( mc_truth, g_truth ));

      // Do the actual checking of the GENIE event here

      // This calculator only applies to coherent pion production events
      // NOTE: Be careful. This function will change its name for GENIE v3.2
      // (which includes coherent elastic neutrino nucleus scattering)
      genie::Interaction* interaction = ev_rec->Summary();
      bool is_coh = interaction->ProcInfo().IsCoherent();

      // Compute weights for each universe
      for ( size_t u = 0u; u < fSigmas.size(); ++u ) {

        // Default to unit weight unless we're dealing with coherent pion
        // production events
        double weight = 1.;
        if ( is_coh ) weight = 1. + fFractionalOneSigma * fSigmas.at(u);

        // Don't let the weight ever go negative (unphysical)
        weight = std::max( weight, 0. );

        // Save the weight for the current universe
        weights[ mc_idx ].push_back( weight );
      }
    }

    return weights;
  }

  REGISTER_WEIGHTCALC(COHNormWeightCalc)
}
