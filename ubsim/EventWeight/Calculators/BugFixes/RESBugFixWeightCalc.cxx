/// @class evwgh::RESBugFixWeightCalc
///
/// @brief Zeros out RES events that produced a P33(1600) or F17(1970) resonance
///        which was mislabeled as a "rootino" (PDG code == 0) in the GENIE event
///        record. This is one way of fixing a bug in GENIE v3.0.6 that was first
///        reported during the October 2019 User Forum (see GENIE docDB #153,
///        https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=153).
///
/// @author Steven Gardiner <gardiner@fnal.gov> (27 December 2019)

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

namespace {
  constexpr int ROOTINO = 0; // PDG code zero
}

namespace evwgh {

  class RESBugFixWeightCalc : public WeightCalc {

    public:

      RESBugFixWeightCalc() {}

      // evwgh::WeightCalc interface
      void Configure(fhicl::ParameterSet const& p,
        CLHEP::HepRandomEngine& engine);

      std::vector<std::vector<double > > GetWeight(art::Event& e);

    private:

     /// Label used for the GENIEGen producer module
     std::string fGenieModuleLabel;

     DECLARE_WEIGHTCALC(RESBugFixWeightCalc)
  };

  void RESBugFixWeightCalc::Configure(fhicl::ParameterSet const& p,
    CLHEP::HepRandomEngine& /*engine*/)
  {
    // Global config
    fGenieModuleLabel= p.get<std::string>("genie_module_label");

    // Config for this weight calculator
    //fhicl::ParameterSet const& pset = p.get<fhicl::ParameterSet>( this->GetName() );
  }

  // TODO: reduce code duplication between this function and SplineWeightCalc::GetWeight()
  // Returns a vector of weights for each neutrino interaction in the event
  std::vector<std::vector<double> > RESBugFixWeightCalc::GetWeight(art::Event & e)
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
      double weight = 1.;
      genie::Interaction* interaction = ev_rec->Summary();

      // The bug only occurs for resonant events
      bool is_res = interaction->ProcInfo().IsResonant();

      // Only events that produce a P33(1600) or F17(1970) are affected
      genie::Resonance_t res_code = interaction->ExclTag().Resonance();
      bool bad_resonance = ( res_code == genie::kP33_1600
        || res_code == genie::kF17_1970 );

      // Either of these two resonances are mislabeled by the bug
      // as a "rootino" (PDG code == 0) in the event record. If we find
      // this rootino (confirming the presence of the bug during event
      // generation), then zero out the weight for the current event.
      genie::GHepParticle* mislabeled_resonance = ev_rec->FindParticle( ROOTINO,
        genie::kIStPreDecayResonantState, 0 );

      if ( is_res && bad_resonance && mislabeled_resonance ) {

        weight = 0.;

        mf::LogInfo("evwgh::RESBugFixWeightCalc") << "Assigned a weight of zero"
          << " to a buggy RES event";
      }

      weights[ mc_idx ].push_back( weight );
    }

    return weights;
  }

  REGISTER_WEIGHTCALC(RESBugFixWeightCalc)
}
