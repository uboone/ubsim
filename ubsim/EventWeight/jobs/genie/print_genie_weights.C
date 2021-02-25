// Make sure to set up the gallery ups product before using this macro
#include <iostream>
#include <string>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

void print_genie_weights( const std::string& filename ) {

  std::vector<std::string> filenames { filename };

  gallery::Event ev( filenames );

  size_t event_count = 0;
  for ( ev.toBegin(); !ev.atEnd(); ++ev ) {

    std::cout << "art event " << event_count << '\n';

    // Use the GENIE producer label here
    auto mctruth_handle = ev.getValidHandle< std::vector<simb::MCTruth> >( "generator" );
    auto gtruth_handle = ev.getValidHandle< std::vector<simb::GTruth> >( "generator" );

    // Use the EventWeight producer label here
    auto weights_handle = ev.getValidHandle< std::vector<evwgh::MCEventWeight> >( "eventweight" );

    // Loop through these objects for each neutrino vertex in the event
    for ( size_t v = 0u; v < mctruth_handle->size(); ++v ) {

      const simb::MCTruth& mc = mctruth_handle->at( v );
      const simb::GTruth& gt = gtruth_handle->at( v );
      const evwgh::MCEventWeight& mc_weights = weights_handle->at( v );

      // Extract desired information
      const auto& nu = mc.GetNeutrino().Nu();
      if ( nu.NumberTrajectoryPoints() > 0 ) {
        double E_nu = nu.E( 0 );
        std::cout << "Neutrino energy = " << E_nu << " GeV\n";
      }

      // Loop over all of the weights in the MCEventWeight object
      std::cout << "Weights\n";
      for ( const auto& pair : mc_weights.fWeight ) {
        std::string knob_name = pair.first;
        std::vector<double> weights = pair.second;

        std::cout << "  " << knob_name << ":\n";
        for ( size_t u = 0u; u < weights.size(); ++u ) {
          double w = weights.at( u );
          std::cout << "    universe #" << u << " has weight = " << w << '\n';
        }
      }
    }

    ++event_count;
  }
}
