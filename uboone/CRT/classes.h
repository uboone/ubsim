/**
  \defgroup crt All things Cosmic Ray Tagger related
**/

#include "canvas/Persistency/Common/Wrapper.h"
#include "uboone/CRT/CRTSimData.hh"
#include "uboone/CRT/CRTAuxFunctions.hh"
#include <vector>

template class std::vector<crt::CRTSimData>;
template class art::Wrapper< std::vector<crt::CRTSimData> >;
template class art::Wrapper< crt::CRTSimData >;

template class std::vector<crt::CRTData::CRTHit>;
template class art::Wrapper< std::vector<crt::CRTData::CRTHit> >;
template class art::Wrapper< crt::CRTData::CRTHit >;


