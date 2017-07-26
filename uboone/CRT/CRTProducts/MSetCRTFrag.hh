#ifndef MSetCRTFrag_hh_
#define MSetCRTFrag_hh_
#include <cstdint>
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"
namespace crt
{
	class MSetCRTFrag
	{
	     public:
		MSetCRTFrag();
		MSetCRTFrag(artdaq::Fragment const& f1, artdaq::Fragment const& f2, artdaq::Fragment const& f3);
                virtual ~MSetCRTFrag();
	};
}

#endif
