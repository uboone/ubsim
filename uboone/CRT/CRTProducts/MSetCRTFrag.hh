#ifndef MSetCRTFrag_hh_
#define MSetCRTFrag_hh_
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"

namespace crt
{
	class MSetCRTFrag
	{
		public:
		MSetCRTFrag();
		MSetCRTFrag(bernfebdaq::BernZMQFragment frag, unsigned tpc_s, unsigned tpc_ns, unsigned crt_frag_beg_ns, unsigned crt_frag_end_ns);
	        virtual ~MSetCRTFrag();
	};
}

#endif
