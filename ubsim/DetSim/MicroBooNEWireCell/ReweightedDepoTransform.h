#ifndef UBSIM_DETSIM_MICROBOONEWIRECELL_REWEIGHTEDDEPOTRANSFORM
#define UBSIM_DETSIM_MICROBOONEWIRECELL_REWEIGHTEDDEPOTRANSFORM

#include "WireCellGen/DepoTransform.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"
#include "WireCellIface/WirePlaneId.h"
#include "WireCellIface/IDepo.h"

#include "TH2F.h"

namespace wcls {

    class ReweightedDepoTransform : public WireCell::Gen::DepoTransform,
			  public IArtEventVisitor
    {
    public:
	ReweightedDepoTransform();
	virtual ~ReweightedDepoTransform();

        /// IArtEventVisitor.
	//
	// Note: we don't actually poke at the event but use this
	// entry to refresh info from services in case they change
	// event-to-event.
        virtual void visit(art::Event & event);
	
        /// IConfigurable.
	//
	// Defer default to parent.  By default this class does not
	// overriding. 
	
    // Defer to parent but override if asked to.
        virtual void configure(const WireCell::Configuration& config);

        /// depo modifier
        virtual WireCell::IDepo::pointer modify_depo(WireCell::WirePlaneId wpid, WireCell::IDepo::pointer depo);
   
    private:
    std::string m_fileName_MC;
    std::vector<std::string> m_histnames;
    
    std::vector<TH2F*> m_hists;
    };

}  // wcls

#endif /* UBSIM_DETSIM_MICROBOONEWIRECELL_REWEIGHTEDDEPOTRANSFORM */

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
