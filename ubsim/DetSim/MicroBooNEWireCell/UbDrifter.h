#ifndef UBSIM_DETSIM_MICROBOONEWIRECELL_UBDRIFTER
#define UBSIM_DETSIM_MICROBOONEWIRECELL_UBDRIFTER

#include "WireCellGen/Drifter.h"
#include "larwirecell/Interfaces/IArtEventVisitor.h"

namespace wcls {

    class UbDrifter : public WireCell::Gen::Drifter,
			  public IArtEventVisitor
    {
    public:
	UbDrifter();
	virtual ~UbDrifter();

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
   
    private:
        double m_lifetime_to_set;
    };

}  // wcls

#endif /* UBSIM_DETSIM_MICROBOONEWIRECELL_UBDRIFTER */

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
