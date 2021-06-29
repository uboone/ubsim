#include "UbDrifter.h"

#include "larevt/CalibrationDBI/Interface/ElectronicsCalibService.h"
#include "larevt/CalibrationDBI/Interface/ElectronicsCalibProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeProvider.h"
#include "ubevt/Database/UbooneElectronLifetimeService.h"
#include "larcore/Geometry/Geometry.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"

#include <memory>

WIRECELL_FACTORY(wclsUbDrifter, wcls::UbDrifter,
		 wcls::IArtEventVisitor, WireCell::IDrifter)


using namespace WireCell;

wcls::UbDrifter::UbDrifter()
    : Drifter(), m_lifetime_to_set(1000.0*units::ms)
{
}

wcls::UbDrifter::~UbDrifter()
{
}

void wcls::UbDrifter::visit(art::Event & event)
{
    /// this function will be executed post WCT configuration !!!
    std::cout << "wcls UbDrifter: electron lifetime read in!\n"; 
    
    /// access lifetime database and fetch the value
    //m_lifetime_to_set = ??? from database 
    const lariov::UBElectronLifetimeProvider& elifetimeCalibProvider
        = art::ServiceHandle<lariov::UBElectronLifetimeService>()->GetProvider();
    float elifetime  = elifetimeCalibProvider.Lifetime(); // [ms]
    /// end 
    if (fELifetimeCorrection) {
      Drifter::set_lifetime(elifetime*units::ms);
      //std::cout << "www: lifetime from database = " << elifetime << std::endl;
      //std::cout << "www: lifetime from database = " << elifetime*units::ms << std::endl;
    }
    else {
      Drifter::set_lifetime(m_lifetime_to_set);
      //std::cout << "www: defaul lifetime = " << m_lifetime_to_set << std::endl;
      //std::cout << "www: defaul lifetime = " << 1000.0*units::ms << std::endl;
    }
}


void wcls::UbDrifter::configure(const WireCell::Configuration& cfg)
{
    Drifter::configure(cfg);
    fELifetimeCorrection = cfg["ELifetimeCorrection"].asBool();
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
