#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/GenericSelector.h"
#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/PATElectronSelector.h"

typedef GenericObjectSelectorFilter<pat::ElectronCollection,PATElectronSelector,pat::ElectronCollection> PATElectronCloneSelectorFilter;
DEFINE_FWK_MODULE( PATElectronCloneSelectorFilter );


