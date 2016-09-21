#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/GenericSelector.h"
#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/PATMuonTriggerMatchSelector.h"

typedef GenericObjectSelectorFilter<pat::MuonCollection,PATMuonTriggerMatchSelector,pat::MuonCollection> PATMuonTriggerMatchCloneSelectorFilter;
DEFINE_FWK_MODULE( PATMuonTriggerMatchCloneSelectorFilter );

