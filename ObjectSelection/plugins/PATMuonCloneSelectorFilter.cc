#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/GenericSelector.h"
#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/PATMuonSelector.h"

//typedef edm::ThinningProducer<pat::MuonCollection, PATMuonSelector> ThinnedPATMuonProducer;
//DEFINE_FWK_MODULE(ThinnedPATMuonProducer);

typedef GenericObjectSelectorFilter<pat::MuonCollection,PATMuonSelector,pat::MuonCollection> PATMuonCloneSelectorFilter;
DEFINE_FWK_MODULE( PATMuonCloneSelectorFilter );


