#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/GenericSelector.h"
#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/PATTauSelector.h"

//typedef edm::ThinningProducer<pat::TauCollection, PATTauSelector> ThinnedPATTauProducer;
//DEFINE_FWK_MODULE(ThinnedPATTauProducer);

typedef GenericObjectSelectorFilter<pat::TauCollection,PATTauSelector,pat::TauCollection> PATTauCloneSelectorFilter;
DEFINE_FWK_MODULE( PATTauCloneSelectorFilter );


