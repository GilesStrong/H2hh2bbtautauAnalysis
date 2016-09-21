#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/GenericSelector.h"
#include "H2hh2bbtautauAnalysis/ObjectSelection/interface/PATJetSelector.h"

//typedef edm::ThinningProducer<pat::JetCollection, PATJetSelector> ThinnedPATJetProducer;
//DEFINE_FWK_MODULE(ThinnedPATJetProducer);

typedef GenericObjectSelectorFilter<pat::JetCollection,PATJetSelector,pat::JetCollection> PATJetCloneSelectorFilter;
DEFINE_FWK_MODULE( PATJetCloneSelectorFilter );
