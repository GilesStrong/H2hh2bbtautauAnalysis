// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelectorStreamProducer.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/UtilAlgos/interface/NonNullNumberSelector.h"
#include "CommonTools/UtilAlgos/interface/StoreManagerTrait.h"
#include "CommonTools/UtilAlgos/interface/StoreContainerTrait.h"
#include "CommonTools/UtilAlgos/interface/SelectedOutputCollectionTrait.h"
#include "CommonTools/UtilAlgos/interface/SelectionAdderTrait.h"
#include "CommonTools/UtilAlgos/interface/NullPostProcessor.h"
#include "CommonTools/UtilAlgos/interface/EventSetupInitTrait.h"

#include "DataFormats/Common/interface/View.h"

//Selector template class (providing full machinery of collection loop...)
template<typename InputCollection,
         typename Selector,
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type,
         typename StoreContainer = typename ::helper::StoreContainerTrait<OutputCollection>::type,
         typename RefAdder = typename ::helper::SelectionAdderTrait<InputCollection, StoreContainer>::type>
struct SelectorTemplate {
  typedef InputCollection collection;
  typedef StoreContainer container;
  typedef typename container::const_iterator const_iterator;

  SelectorTemplate ( const edm::ParameterSet & cfg, edm::ConsumesCollector&& cc ) :
    selector_(cfg,cc) { }
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle<collection> & c , const edm::Event & ed, const edm::EventSetup & es) {
    selected_.clear();
    selector_.preChoose(ed, es);
    for(size_t idx = 0; idx < c->size(); ++ idx) {
      if(selector_.choose(idx, (*c)[idx]))
        addRef_(selected_, c, idx);
    }
  }
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  RefAdder addRef_;
  Selector selector_;
};


template<typename InputCollection,
         typename Selector,
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type,
         typename StoreContainer = typename ::helper::StoreContainerTrait<OutputCollection>::type,
         typename RefAdder = typename ::helper::SelectionAdderTrait<InputCollection, StoreContainer>::type>
struct SelectorRefTemplate {
  typedef InputCollection collection;
  typedef StoreContainer container;
  typedef typename container::const_iterator const_iterator;

  SelectorRefTemplate ( const edm::ParameterSet & cfg, edm::ConsumesCollector&& cc ) :
    selector_(cfg,cc) { }
  const_iterator begin() const { return selected_.begin(); }
  const_iterator end() const { return selected_.end(); }
  void select(const edm::Handle<collection> & c , const edm::Event & ed, const edm::EventSetup & es) {
    selected_.clear();
    selector_.preChoose(ed, es);
    for(size_t idx = 0; idx < c->size(); ++ idx) {
      if(selector_.choose(idx, *((*c)[idx])))
        addRef_(selected_, c, idx);
    }
  }
  size_t size() const { return selected_.size(); }
private:
  container selected_;
  RefAdder addRef_;
  Selector selector_;
};


//little redefinition to make the use easier
template<typename InputCollection,
         typename Selector,
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type>
using GenericObjectSelectorFilter = ObjectSelector<SelectorTemplate<InputCollection,Selector,OutputCollection>,OutputCollection>;

template<typename InputCollection,
         typename Selector,
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type>
using GenericRefObjectSelectorFilter = ObjectSelector<SelectorRefTemplate<InputCollection,Selector,OutputCollection>,OutputCollection>;

template<typename InputCollection,
         typename Selector,
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type>
using GenericObjectSelectorProducer = ObjectSelectorStreamProducer<SelectorTemplate<InputCollection,Selector,OutputCollection>,OutputCollection>;

template<typename InputCollection,
         typename Selector,
         typename OutputCollection = typename ::helper::SelectedOutputCollectionTrait<InputCollection>::type>
using GenericRefObjectSelectorProducer = ObjectSelectorStreamProducer<SelectorRefTemplate<InputCollection,Selector,OutputCollection>,OutputCollection>;
