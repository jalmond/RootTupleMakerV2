#ifndef RootTupleMakerV2PFJets
#define RootTupleMakerV2PFJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_PFJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PFJets(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag , inputTagL1Offset;
  const edm::InputTag   inputTagSmearedUp, inputTagSmearedDown;
  const edm::InputTag   inputTagScaledUp, inputTagScaledDown;	 
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const std::string     jecUncPath;
  const bool            readJECuncertainty;
  const edm::InputTag   vtxInputTag;

  //OLD
  /*   const bool            applyResJEC; */
  /*   const std::string     resJEC; */
};

#endif
