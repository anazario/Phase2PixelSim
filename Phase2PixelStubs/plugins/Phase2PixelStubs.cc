// -*- C++ -*-
//
// Package:    Phase2PixelSim/Phase2PixelStubs
// Class:      Phase2PixelStubs
// 
/**\class Phase2PixelStubs Phase2PixelStubs.cc Phase2PixelSim/Phase2PixelStubs/plugins/Phase2PixelStubs.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  andres abreu nazario
//         Created:  Thu, 06 Jul 2017 21:24:55 GMT
//
//


// system include files
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "DQM/Phase2OuterTracker/interface/OuterTrackerMonitorTTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

typedef edm::Ref< edm::DetSetVector< Phase2TrackerDigi >, Phase2TrackerDigi > Ref_Phase2TrackerDigi_;

class Phase2PixelStubs : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Phase2PixelStubs(const edm::ParameterSet&);
      ~Phase2PixelStubs();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  TH1F * h_stub;

      // ----------member data ---------------------------
  edm::InputTag stubSrc_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >  StubTok_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Phase2PixelStubs::Phase2PixelStubs(const edm::ParameterSet& iConfig)

{
  edm::Service<TFileService> fs;
  h_stub = fs->make<TH1F>("StubAccepted", "StubAccepted", 5, 0, 5);

  stubSrc_ = iConfig.getParameter<edm::InputTag>("TTStubs");

  StubTok_ = consumes<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(stubSrc_);

  //produces<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >("StubAccepted");
}


Phase2PixelStubs::~Phase2PixelStubs()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Phase2PixelStubs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > Phase2TrackerDigiTTStubHandle;
  iEvent.getByToken(StubTok_, Phase2TrackerDigiTTStubHandle);
 
  /// Loop over input Stubs
  typename edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  typename edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator contentIter;
  //Adding protection
  if ( !Phase2TrackerDigiTTStubHandle.isValid() )  return;

  TH1* h1 = new TH1I("h1", "Number of Stubs", 8, 0.0, 8.0);

  int temp;//, counter = 0;
  std::vector<int> stubPerEvent;

  for ( inputIter = Phase2TrackerDigiTTStubHandle->begin();
        inputIter != Phase2TrackerDigiTTStubHandle->end();
        ++inputIter )
    {  
      temp = 0;
      for ( contentIter = inputIter->begin(); contentIter != inputIter->end(); ++contentIter )
	{
	  temp++;
	}
      stubPerEvent.push_back(temp);
    }
  int vecSize = stubPerEvent.size();
  std::cout << vecSize << std::endl;
  
  int i = 0;
  for ( inputIter = Phase2TrackerDigiTTStubHandle->begin();
        inputIter != Phase2TrackerDigiTTStubHandle->end();
        ++inputIter )
    {
      //std::cout << stubPerEvent[i];
      if (stubPerEvent[i] == 1)
	h1->Fill(stubPerEvent[i],1);
      if (stubPerEvent[i] == 2)
        h1->Fill(stubPerEvent[i],2);
      if (stubPerEvent[i] == 3)
        h1->Fill(stubPerEvent[i],3);
      if (stubPerEvent[i] == 4)
        h1->Fill(stubPerEvent[i],4);
      if (stubPerEvent[i] == 5)
        h1->Fill(stubPerEvent[i],5);
      if (stubPerEvent[i] == 6)
        h1->Fill(stubPerEvent[i],6);
      i++;
    }

  std::string intstr = std::to_string(vecSize);
  TString name = intstr + "NStubs.pdf";
  TCanvas* hcanvas = new TCanvas("hcanvas","Canvas 1",100,100,800,800);
  h1->Draw("HIST");
  hcanvas->SaveAs("../plugins/plots/"+name);
  delete hcanvas;
}


// ------------ method called once each job just before starting event loop  ------------
void 
Phase2PixelStubs::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Phase2PixelStubs::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Phase2PixelStubs::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2PixelStubs);
