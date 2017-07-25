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

//root
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TTree.h>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
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

  // ----------Tree and branches for mini-tuple ---------------------------
  TTree* eventTree;

  std::vector<float>* stub_pt;        // pt
  std::vector<float>* stub_eta;       // eta
  std::vector<float>* stub_phi;       // phi

  std::vector<int>*   nstub; //number of stubs
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
  stubSrc_ = iConfig.getParameter<edm::InputTag>("TTStubs");
  StubTok_ = consumes<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(stubSrc_);
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

  stub_pt->clear();
  stub_eta->clear();
  stub_phi->clear();
  nstub->clear();

  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > Phase2TrackerDigiTTStubHandle;
  iEvent.getByToken(StubTok_, Phase2TrackerDigiTTStubHandle);
  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;

  // Geometry
  /*  edm::ESHandle< TrackerGeometry > tGeometryHandle;
  const TrackerGeometry* theTrackerGeometry;
  iSetup.get< TrackerDigiGeometryRecord >().get( tGeometryHandle );
  theTrackerGeometry = tGeometryHandle.product();
  */
  /// Loop over input Stubs
  typename edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  typename edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator contentIter;
  //Adding protection
  if ( !Phase2TrackerDigiTTStubHandle.isValid() )  return;

  TH1* h1 = new TH1I("h1", "Number of Stubs", 11, 0.0, 11.0);
  TH1* h2 = new TH1I("h2", "Number of Stubs", 3, 0.0, 3.0);

  int temp, temp1 = 0;//, counter = 0;
  std::vector<int> stubPerEvent, Nstubs;

  for ( inputIter = Phase2TrackerDigiTTStubHandle->begin();
        inputIter != Phase2TrackerDigiTTStubHandle->end();
        ++inputIter )
    {  
      temp = 0;
      for ( contentIter = inputIter->begin(); contentIter != inputIter->end(); ++contentIter )
	{
	  /// Make reference stub
	  edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = edmNew::makeRefTo( Phase2TrackerDigiTTStubHandle, contentIter );
	  /*
	  /// Get det ID (place of the stub)
	  //  tempStubRef->getDetId() gives the stackDetId, not rawId
	  DetId detIdStub = theTrackerGeometry->idToDet( (tempStubRef->getClusterRef(0))->getDetId() )->geographicalId();

	  /// Define position stub by position inner cluster
	  MeasurementPoint mp = (tempStubRef->getClusterRef(0))->findAverageLocalCoordinates();
	  const GeomDet* theGeomDet = theTrackerGeometry->idToDet(detIdStub);
	  Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(mp) );
	 
	  double eta = posStub.eta();

	  stub_eta->push_back(eta);
	  */
	  temp++;
	  temp1 = 0;
	  temp1++;
	  Nstubs.push_back(temp1);
	}
      stubPerEvent.push_back(temp);
    }
  int vecSize = stubPerEvent.size();
  int vecSize2 = Nstubs.size();
    
  int i = 0;
  for ( inputIter = Phase2TrackerDigiTTStubHandle->begin();
        inputIter != Phase2TrackerDigiTTStubHandle->end();
        ++inputIter )
    {
      for (int j = 1; j < 10; j++) {
	if (stubPerEvent[i] == j)                                                                       
	  h1->Fill(stubPerEvent[i],j);
      } 
      i++;
      nstub->push_back(vecSize);
    }

  for (int k = 0; k <= vecSize2; k++)
    h2->Fill(Nstubs[k],1);
  
  //std::string intstr = std::to_string(vecSize);
  std::string intstr = std::to_string(vecSize2);
  TString name = intstr + "NStubs.pdf";
  TCanvas* hcanvas = new TCanvas("hcanvas","Canvas 1",100,100,800,800);
  h2->Draw("HIST");
  //gPad->SetLogy();
  hcanvas->SaveAs("../plugins/plots/"+name);
  delete hcanvas;

  eventTree->Fill();
}



// ------------ method called once each job just before starting event loop  ------------
void 
Phase2PixelStubs::beginJob()
{
  edm::Service<TFileService> fs;

  stub_pt = new std::vector<float>;
  stub_eta = new std::vector<float>;    
  stub_phi = new std::vector<float>;    

  nstub = new std::vector<int>;

  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");

  eventTree->Branch("stub_pt",     &stub_pt);
  eventTree->Branch("stub_eta",    &stub_eta);
  eventTree->Branch("stub_phi",    &stub_phi);
  eventTree->Branch("nstub",     &nstub);
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
