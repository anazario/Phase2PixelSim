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



// system include files
#include <memory>
#include <iostream>

// user include files
#include "CondFormats/HcalObjects/interface/HcalPedestals.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/typelookup.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerDetUnit.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DQM/Phase2OuterTracker/interface/OuterTrackerMonitorTTStub.h"

//root
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <TTree.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

typedef edm::Ref< edm::DetSetVector< Phase2TrackerDigi >, Phase2TrackerDigi > Ref_Phase2TrackerDigi_;

class DQMStore;

class Phase2PixelStubs : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Phase2PixelStubs(const edm::ParameterSet&);
      ~Phase2PixelStubs();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
 
  // ----------member data ---------------------------
  edm::InputTag stubSrc_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >  StubTok_;

  edm::InputTag TrackingParticleSrc_;
  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleTok_;

  edm::InputTag TrackingVertexSrc_;
  edm::EDGetTokenT< std::vector< TrackingVertex > > TrackingVertexTok_;

  edm::InputTag generalTracksSrc_;
  edm::EDGetTokenT< std::vector< reco::Track > > GeneralTracksTok_;

  // ----------Tree and branches for mini-tuple ---------------------------
  TTree* eventTree;

  std::vector<float>*   stub_pt;        // pt
  std::vector<float>*   stub_eta;       // eta
  std::vector<float>*   stub_phi;       // phi
  std::vector<float>*   trig_bend;
  std::vector<float>*   gentracks_eta;
  //number of stubs
  std::vector<int>*     nstub; //number of stubs per event
  std::vector<int>*     stub_barrel;
  std::vector<int>*     stub_endcap_disc;
  std::vector<int>*     stub_endcap_ring;
  std::vector<int>*     stub_endcap_ring_Fw;
  std::vector<int>*     stub_endcap_disc_Fw;
  std::vector<int>*     stub_endcap_ring_Bw;
  std::vector<int>*     stub_endcap_disc_Bw;
  //Global position of the stubs
  std::vector<float>*   stub_barrel_x;//stub x position TOB
  std::vector<float>*   stub_barrel_y;//stub y position TOB
  std::vector<float>*   stub_endcap_Fw_x;//stub x position TEC z>0
  std::vector<float>*   stub_endcap_Fw_y;//stub y position TEC z>0
  std::vector<float>*   stub_endcap_Fw_z;//stub z position TEC z>0
  std::vector<float>*   stub_endcap_Fw_r;//stub r position TEC z>0
  std::vector<float>*   stub_endcap_Bw_x;//stub x position TEC z<0                                                             
  std::vector<float>*   stub_endcap_Bw_y;//stub y position TEC z<0                                                                                       
  std::vector<float>*   stub_endcap_Bw_z;//stub z position TEC z<0                                                                                        
  std::vector<float>*   stub_endcap_Bw_r;//stub r position TEC z<0 
  //Stub displacement - offset
  std::vector<float>*   stub_barrel_w;//stub displacement - offset
  std::vector<float>*   stub_barrel_o;//stub offset TOB
  std::vector<float>*   stub_endcap_disc_w;
  std::vector<float>*   stub_endcap_disc_o;
  std::vector<float>*   stub_endcap_ring_w;
  std::vector<float>*   stub_endcap_ring_o;
  std::vector<float>*   stub_endcap_ring_w_Fw;
  std::vector<float>*   stub_endcap_ring_o_Fw;
  std::vector<float>*   stub_endcap_ring_w_Bw;
  std::vector<float>*   stub_endcap_ring_o_Bw;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//-------------------------------------------------------------------------------------------------
// constructors and destructor
//------------------------------------------------------------------------------------------------- 
Phase2PixelStubs::Phase2PixelStubs(const edm::ParameterSet& iConfig)
{
  stubSrc_ = iConfig.getParameter<edm::InputTag>("TTStubs");
  StubTok_ = consumes<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(stubSrc_);

  TrackingParticleSrc_ = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  TrackingParticleTok_ = consumes< std::vector< TrackingParticle > >(TrackingParticleSrc_);

  TrackingVertexSrc_ = iConfig.getParameter<edm::InputTag>("TrackingVertexInputTag");
  TrackingVertexTok_ = consumes< std::vector< TrackingVertex > >(TrackingVertexSrc_);

  generalTracksSrc_ = iConfig.getParameter<edm::InputTag>("GeneralTracksInputTag");
  GeneralTracksTok_ = consumes< std::vector< reco::Track > >(generalTracksSrc_);
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
  stub_barrel->clear();
  trig_bend->clear();
  nstub->clear();
  gentracks_eta->clear();
  stub_barrel->clear();
  stub_endcap_disc->clear();
  stub_endcap_ring->clear();
  stub_endcap_disc_Fw->clear();
  stub_endcap_ring_Fw->clear();
  stub_endcap_disc_Bw->clear();
  stub_endcap_ring_Bw->clear();
  stub_barrel_x->clear();
  stub_barrel_y->clear();
  stub_endcap_Fw_x->clear();
  stub_endcap_Fw_y->clear();
  stub_endcap_Fw_z->clear();
  stub_endcap_Fw_r->clear();
  stub_endcap_Bw_x->clear();
  stub_endcap_Bw_y->clear();
  stub_endcap_Bw_z->clear();
  stub_endcap_Bw_r->clear();
  stub_barrel_w->clear();
  stub_barrel_o->clear();
  stub_endcap_disc_w->clear();
  stub_endcap_disc_o->clear();
  stub_endcap_ring_w->clear();
  stub_endcap_ring_o->clear();
  stub_endcap_ring_w_Fw->clear();
  stub_endcap_ring_o_Fw->clear();
  stub_endcap_ring_w_Bw->clear();
  stub_endcap_ring_o_Bw->clear();  

  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > Phase2TrackerDigiTTStubHandle;
  iEvent.getByToken(StubTok_, Phase2TrackerDigiTTStubHandle);
  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
  
  /// Geometry
  edm::ESHandle<TrackerTopology> tTopoHandle;
  const TrackerTopology* tTopo;
  iSetup.get< TrackerTopologyRcd >().get(tTopoHandle);
  tTopo = tTopoHandle.product();
  
  //general tracks**
  //  edm::Handle< std::vector< reco::Track > > GeneralTracksHandle;
  //iEvent.getByToken(GeneralTracksTok_, GeneralTracksHandle);

  //tracking particles 
  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  edm::Handle< std::vector< TrackingVertex > > TrackingVertexHandle;
  iEvent.getByToken(TrackingParticleTok_, TrackingParticleHandle);
  iEvent.getByToken(TrackingVertexTok_, TrackingVertexHandle);
  
  edm::ESHandle< TrackerGeometry > tGeomHandle;
  iSetup.get< TrackerDigiGeometryRecord >().get(tGeomHandle);
  const TrackerGeometry* const theTrackerGeometry = tGeomHandle.product();
  
  /// Loop over input Stubs
  typename edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  typename edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator contentIter;
  //Adding protection
  if ( !Phase2TrackerDigiTTStubHandle.isValid() )  return;

  TH1* h1 = new TH1I("h1", "Number of Stubs", 11, 0.0, 11.0);
  TH1* h2 = new TH1I("h2", "Number of Stubs", 3, 0.0, 3.0);

  int temp, temp1 = 0;//, counter = 0;
  std::vector<int> stubPerEvent, Nstubs;

  // ----------------------------------------------------------------------------------------------
  // loop over L1 stubs
  // ----------------------------------------------------------------------------------------------

  for ( inputIter = Phase2TrackerDigiTTStubHandle->begin();
        inputIter != Phase2TrackerDigiTTStubHandle->end();
        ++inputIter )
    {  
      temp = 0;
      for ( contentIter = inputIter->begin(); contentIter != inputIter->end(); ++contentIter )
      {
        /// Make reference stub
	edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = edmNew::makeRefTo( Phase2TrackerDigiTTStubHandle, contentIter );

	/// Get det ID (place of the stub)
	//  tempStubRef->getDetId() gives the stackDetId, not rawId
	DetId detIdStub = theTrackerGeometry->idToDet( (tempStubRef->getClusterRef(0))->getDetId() )->geographicalId();
	
	/// Define position stub by position inner cluster
	MeasurementPoint mp = (tempStubRef->getClusterRef(0))->findAverageLocalCoordinates();
	const GeomDet* theGeomDet = theTrackerGeometry->idToDet(detIdStub);
	Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(mp) );
	  
	/// Get trigger displacement/offset
	double displStub = tempStubRef->getTriggerDisplacement();
	double offsetStub = tempStubRef->getTriggerOffset();

	double eta = posStub.eta();
	double phi = posStub.phi();

	float trigBend = tempStubRef->getTriggerBend();

	stub_phi->push_back(phi);		  
	stub_eta->push_back(eta);
	trig_bend->push_back(trigBend);

	//std::cout << "detIdStub.subdetId() = " << detIdStub.subdetId() << std::endl;

	if ( detIdStub.subdetId() == static_cast<int>(StripSubdetector::TOB) )  // Phase 2 Outer Tracker Barrel
	  {
	    stub_barrel->push_back(tTopo->layer(detIdStub)); 
	    stub_barrel_x->push_back(posStub.x());
	    stub_barrel_y->push_back(posStub.y());
	    stub_barrel_w->push_back(displStub - offsetStub); 
	    stub_barrel_o->push_back(offsetStub);
	  }

	else if ( detIdStub.subdetId() == static_cast<int>(StripSubdetector::TID) )  // Phase 2 Outer Tracker Endcap
	  {
	    int disc = tTopo->layer(detIdStub); // returns wheel
	    int ring = tTopo->tidRing(detIdStub);
	    stub_endcap_disc->push_back(disc);
	    stub_endcap_ring->push_back(ring);
	    stub_endcap_disc_w->push_back(displStub - offsetStub);
	    stub_endcap_ring_w->push_back(displStub - offsetStub);
	    stub_endcap_disc_o->push_back(offsetStub);
	    stub_endcap_ring_o->push_back(offsetStub);

	    if ( posStub.z() > 0 )
	      {
		stub_endcap_Fw_x->push_back( posStub.x());
		stub_endcap_Fw_y->push_back( posStub.y());
		stub_endcap_Fw_z->push_back( posStub.z());
		stub_endcap_Fw_r->push_back( posStub.perp() );
		stub_endcap_disc_Fw->push_back(disc);
		stub_endcap_ring_Fw->push_back(ring);
		stub_endcap_ring_w_Fw->push_back(displStub - offsetStub);
		stub_endcap_ring_o_Fw->push_back(offsetStub);
	      }
	    else
	      {
		stub_endcap_Bw_x->push_back( posStub.x());
		stub_endcap_Bw_y->push_back( posStub.y());
		stub_endcap_Bw_z->push_back( posStub.z());
		stub_endcap_Bw_r->push_back( posStub.perp() );
		stub_endcap_disc_Bw->push_back(disc);
		stub_endcap_ring_Bw->push_back(ring);
		stub_endcap_ring_w_Bw->push_back(displStub - offsetStub);
		stub_endcap_ring_o_Bw->push_back(offsetStub);
		}
	      }
	temp++;

	temp1 = 0;
	temp1++;
	Nstubs.push_back(temp1); //Actual number of stubs per event
      }
      stubPerEvent.push_back(temp); //more research needed
    }

  // ----------------------------------------------------------------------------------------------                     
  // loop over general tracks                                                                                             
  // ----------------------------------------------------------------------------------------------  
  /*
  std::vector< reco::Track >::const_iterator iterGTrk;

  for (iterGTrk = GeneralTracksHandle->begin(); iterGTrk != GeneralTracksHandle->end(); ++iterGTrk) {
    
    //edm::Ptr< reco::Track > gtrk_ptr(GeneralTracksHandle, this_gtrk);                                                  
    //this_gtrk++;

    float tmp_gtrk_eta = iterGTrk->eta();
    gentracks_eta->push_back(tmp_gtrk_eta);
  }
  */

  // ----------------------------------------------------------------------------------------------
  // loop over tracking particles
  // ----------------------------------------------------------------------------------------------
  /*
  int this_tp = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;

  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
 
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, this_tp);
    this_tp++;

    float track_pt = iterTP->pt();
    //if (track_pt < 2.0) continue;

    stub_pt->push_back(track_pt);
  }
  */

  //int vecSize = stubPerEvent.size();
  int vecSize2 = Nstubs.size();
  
  //For loops for local histogram creation (not on TTree)
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
      //nstub->push_back(vecSize);
    }

  for (int k = 0; k <= vecSize2; k++)
    h2->Fill(Nstubs[k],1);
  
  nstub->push_back(vecSize2); //fill stub per event in TTree
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
  trig_bend = new std::vector<float>;
  gentracks_eta = new std::vector<float>;
  nstub = new std::vector<int>;
  stub_barrel = new std::vector<int>;
  stub_endcap_disc = new std::vector<int>;
  stub_endcap_ring = new std::vector<int>;
  stub_endcap_disc_Fw = new std::vector<int>;
  stub_endcap_ring_Fw = new std::vector<int>;
  stub_endcap_disc_Bw = new std::vector<int>;
  stub_endcap_ring_Bw = new std::vector<int>;
  stub_barrel_x = new std::vector<float>;
  stub_barrel_y = new std::vector<float>;
  stub_endcap_Fw_x = new std::vector<float>;
  stub_endcap_Fw_y = new std::vector<float>;
  stub_endcap_Fw_z = new std::vector<float>;
  stub_endcap_Fw_r = new std::vector<float>;
  stub_endcap_Bw_x = new std::vector<float>;
  stub_endcap_Bw_y = new std::vector<float>;
  stub_endcap_Bw_z = new std::vector<float>;
  stub_endcap_Bw_r = new std::vector<float>;
  stub_barrel_w = new std::vector<float>;
  stub_barrel_o = new std::vector<float>;
  stub_endcap_disc_w = new std::vector<float>;
  stub_endcap_disc_o = new std::vector<float>;
  stub_endcap_ring_w = new std::vector<float>;
  stub_endcap_ring_o = new std::vector<float>;
  stub_endcap_ring_w_Fw = new std::vector<float>;
  stub_endcap_ring_o_Fw= new std::vector<float>;
  stub_endcap_ring_w_Bw = new std::vector<float>;
  stub_endcap_ring_o_Bw= new std::vector<float>;  
  // ntuple
  eventTree = fs->make<TTree>("eventTree", "Event tree");

  eventTree->Branch("stub_pt",     &stub_pt);
  eventTree->Branch("stub_eta",    &stub_eta);
  eventTree->Branch("stub_phi",    &stub_phi);
  eventTree->Branch("trig_bend",   &trig_bend);
  eventTree->Branch("nstub",     &nstub);
  eventTree->Branch("gentracks_eta",     &gentracks_eta);
  eventTree->Branch("stub_barrel", &stub_barrel);
  eventTree->Branch("stub_endcap_disc", &stub_endcap_disc);
  eventTree->Branch("stub_endcap_ring", &stub_endcap_ring);
  eventTree->Branch("stub_endcap_disc_Fw", &stub_endcap_disc_Fw);
  eventTree->Branch("stub_endcap_ring_Fw", &stub_endcap_ring_Fw);
  eventTree->Branch("stub_endcap_disc_Bw", &stub_endcap_disc_Bw);
  eventTree->Branch("stub_endcap_ring_Bw", &stub_endcap_ring_Bw);
  eventTree->Branch("stub_barrel_x", &stub_barrel_x);
  eventTree->Branch("stub_barrel_y", &stub_barrel_y);
  eventTree->Branch("stub_endcap_Fw_x", &stub_endcap_Fw_x);
  eventTree->Branch("stub_endcap_Fw_y", &stub_endcap_Fw_y);
  eventTree->Branch("stub_endcap_Fw_z", &stub_endcap_Fw_z);
  eventTree->Branch("stub_endcap_Fw_r", &stub_endcap_Fw_r);
  eventTree->Branch("stub_endcap_Bw_x", &stub_endcap_Bw_x);
  eventTree->Branch("stub_endcap_Bw_y", &stub_endcap_Bw_y);
  eventTree->Branch("stub_endcap_Bw_z", &stub_endcap_Bw_z);
  eventTree->Branch("stub_endcap_Bw_r", &stub_endcap_Bw_r);
  eventTree->Branch("stub_barrel_w", &stub_barrel_w);
  eventTree->Branch("stub_barrel_o", &stub_barrel_o);
  eventTree->Branch("stub_endcap_disc_w", &stub_endcap_disc_w);
  eventTree->Branch("stub_endcap_disc_o", &stub_endcap_disc_o);
  eventTree->Branch("stub_endcap_ring_w", &stub_endcap_ring_w);
  eventTree->Branch("stub_endcap_ring_o", &stub_endcap_ring_o);
  eventTree->Branch("stub_endcap_ring_w_Fw", &stub_endcap_ring_w_Fw);
  eventTree->Branch("stub_endcap_ring_o_Fw", &stub_endcap_ring_o_Fw);
  eventTree->Branch("stub_endcap_ring_w_Bw", &stub_endcap_ring_w_Bw);
  eventTree->Branch("stub_endcap_ring_o_Bw", &stub_endcap_ring_o_Bw);
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
