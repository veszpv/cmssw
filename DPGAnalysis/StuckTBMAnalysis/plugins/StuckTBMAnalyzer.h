#ifndef STUCKTBMANALYZER_H
#define STUCKTBMANALYZER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/SiPixelDetId/interface/PixelFEDChannel.h"
#include "DQM/SiPixelPhase1Common/interface/SiPixelCoordinates.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

#ifndef NOVAL_I
#define NOVAL_I -9999
#endif
#ifndef NOVAL_F
#define NOVAL_F -9999.0
#endif

#define MAX_N_BAD_CHANNELS 5000  // should be something like 4736...

class StuckTBMAnalyzer : public edm::EDAnalyzer {
 public:
  StuckTBMAnalyzer(edm::ParameterSet const& iConfig);
  virtual ~StuckTBMAnalyzer();
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void getEventData(const edm::Event&, edm::Handle<PixelFEDChannelCollection>&);

 private:
  edm::ParameterSet iConfig_;

  SiPixelCoordinates coord_;

  std::string lumiFilename_;
  typedef std::map<int,std::pair<float,float> > LumiMap;
  typedef std::map<int,LumiMap> RunMap;
  RunMap puLumi_;

  std::string ntupleOutputFilename_;
  TFile* ntupleOutputFile_;
  
  int eventSaveDownscaling_;
  int nEvent_;

  TTree* eventTree_;
  TTree *permanentTree_;
  TTree *recoveredTree_;
  TTree *breakingTree_;

  class BadChannelData {
  public:
    unsigned int detid[MAX_N_BAD_CHANNELS];
    unsigned int fed[MAX_N_BAD_CHANNELS];
    unsigned int link[MAX_N_BAD_CHANNELS];
    unsigned int roc_first[MAX_N_BAD_CHANNELS];
    unsigned int roc_last[MAX_N_BAD_CHANNELS];

    int          det[MAX_N_BAD_CHANNELS];

    int          layer[MAX_N_BAD_CHANNELS];
    int          sec[MAX_N_BAD_CHANNELS];
    int          ladder[MAX_N_BAD_CHANNELS];
    int          module[MAX_N_BAD_CHANNELS];
    int          flipped[MAX_N_BAD_CHANNELS];

    int          disk[MAX_N_BAD_CHANNELS];
    int          ring[MAX_N_BAD_CHANNELS];
    int          blade[MAX_N_BAD_CHANNELS];

    int          eventInLS[MAX_N_BAD_CHANNELS]; // event num - first recorded event num
    int          timeInLS[MAX_N_BAD_CHANNELS];

    size_t       N;

    BadChannelData() { init(); }

    void init() {
      N=0;
    }

    int last() { return N==0 ? -1 : N-1; }

    void addBranches(TTree* tree, std::string sizevariable="") {
      tree->Branch("detid", detid, Form("detid[%s]/i", sizevariable.c_str()));
      tree->Branch("fed", fed, Form("fed[%s]/i", sizevariable.c_str()));
      tree->Branch("link", link, Form("link[%s]/i", sizevariable.c_str()));
      tree->Branch("roc_first", roc_first, Form("roc_first[%s]/i", sizevariable.c_str()));
      tree->Branch("roc_last", roc_last, Form("roc_last[%s]/i", sizevariable.c_str()));
      tree->Branch("det", det, Form("det[%s]/I", sizevariable.c_str()));
      tree->Branch("layer", layer, Form("layer[%s]/I", sizevariable.c_str()));
      tree->Branch("sec", sec, Form("sec[%s]/I", sizevariable.c_str()));
      tree->Branch("ladder", ladder, Form("ladder[%s]/I", sizevariable.c_str()));
      tree->Branch("module", module, Form("module[%s]/I", sizevariable.c_str()));
      tree->Branch("flipped", flipped, Form("flipped[%s]/I", sizevariable.c_str()));
      tree->Branch("disk", disk, Form("disk[%s]/I", sizevariable.c_str()));
      tree->Branch("ring", ring, Form("ring[%s]/I", sizevariable.c_str()));
      tree->Branch("blade", blade, Form("blade[%s]/I", sizevariable.c_str()));      
      tree->Branch("eventInLS", eventInLS, Form("eventInLS[%s]/I", sizevariable.c_str()));
      tree->Branch("timeInLS", timeInLS, Form("timeInLS[%s]/I", sizevariable.c_str()));
    }

    void insert(const BadChannelData& bc, unsigned int i) {
      if (i>=bc.N) {
	std::cout<<" *** ERROR: overindexing in source bad channel array\n";
	return;
      }
      if (N>=MAX_N_BAD_CHANNELS) {
	std::cout<<" *** ERROR: destination bad channel array is full\n";
	return;
      }
      detid[N]     = bc.detid[i];
      fed[N]       = bc.fed[i];
      link[N]      = bc.link[i];
      roc_first[N] = bc.roc_first[i];
      roc_last[N]  = bc.roc_last[i];
      det[N]       = bc.det[i];
      layer[N]     = bc.layer[i];
      sec[N]       = bc.sec[i];
      ladder[N]    = bc.ladder[i];
      module[N]    = bc.module[i];
      flipped[N]   = bc.flipped[i];
      disk[N]      = bc.disk[i];
      ring[N]      = bc.ring[i];
      blade[N]     = bc.blade[i];
      eventInLS[N] = bc.eventInLS[i];
      timeInLS[N]  = bc.timeInLS[i];
      N++;
    }
  };
  
  BadChannelData badchannels_;

  class EventData {
   public:
    
    int            fill;
    int            run;
    int            ls;
    int            bx;
    unsigned int   time;
    float          avgpu;
    float          lslumi;
    unsigned int   nbadchannels;
    int            eventInLS;
    // not saved:
    unsigned long long event;

    const std::string list = 
      "fill/I:run:ls:bx/I:time/i:avgpu/F:lslumi/F:nbadchannels/i:eventInLS/I";

    EventData() {
      init();
    }
    
    void init() {
      fill         = NOVAL_I;
      run          = NOVAL_I;
      ls           = NOVAL_I;
      bx           = NOVAL_I;
      time         = abs(NOVAL_I);
      avgpu        = NOVAL_F;
      lslumi       = NOVAL_F;
      nbadchannels = 0;
      eventInLS    = 0;
      event        = 0;
    }
    
  };

  EventData event_;

  class LSData {
   public:
    
    int          fill;
    int          run;
    int          ls;
    unsigned int time;
    float        avgpu;
    float        lslumi;
    unsigned int npermanent;
    unsigned int nrecovered;
    unsigned int nbreaking;
    // not saved:
    unsigned long long firstEventInLS;

    const std::string list = 
      "fill/I:run:ls:time/i:avgpu/F:lslumi/F:npermanent/i:nrecovered:nbreaking";

    LSData() {
      init();
    }
    
    void init() {
      fill         = NOVAL_I;
      run          = NOVAL_I;
      ls           = NOVAL_I;
      time         = abs(NOVAL_I);
      avgpu        = NOVAL_F;
      lslumi       = NOVAL_F;
      npermanent   = 0;
      nrecovered   = 0;
      nbreaking    = 0;
    }
    
  };
  
  LSData ls_;
  std::map<unsigned long long,std::pair<TH1I,TH1I> > badchannelsLSHists_;
  std::map<unsigned long long,BadChannelData> badchannelsLS_;
  BadChannelData permanentLS_;
  BadChannelData recoveredLS_;
  BadChannelData breakingLS_;

  const int n_of_FEDs_ = 139;
  const int n_of_Channels_ = 48;  
  const int fedId0_=1200;

  edm::EDGetTokenT<PixelFEDChannelCollection> pixelFEDChannelCollectionToken_;

};

#endif
