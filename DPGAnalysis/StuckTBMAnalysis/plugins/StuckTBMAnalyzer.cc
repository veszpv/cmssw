#include "StuckTBMAnalyzer.h"


StuckTBMAnalyzer::StuckTBMAnalyzer(edm::ParameterSet const& iConfig) : 
  iConfig_(iConfig),
  lumiFilename_(iConfig.getUntrackedParameter<std::string>
		("lumiFileName", "")),
  ntupleOutputFilename_(iConfig.getUntrackedParameter<std::string>
			("outputFileName", "Ntuple.root")),
  eventSaveDownscaling_(iConfig.getUntrackedParameter<int>("eventSaveDownscaleFactor",1)) {

  nEvent_=0;
  pixelFEDChannelCollectionToken_ = 
    consumes<PixelFEDChannelCollection>(edm::InputTag("siPixelDigis"));

  std::ifstream infile(lumiFilename_);
  std::string line;
  while (std::getline(infile, line)) {
    std::cout<<line<<std::endl;
    if (std::count(line.begin(), line.end(), ',')!=8) continue;
    //run:fill,ls,time,beamstatus,E(GeV),delivered(/ub),recorded(/ub),avgpu,source
    //300806:6061,32:32,08/09/17 21:45:48,STABLE BEAMS,6500,403820.866,262944.928,48.5,HFET

    int run = atoi(line.substr(0, line.find_first_of(':')).c_str());
    size_t p1 = line.find_first_of(',')+1;
    size_t p2 = line.find_first_of(':', p1);
    int ls = atoi(line.substr(p1, p2-p1).c_str());
    for (int i=0; i<4; i++) p1 = line.find_first_of(',', p1)+1;
    p2 = line.find_first_of(',', p1);
    float lumi = atof(line.substr(p1, p2-p1).c_str());
    p1 = line.find_first_of(',', p2+1)+1;
    float avgpu = atof(line.substr(p1, line.find_first_of(',', p1)-p1).c_str());

    std::cout<<"Run: "<<run<<" LS: "<<ls<<" Lumi: "<<lumi<<" PU: "<<avgpu<<std::endl;
    
    RunMap::iterator it=puLumi_.find(run);
    if (it==puLumi_.end()) {
      LumiMap lsmap;
      lsmap[ls]=std::pair<float,float>(avgpu,lumi);
      puLumi_[run]=lsmap;
    } else {
      LumiMap::iterator it2=it->second.find(ls);
      if (it2==it->second.end()) {
	it->second[ls]=std::pair<float,float>(avgpu,lumi);
      } else {
	std::cout<<"*** ERROR: run+LS has already been recorded...\n";
      }
    }
  }

}


StuckTBMAnalyzer::~StuckTBMAnalyzer() {}


void StuckTBMAnalyzer::beginJob() {

  ntupleOutputFile_ = new TFile(ntupleOutputFilename_.c_str(), "RECREATE");

  if(!(ntupleOutputFile_ -> IsOpen())) {
    std::cout << "Cannot open output file: \"" << ntupleOutputFilename_ << "\""<< std::endl;
  } else {
    std::cout << "Output file: \"" << ntupleOutputFilename_ << "\" created." << std::endl;
  }

  eventTree_= new TTree("eventTree", "The events");
  eventTree_->Branch("event", &event_, event_.list.c_str());
  badchannels_.addBranches(eventTree_, "event.nbadchannels");
  
  permanentTree_= new TTree("permanentTree", "Channels bad for the whole lumisection");
  permanentTree_->Branch("ls", &ls_, ls_.list.c_str());
  permanentLS_.addBranches(permanentTree_, "ls.npermanent");

  recoveredTree_= new TTree("recoveredTree", "Channels recovered during the lumisection");
  recoveredTree_->Branch("ls", &ls_, ls_.list.c_str());
  recoveredLS_.addBranches(recoveredTree_, "ls.nrecovered");

  breakingTree_= new TTree("breakingTree", "Channels breaking during the lumisection");
  breakingTree_->Branch("ls", &ls_, ls_.list.c_str());
  breakingLS_.addBranches(breakingTree_, "ls.nbreaking");

}


void StuckTBMAnalyzer::endJob() {
  std::cout << "Ntuplizer endjob step with outputFileName: \"" 
	    << ntupleOutputFilename_ << "\"." << std::endl;
  ntupleOutputFile_ -> Write();
  std::cout << "Closing file: \"" << ntupleOutputFilename_ << "\"." << std::endl;
  ntupleOutputFile_ -> Close();
}


void StuckTBMAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi,
						edm::EventSetup const& iSetup) {
  ls_.init();
  breakingLS_.init();
  recoveredLS_.init();
  permanentLS_.init();
  badchannelsLSHists_.clear();
  badchannelsLS_.clear();
}


void StuckTBMAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& iLumi,
					      edm::EventSetup const& iSetup) {

  TH1I permanent("permanent", "Each bin is a FED link, links always bad in LS", 
		 n_of_FEDs_*n_of_Channels_, 0,  n_of_FEDs_*n_of_Channels_);

  std::vector<int> breaking(n_of_FEDs_*n_of_Channels_, 0);
  std::vector<unsigned long long> breaking_event(n_of_FEDs_*n_of_Channels_, 0);

//   TH1I breaking("breaking", "Each bin is a FED link, links breaking and never recovering", 
// 		n_of_FEDs_*n_of_Channels_, 0,  n_of_FEDs_*n_of_Channels_);
  
//   TH1I breaking_event("breaking_event", "Eqach bin is a FED link, when channel broke", 
// 		      n_of_FEDs_*n_of_Channels_, 0,  n_of_FEDs_*n_of_Channels_);

  std::map<unsigned long long,std::pair<TH1I,TH1I> >::iterator it_hists = badchannelsLSHists_.begin();
  std::map<unsigned long long,std::pair<TH1I,TH1I> >::iterator it_hists_prev = it_hists;

  // collecting all bad channels in LS
  for (; it_hists!=badchannelsLSHists_.end(); it_hists++) {
    if (it_hists==badchannelsLSHists_.begin()) {
      permanent=it_hists->second.first;
      continue; // rest of for loop reached only from the second event in the LS
    } else {
      permanent.Multiply(&it_hists->second.first);
    }

    //int eventInLS=0;
    if (it_hists->first<=it_hists_prev->first) { // the "<=" is OK since the "continue" above...
      std::cout<<" *** ERROR: events are not ordered by number in maps\n";
    } else {
      //eventInLS = it_hists->first - badchannelsLSHists_.begin()->first;
    }

    TH1I h=it_hists->second.first;
    h.Add(&it_hists_prev->second.first, -1.0);

    const TH1I& ind=it_hists->second.second;
    const BadChannelData& bcd=badchannelsLS_[it_hists->first];

    for (int i=1; i<=h.GetNbinsX(); i++) {
      int diff=h.GetBinContent(i);
      if (diff==0) continue;

      //int n=breaking.GetBinContent(i);
      int& n=breaking[i-1];
      if (diff==1) { // just breaking
	if (n!=0) std::cout<<" *** ERROR: channel breaking twice\n";
// 	breaking.SetBinContent(i, n+1);
// 	breaking_event.SetBinContent(i, it_hists-badchannelsLSHists_.begin());
	n++;
	breaking_event[i-1]=it_hists->first;
	//breakingLS_.insert(bcd, ind.GetBinContent(i));
	//ls_.nbreaking++;
	continue;
      }
      if (diff==-1) { // just recovering
	if (n==0) { // broke in earlier LS
	} else { // broke in this LS
	  if (n>1) std::cout<<" *** ERROR: channel was broken twice, but now recovering\n";
	  //breaking.SetBinContent(i, n-1);
	  n--;
	}
	recoveredLS_.insert(bcd, ind.GetBinContent(i));
	//recoveredLS_.eventInLS[recoveredLS_.last()]=eventInLS;
	ls_.nrecovered++;
	continue;
      }
      std::cout<<" *** ERROR: diffing, duplicate bad channel in event? \n";
    }

    it_hists_prev=it_hists;
  }

  // saving permanently bad and just broken but not recovered channels in LS
  for (int i=1; i<=permanent.GetNbinsX(); i++) {
    if (permanent.GetBinContent(i)!=0) {
      permanentLS_.insert(badchannelsLS_.begin()->second, // 1st event OK, was bad in entire LS
			  badchannelsLSHists_.begin()->second.second.GetBinContent(i));
      ls_.npermanent++;
    }
    if (breaking[i-1]!=0) {
      breakingLS_.insert(badchannelsLS_[breaking_event[i-1]], 
			 badchannelsLSHists_[breaking_event[i-1]].second.GetBinContent(i));
//       breakingLS_.eventInLS[breakingLS_.last()] = 
// 	breaking_event[i-1] - badchannelsLSHists_.begin()->first;
      ls_.nbreaking++;
    }
  }

  assert(breakingLS_.N==ls_.nbreaking);
  assert(recoveredLS_.N==ls_.nrecovered);
  assert(permanentLS_.N==ls_.npermanent);

  // Fill
  permanentTree_->Fill();
  recoveredTree_->Fill();
  breakingTree_->Fill();

}


void StuckTBMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (++nEvent_ % eventSaveDownscaling_ != 0) return;
  
  edm::Handle<PixelFEDChannelCollection> pixelFEDChannelCollectionHandle;
  iEvent.getByToken(pixelFEDChannelCollectionToken_, pixelFEDChannelCollectionHandle);

  coord_.init(iSetup);

  getEventData(iEvent, pixelFEDChannelCollectionHandle);
}


void StuckTBMAnalyzer::getEventData(const edm::Event& iEvent, 
				    edm::Handle<PixelFEDChannelCollection>& 
				    pixelFEDChannelCollectionHandle) {

  event_.init();
  // badchannel_.init(); not needed as long as event_.nbadchannels is set to zero

  event_.fill         = iEvent.eventAuxiliary().storeNumber();
  event_.run          = iEvent.id().run();
  event_.ls           = iEvent.luminosityBlock();
  event_.event        = iEvent.id().event();
  //std::cout << "Event number: " << event_.event << std::endl;
  event_.bx           = iEvent.bunchCrossing();
  event_.time         = iEvent.time().unixTime();

  if (ls_.ls==NOVAL_I) { // in first event of the LS
    RunMap::const_iterator it=puLumi_.find(event_.run);
    if (it!=puLumi_.end()) {
      LumiMap::const_iterator it2=it->second.find(event_.ls);
      if (it2!=it->second.end()) {
	event_.avgpu        = it2->second.first;
	event_.lslumi       = it2->second.second;
      }
    }
    ls_.fill=event_.fill;
    ls_.run=event_.run;
    ls_.ls=event_.ls;
    ls_.time=event_.time;
    ls_.avgpu=event_.avgpu;
    ls_.lslumi=event_.lslumi;
    ls_.firstEventInLS=event_.event;
  } else {
    event_.avgpu = ls_.avgpu;
    event_.lslumi = ls_.lslumi;
    event_.eventInLS = event_.event - ls_.firstEventInLS;
  }

  TH1I channels("channels", "Each bin is a FED link", 
		n_of_FEDs_*n_of_Channels_, 0,  n_of_FEDs_*n_of_Channels_);
  
  TH1I indices("indices", "Each bin is a FED link", 
	       n_of_FEDs_*n_of_Channels_, 0,  n_of_FEDs_*n_of_Channels_);

  // Loon on bad channels
  event_.nbadchannels=0;
  badchannels_.N=0;
  //std::cout<<"Loop on bad channels\n";
  for(const auto& disabledOnDetId: *pixelFEDChannelCollectionHandle) {
    for(const auto& disabledChannels: disabledOnDetId) {
//       std::cout<<" DetId: " <<disabledOnDetId.detId()
// 	       <<" fed: "<<disabledChannels.fed
// 	       <<" link: "<<disabledChannels.link
// 	       <<" roc_first: "<< disabledChannels.roc_first
//	       <<" roc_last: "<<disabledChannels.roc_last<<std::endl;
      badchannels_.detid[event_.nbadchannels]=disabledOnDetId.detId();
      badchannels_.fed[event_.nbadchannels]=disabledChannels.fed;
      badchannels_.link[event_.nbadchannels]=disabledChannels.link;
      badchannels_.roc_first[event_.nbadchannels]=disabledChannels.roc_first;
      badchannels_.roc_last[event_.nbadchannels]=disabledChannels.roc_last;

      Float_t bin=float(disabledChannels.fed-fedId0_)*n_of_Channels_+disabledChannels.link;

      if (channels.GetBinContent(channels.FindBin(bin))!=0) {
	std::cout<<"*** Error: FED link in error is listed multiple times";
      } else {
	channels.Fill(bin);
	indices.Fill(bin, event_.nbadchannels);
      }

      DetId detId(disabledOnDetId.detId());

      if (detId.subdetId() == PixelSubdetector::PixelBarrel) {
	badchannels_.det[event_.nbadchannels]=0;
	badchannels_.layer[event_.nbadchannels]=coord_.layer(detId);
	badchannels_.sec[event_.nbadchannels]=coord_.sector(detId);
	badchannels_.ladder[event_.nbadchannels]=coord_.signed_ladder(detId);
	badchannels_.module[event_.nbadchannels]=coord_.signed_module(detId);
	badchannels_.flipped[event_.nbadchannels]=coord_.flipped(detId);
	badchannels_.disk[event_.nbadchannels]=NOVAL_I;
	badchannels_.ring[event_.nbadchannels]=NOVAL_I;
	badchannels_.blade[event_.nbadchannels]=NOVAL_I;
      } else {
	badchannels_.det[event_.nbadchannels]=1;
	badchannels_.layer[event_.nbadchannels]=NOVAL_I;
	badchannels_.sec[event_.nbadchannels]=NOVAL_I;
	badchannels_.ladder[event_.nbadchannels]=NOVAL_I;
	badchannels_.module[event_.nbadchannels]=coord_.module(detId);
	badchannels_.flipped[event_.nbadchannels]=NOVAL_I;
	badchannels_.disk[event_.nbadchannels]=coord_.signed_disk(detId);
	badchannels_.ring[event_.nbadchannels]=coord_.ring(detId);
	badchannels_.blade[event_.nbadchannels]=coord_.signed_blade(detId);;
      }

      badchannels_.eventInLS[event_.nbadchannels]=event_.eventInLS;
      badchannels_.timeInLS[event_.nbadchannels]=ls_.time-event_.time;

      badchannels_.N++;
      event_.nbadchannels++;
    }
  }

//   for (int i=0; i<event_.nbadchannels; i++) {
//     std::cout<<" DetId: " <<badchannels_.detid[i]
// 	     <<" fed: "<<badchannels_.fed[i]
// 	     <<" link: "<<badchannels_.link[i]
// 	     <<" roc_first: "<<badchannels_.roc_first[i]
// 	     <<" roc_last: "<<badchannels_.roc_last[i]<<std::endl;
//   }

  eventTree_ -> Fill();


  if (badchannelsLSHists_.find(event_.event)!=badchannelsLSHists_.end()) {
    std::cout<<" *** ERROR: this event has already been processed.\n";
  }

  badchannelsLSHists_[event_.event]=std::pair<TH1I,TH1I>(channels, indices);
  badchannelsLS_[event_.event]=badchannels_;

}

DEFINE_FWK_MODULE(StuckTBMAnalyzer);
