/** \class MuonSeedsAnalyzer
 */
      
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include <map>
#include <string>
#include <iomanip>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TMath.h"
#include <memory>
#include "DisplacedMuonReco/Analyzer/include/MuonEvent.h"


class MuonSeedsAnalyzer : public edm::EDAnalyzer {

 public:
    explicit MuonSeedsAnalyzer(const edm::ParameterSet& cfg);
    ~MuonSeedsAnalyzer();


    //~ virtual void beginRun(const edm::Run & run,    const edm::EventSetup & eventSetup);
    //~ virtual void endRun  (const edm::Run & run,    const edm::EventSetup & eventSetup);

 private:

    virtual void analyze (const edm::Event& event, const edm::EventSetup & eventSetup) override;
    virtual void beginJob() override;
    virtual void endJob() override;
    virtual void beginEvent();
//   virtual void endEvent();

    void fillMuonSeeds(const edm::Handle<std::vector<int>>  &,
                        const edm::Event    &);
    void fillMuonSeedsDisplaced(const edm::Handle<std::vector<int>>  &,
                        const edm::Event    &);
                        
    bool FoundGENtoRECOmatch(const edm::Handle<reco::MuonCollection>      &, 
                        const std::vector<reco::GenParticle>::const_iterator    &,
                        const double             &  );
    bool FoundGENtoRECOmatch(const edm::Handle<reco::TrackCollection>      & ,
                        const std::vector<reco::GenParticle>::const_iterator & ,
                        const double                                  &  );
    void printInfo(reco::TrackCollection::const_iterator & );
    void printGENInfo(reco::GenParticleCollection::const_iterator & );
    
    void fillEff();

    edm::InputTag muonSeededTracksOutInDisplacedTag_;
    edm::EDGetTokenT<reco::TrackCollection>  muonSeededTracksOutInDisplacedToken_;

    edm::InputTag muonSeededTrackCandidatesOutInDisplacedTag_;
    edm::EDGetTokenT<TrackCandidateCollection>  muonSeededTrackCandidatesOutInDisplacedToken_;

    //~ edm::InputTag muonSeededTracksOutInDisplacedClassifierTag_;
    //~ edm::EDGetTokenT<reco::TrackCollection>  muonSeededTracksOutInDisplacedClassifierToken_;

    //~ edm::InputTag displacedtracksTag_;
    //~ edm::EDGetTokenT<reco::TrackCollection>  displacedtracksToken_;

    edm::InputTag earlyDisplacedMuonsTag_;
    edm::EDGetTokenT<std::vector<reco::Muon>>  earlyDisplacedMuonsToken_;
    
    edm::InputTag genTag_;
    edm::EDGetTokenT<reco::GenParticleCollection> genToken_;

    edm::InputTag tkptTag_;
    edm::EDGetTokenT<std::vector<double>>  tkptToken_;

    edm::InputTag tketaTag_;
    edm::EDGetTokenT<std::vector<double>>  tketaToken_;

    edm::InputTag tkphiTag_;
    edm::EDGetTokenT<std::vector<double>>  tkphiToken_;

    edm::InputTag tkd0Tag_;
    edm::EDGetTokenT<std::vector<double>>  tkd0Token_;

    edm::InputTag tkdxyTag_;
    edm::EDGetTokenT<std::vector<double>>  tkdxyToken_;

    edm::InputTag tkDisplacedptTag_;
    edm::EDGetTokenT<std::vector<double>>  tkDisplacedptToken_;

    edm::InputTag tkDisplacedetaTag_;
    edm::EDGetTokenT<std::vector<double>>  tkDisplacedetaToken_;

    edm::InputTag tkDisplacedphiTag_;
    edm::EDGetTokenT<std::vector<double>>  tkDisplacedphiToken_;

    edm::InputTag tkDisplacedd0Tag_;
    edm::EDGetTokenT<std::vector<double>>  tkDisplacedd0Token_;

    edm::InputTag tkDisplaceddxyTag_;
    edm::EDGetTokenT<std::vector<double>>  tkDisplaceddxyToken_;


    edm::InputTag MuonSeededSeedsOutInDisplacedNumSeedsTag_;
    edm::EDGetTokenT<std::vector<int>>  MuonSeededSeedsOutInDisplacedNumSeedsToken_;    

    edm::InputTag MuonSeededSeedsOutInNumSeedsTag_;
    edm::EDGetTokenT<std::vector<int>>  MuonSeededSeedsOutInNumSeedsToken_;    

    edm::Service<TFileService> outfile_;

    double deltaR_GENToRECO;

    MuonEvent event_;
    std::map<std::string,TTree*> tree_;
    
    TH1D* h1_NumSeedsDisplaced;
    TH1D* h1_NumSeeds;
    
    TH1D* h1_pt;
    TH1D* h1_phi;
    TH1D* h1_eta;
    TH1D* h1_d0;
    TH1D* h1_dxy;

    TH1D* h1_GenInfo_num;
    TH1D* h1_GenInfo_pt;
    TH1D* h1_GenInfo_eta;
    TH1D* h1_GenInfo_phi;

    TH1D* h1_Displaced_pt;
    TH1D* h1_Displaced_phi;
    TH1D* h1_Displaced_eta;
    TH1D* h1_Displaced_d0;
    TH1D* h1_Displaced_dxy;

    TH2D* h2_NumSeeds_vs_abseta;
    TH2D* h2_NumSeeds_vs_eta;
    TH2D* h2_NumSeeds_vs_pt;
    TH2D* h2_NumSeeds_vs_phi;
    TH2D* h2_NumSeeds_vs_dxy;
    TH2D* h2_NumSeeds_vs_d0;

    TH2D* h2_Displaced_NumSeeds_vs_abseta;
    TH2D* h2_Displaced_NumSeeds_vs_eta;
    TH2D* h2_Displaced_NumSeeds_vs_pt;
    TH2D* h2_Displaced_NumSeeds_vs_phi;
    TH2D* h2_Displaced_NumSeeds_vs_dxy;
    TH2D* h2_Displaced_NumSeeds_vs_d0;
    
    TH1D* h1_MuonSeededTracksOutInDisplaced_eta;
    TH1D* h1_MuonSeededTracksOutInDisplaced_pt;
    TH1D* h1_MuonSeededTracksOutInDisplaced_phi;
    TH1D* h1_MuonSeededTracksOutInDisplaced_num;
    TH1D* h1_MuonSeededTracksOutInDisplaced_HitPattern_numberOfValidStripLayersWithMonoAndStereo;
    TH1D* h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement;
    
    TH1D* h1_EarlyDisplacedMuons_eta;
    TH1D* h1_EarlyDisplacedMuons_pt;
    TH1D* h1_EarlyDisplacedMuons_phi;
    TH1D* h1_EarlyDisplacedMuons_num;

    //~ TH1D* h1_MuonSeededTracksCandidatesOutInDisplaced_eta;
    //~ TH1D* h1_MuonSeededTracksCandidatesOutInDisplaced_pt;
    //~ TH1D* h1_MuonSeededTracksCandidatesOutInDisplaced_phi;
    
    TH1D* h1_MuonSeededTracksCandidatesOutInDisplaced_num;
    TH1D* h1_MuonSeededTracksCandidatesOutInDisplaced_size;
    
    TEfficiency* eff1_EarlyDisplacedMuons_vs_dxy;
    TEfficiency* eff1_EarlyDisplacedMuons_vs_eta;
    
    TEfficiency* eff1_MuonSeededTracksOutInDisplaced_vs_eta;
    TEfficiency* eff1_MuonSeededTracksOutInDisplaced_vs_dxy;

};

MuonSeedsAnalyzer::~MuonSeedsAnalyzer()
{
    // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


/// default constructor
MuonSeedsAnalyzer::MuonSeedsAnalyzer(const edm::ParameterSet& cfg) 
  //~ offlinePVTag_           (cfg.getParameter<edm::InputTag>("offlineVtx")), 
    //~ offlinePVToken_         (consumes<reco::VertexCollection>(offlinePVTag_)), 
  //~ offlineMuonTag_         (cfg.getParameter<edm::InputTag>("offlineMuons")), 
    //~ offlineMuonToken_       (consumes<std::vector<reco::Muon>>(offlineMuonTag_)), 
//~ 
{

    MuonSeededSeedsOutInDisplacedNumSeedsTag_ = cfg.getParameter<edm::InputTag >("MuonSeededSeedsOutInDisplacedNumSeeds");
    MuonSeededSeedsOutInDisplacedNumSeedsToken_= consumes<std::vector<int>>(MuonSeededSeedsOutInDisplacedNumSeedsTag_);

    MuonSeededSeedsOutInNumSeedsTag_ = cfg.getParameter<edm::InputTag >("MuonSeededSeedsOutInNumSeeds");
    MuonSeededSeedsOutInNumSeedsToken_= consumes<std::vector<int>>(MuonSeededSeedsOutInNumSeedsTag_);

    tkptTag_ = cfg.getParameter<edm::InputTag >("tkpt");
    tkptToken_ = consumes<std::vector<double>>(tkptTag_);

    tketaTag_ = cfg.getParameter<edm::InputTag >("tketa");
    tketaToken_ = consumes<std::vector<double>>(tketaTag_);

    tkphiTag_ = cfg.getParameter<edm::InputTag >("tkphi");
    tkphiToken_ = consumes<std::vector<double>>(tkphiTag_);

    tkd0Tag_ = cfg.getParameter<edm::InputTag >("tkd0");
    tkd0Token_ = consumes<std::vector<double>>(tkd0Tag_);

    tkdxyTag_ = cfg.getParameter<edm::InputTag >("tkdxy");
    tkdxyToken_ = consumes<std::vector<double>>(tkdxyTag_);

    tkDisplacedptTag_ = cfg.getParameter<edm::InputTag >("tkDisplacedpt");
    tkDisplacedptToken_ = consumes<std::vector<double>>(tkDisplacedptTag_);

    tkDisplacedetaTag_ = cfg.getParameter<edm::InputTag >("tkDisplacedeta");
    tkDisplacedetaToken_ = consumes<std::vector<double>>(tkDisplacedetaTag_);

    tkDisplacedphiTag_ = cfg.getParameter<edm::InputTag >("tkDisplacedphi");
    tkDisplacedphiToken_ = consumes<std::vector<double>>(tkDisplacedphiTag_);

    tkDisplacedd0Tag_ = cfg.getParameter<edm::InputTag >("tkDisplacedd0");
    tkDisplacedd0Token_ = consumes<std::vector<double>>(tkDisplacedd0Tag_);

    tkDisplaceddxyTag_ = cfg.getParameter<edm::InputTag >("tkDisplaceddxy");
    tkDisplaceddxyToken_ = consumes<std::vector<double>>(tkDisplaceddxyTag_);
    
    earlyDisplacedMuonsTag_ = cfg.getParameter<edm::InputTag>("earlyDisplacedMuons");
    earlyDisplacedMuonsToken_ = consumes<std::vector<reco::Muon>>(earlyDisplacedMuonsTag_);

    genTag_ = cfg.getParameter<edm::InputTag>("genParticles");
    genToken_ = consumes<reco::GenParticleCollection>(genTag_);

    muonSeededTracksOutInDisplacedTag_ = cfg.getParameter<edm::InputTag>("MuonSeededTracksOutInDisplaced");
    muonSeededTracksOutInDisplacedToken_ = consumes<reco::TrackCollection>(muonSeededTracksOutInDisplacedTag_);

    muonSeededTrackCandidatesOutInDisplacedTag_ = cfg.getParameter<edm::InputTag>("MuonSeededTrackCandidatesOutInDisplaced");
    muonSeededTrackCandidatesOutInDisplacedToken_ = consumes<TrackCandidateCollection>(muonSeededTrackCandidatesOutInDisplacedTag_);

    //~ muonSeededTracksOutInDisplacedClassifierTag_ = cfg.getParameter<edm::InputTag>("MuonSeededTracksOutInDisplacedClassifier");
    //~ muonSeededTracksOutInDisplacedClassifierToken_ = consumes<reco::TrackCollection>(muonSeededTracksOutInDisplacedClassifierTag_);

    
    deltaR_GENToRECO = 0.1;
}


void MuonSeedsAnalyzer::beginJob() {

    TH1::SetDefaultSumw2() ;
    //~ tree_["muonTree"] = outfile_-> make<TTree>("muonTree","muonTree");
    //~ tree_["muonTree"] -> Branch("event.runNumber" ,&event_.runNumber, 64000,2);
    
    TFileDirectory DisplacedReco = outfile_->mkdir("DisplacedReco");
    TFileDirectory StandardReco = outfile_->mkdir("StandardReco");
    TFileDirectory TracksOutInDisplaced = outfile_->mkdir("TracksOutInDisplaced");
    TFileDirectory Effs = outfile_->mkdir("Efficiencies");
    TFileDirectory GenInfo = outfile_->mkdir("GENInfo");
    
    
    h1_GenInfo_num = GenInfo.make<TH1D>("h1_GenInfo_num","h1_GenInfo_num",2,0.,2.);
    h1_GenInfo_pt = GenInfo.make<TH1D>("h1_GenInfo_pt","h1_GenInfo_pt",1000,0.,300.);
    h1_GenInfo_eta = GenInfo.make<TH1D>("h1_GenInfo_eta","h1_GenInfo_eta",100, -4.,4.);
    h1_GenInfo_phi = GenInfo.make<TH1D>("h1_GenInfo_phi","h1_GenInfo_phi",100,-3.14,3.14);
    
    h1_NumSeedsDisplaced = DisplacedReco.make<TH1D>("h1_NumSeedsDisplaced", "h1_NumSeedsDisplaced", 100, 0, 100);
    h1_NumSeeds = StandardReco.make<TH1D>("h1_NumSeeds", "h1_NumSeeds", 100, 0, 100);
    
    h1_pt = StandardReco.make<TH1D>("h1_pt", "h1_pt", 1000, 0, 500);
    h1_eta = StandardReco.make<TH1D>("h1_eta", "h1_eta", 100, -4.,4.);
    h1_phi = StandardReco.make<TH1D>("h1_phi", "h1_phi", 100,-3.14,3.14);
    h1_d0 = StandardReco.make<TH1D>("h1_d0", "h1_d0", 100, 0, 100);
    h1_dxy = StandardReco.make<TH1D>("h1_dxy", "h1_dxy", 100, 0, 100);

    h1_Displaced_pt = DisplacedReco.make<TH1D>("h1_Displaced_pt", "h1_Displaced_pt", 1000, 0, 500);
    h1_Displaced_eta = DisplacedReco.make<TH1D>("h1_Displaced_eta", "h1_Displaced_eta", 100, -4.,4.);
    h1_Displaced_phi = DisplacedReco.make<TH1D>("h1_Displaced_phi", "h1_Displaced_phi", 100,-3.14,3.14);
    h1_Displaced_d0 = DisplacedReco.make<TH1D>("h1_Displaced_d0", "h1_Displaced_d0", 100, 0, 100);
    h1_Displaced_dxy = DisplacedReco.make<TH1D>("h1_Displaced_dxy", "h1_Displaced_dxy", 100, 0, 100);

    h2_NumSeeds_vs_abseta = StandardReco.make<TH2D>("h2_NumSeeds_vs_abseta", "h2_NumSeeds_vs_abseta",100,0.,4., 100, 0, 100);
    h2_NumSeeds_vs_eta = StandardReco.make<TH2D>("h2_NumSeeds_vs_eta", "h2_NumSeeds_vs_eta",100,-4.,4., 100, 0, 100);
    h2_NumSeeds_vs_pt = StandardReco.make<TH2D>("h2_NumSeeds_vs_pt", "h2_NumSeeds_vs_pt",1000,0,1000., 100, 0, 100);
    h2_NumSeeds_vs_phi = StandardReco.make<TH2D>("h2_NumSeeds_vs_phi", "h2_NumSeeds_vs_phi",100,-3.14,3.14, 100, 0, 100);
    h2_NumSeeds_vs_dxy = StandardReco.make<TH2D>("h2_NumSeeds_vs_dxy", "h2_NumSeeds_vs_dxy",100,0,1000, 100, 0, 100);
    h2_NumSeeds_vs_d0 = StandardReco.make<TH2D>("h2_NumSeeds_vs_d0", "h2_NumSeeds_vs_d0",100,0,1000, 100, 0, 100);

    h2_Displaced_NumSeeds_vs_abseta = DisplacedReco.make<TH2D>("h2_Displaced_NumSeeds_vs_abseta", "h2_Displaced_NumSeeds_vs_abseta",100,0.,4., 100, 0, 100);
    h2_Displaced_NumSeeds_vs_eta = DisplacedReco.make<TH2D>("h2_Displaced_NumSeeds_vs_eta", "h2_Displaced_NumSeeds_vs_eta",100,-4.,4., 100, 0, 100);
    h2_Displaced_NumSeeds_vs_pt = DisplacedReco.make<TH2D>("h2_Displaced_NumSeeds_vs_pt", "h2_Displaced_NumSeeds_vs_pt",1000,0,1000., 100, 0, 100);
    h2_Displaced_NumSeeds_vs_phi = DisplacedReco.make<TH2D>("h2_Displaced_NumSeeds_vs_phi", "h2_Displaced_NumSeeds_vs_phi",100,-3.14,3.14, 100, 0, 100);
    h2_Displaced_NumSeeds_vs_dxy = DisplacedReco.make<TH2D>("h2_Displaced_NumSeeds_vs_dxy", "h2_Displaced_NumSeeds_vs_dxy",100,0,1000, 100, 0, 100);
    h2_Displaced_NumSeeds_vs_d0 = DisplacedReco.make<TH2D>("h2_Displaced_NumSeeds_vs_d0", "h2_Displaced_NumSeeds_vs_d0",100,0,1000, 100, 0, 100);
    
    eff1_EarlyDisplacedMuons_vs_eta =  Effs.make<TEfficiency>("eff1_EarlyDisplacedMuons_vs_eta","eff1_EarlyDisplacedMuons_vs_eta; #eta; #epsilon",100,-4.,4.);
    eff1_EarlyDisplacedMuons_vs_dxy =  Effs.make<TEfficiency>("eff1_EarlyDisplacedMuons_vs_dxy","eff1_EarlyDisplacedMuons_vs_dxy; d_{xy} (cm); #epsilon",100,0,50.);

    eff1_MuonSeededTracksOutInDisplaced_vs_eta =  Effs.make<TEfficiency>("eff1_MuonSeededTracksOutInDisplaced_vs_eta","eff1_MuonSeededTracksOutInDisplaced_vs_eta; #eta; #epsilon",100,-4.,4.);
    eff1_MuonSeededTracksOutInDisplaced_vs_dxy =  Effs.make<TEfficiency>("eff1_MuonSeededTracksOutInDisplaced_vs_dxy","eff1_MuonSeededTracksOutInDisplaced_vs_dxy; d_{xy} (cm); #epsilon",100,0,50.);

    h1_MuonSeededTracksOutInDisplaced_pt = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksOutInDisplaced_pt","h1_MuonSeededTracksOutInDisplaced_pt",1000,0,500);
    h1_MuonSeededTracksOutInDisplaced_eta = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksOutInDisplaced_eta","h1_MuonSeededTracksOutInDisplaced_eta",100, -4.,4.);
    h1_MuonSeededTracksOutInDisplaced_phi = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksOutInDisplaced_phi","h1_MuonSeededTracksOutInDisplaced_phi",100,-3.14,3.14);
    h1_MuonSeededTracksOutInDisplaced_num = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksOutInDisplaced_num","h1_MuonSeededTracksOutInDisplaced_num",2,0.,2.);
    h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement","h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement",50,0.,50.);
    h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement","h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement",50,0.,50.);
    
    
    h1_EarlyDisplacedMuons_pt = TracksOutInDisplaced.make<TH1D>("h1_EarlyDisplacedMuons_pt","h1_EarlyDisplacedMuons_pt",1000,0,500);
    h1_EarlyDisplacedMuons_eta = TracksOutInDisplaced.make<TH1D>("h1_EarlyDisplacedMuons_eta","h1_EarlyDisplacedMuons_eta",100, -4.,4.);
    h1_EarlyDisplacedMuons_phi = TracksOutInDisplaced.make<TH1D>("h1_EarlyDisplacedMuons_phi","h1_EarlyDisplacedMuons_phi",100,-3.14,3.14);
    h1_EarlyDisplacedMuons_num = TracksOutInDisplaced.make<TH1D>("h1_EarlyDisplacedMuons_num","h1_EarlyDisplacedMuons_num",2,0.,2.);


    h1_MuonSeededTracksCandidatesOutInDisplaced_num = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksCandidatesOutInDisplaced_num","h1_MuonSeededTracksCandidatesOutInDisplaced_num",2,0.,2.);
    h1_MuonSeededTracksCandidatesOutInDisplaced_size = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksCandidatesOutInDisplaced_size","h1_MuonSeededTracksCandidatesOutInDisplaced_size",10,0.,10.);
    //~ h1_MuonSeededTracksCandidatesOutInDisplaced_pt = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksCandidatesOutInDisplaced_pt","h1_MuonSeededTracksCandidatesOutInDisplaced_pt",1000,0,500);
    //~ h1_MuonSeededTracksCandidatesOutInDisplaced_eta = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksCandidatesOutInDisplaced_eta","h1_MuonSeededTracksCandidatesOutInDisplaced_eta",100, -4.,4.);
    //~ h1_MuonSeededTracksCandidatesOutInDisplaced_phi = TracksOutInDisplaced.make<TH1D>("h1_MuonSeededTracksCandidatesOutInDisplaced_phi","h1_MuonSeededTracksCandidatesOutInDisplaced_phi",100,-3.14,3.14);

}    

void MuonSeedsAnalyzer::endJob() {}

//~ void MuonSeedsAnalyzer::beginRun(const edm::Run & run, const edm::EventSetup & eventSetup) {}
//~ 
//~ void MuonSeedsAnalyzer::endRun  (const edm::Run & run, const edm::EventSetup & eventSetup) {}
 
void MuonSeedsAnalyzer::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) {

    beginEvent();

    // Fill general info
    event_.runNumber             = event.id().run();
    event_.luminosityBlockNumber = event.id().luminosityBlock();
    event_.eventNumber           = event.id().event();
   
    edm::Handle<std::vector<int>> MuonSeededSeedsOutInDisplacedNumSeeds_;
    event.getByToken(MuonSeededSeedsOutInDisplacedNumSeedsToken_,MuonSeededSeedsOutInDisplacedNumSeeds_);
    
    edm::Handle<std::vector<int>> MuonSeededSeedsOutInNumSeeds_;
    event.getByToken(MuonSeededSeedsOutInNumSeedsToken_,MuonSeededSeedsOutInNumSeeds_);
    
    edm::Handle<std::vector<double>> tkpt_;
    event.getByToken(tkptToken_,tkpt_);

    edm::Handle<std::vector<double>> tketa_;
    event.getByToken(tketaToken_,tketa_);
    
    edm::Handle<std::vector<double>> tkphi_;
    event.getByToken(tkphiToken_,tkphi_);

    edm::Handle<std::vector<double>> tkdxy_;
    event.getByToken(tkdxyToken_,tkdxy_);

    edm::Handle<std::vector<double>> tkd0_;
    event.getByToken(tkd0Token_,tkd0_);

    edm::Handle<std::vector<double>> tkDisplacedpt_;
    event.getByToken(tkDisplacedptToken_,tkDisplacedpt_);

    edm::Handle<std::vector<double>> tkDisplacedeta_;
    event.getByToken(tkDisplacedetaToken_,tkDisplacedeta_);
    
    edm::Handle<std::vector<double>> tkDisplacedphi_;
    event.getByToken(tkDisplacedphiToken_,tkDisplacedphi_);

    edm::Handle<std::vector<double>> tkDisplaceddxy_;
    event.getByToken(tkDisplaceddxyToken_,tkDisplaceddxy_);

    edm::Handle<std::vector<double>> tkDisplacedd0_;
    event.getByToken(tkDisplacedd0Token_,tkDisplacedd0_);
    
    edm::Handle<std::vector<reco::Muon>> earlyDisplacedMuons_;
    event.getByToken(earlyDisplacedMuonsToken_, earlyDisplacedMuons_);

    edm::Handle<reco::GenParticleCollection> genParticles_;
    event.getByToken(genToken_, genParticles_);

    //~ edm::Handle<reco::TrackCollection> displacedtracks_;
    //~ event.getByToken(displacedtracksToken_, displacedtracks_);

    edm::Handle<reco::TrackCollection> muonSeededTracksOutInDisplaced_;
    event.getByToken(muonSeededTracksOutInDisplacedToken_, muonSeededTracksOutInDisplaced_);

    edm::Handle<TrackCandidateCollection> muonSeededTrackCandidatesOutInDisplaced_;
    event.getByToken(muonSeededTrackCandidatesOutInDisplacedToken_, muonSeededTrackCandidatesOutInDisplaced_);

    //~ edm::Handle<reco::TrackCollection> muonSeededTracksOutInDisplacedClassifier_;
    //~ event.getByToken(muonSeededTracksOutInDisplacedClassifierToken_, muonSeededTracksOutInDisplacedClassifier_);
    
    //~ std::cout << "Size of: " << muonSeededTrackCandidatesOutInDisplaced_->size() << std::endl;
    
    for (reco::GenParticleCollection::const_iterator genPart = genParticles_->begin(); genPart != genParticles_->end() ; ++genPart){
        printGENInfo(genPart);
        
        h1_GenInfo_num->Fill(0);
        h1_GenInfo_eta->Fill(genPart->eta());
        h1_GenInfo_phi->Fill(genPart->phi());
        h1_GenInfo_pt->Fill(genPart->pt());
        
        
        double dxy = TMath::Abs( (genPart->p4().y()*genPart->vertex().x() -  genPart->p4().x()*genPart->vertex().y())/genPart->pt() );

        bool MatchedGENtoRECO_EarlyDisplacedMuons = FoundGENtoRECOmatch(earlyDisplacedMuons_,genPart,deltaR_GENToRECO);
        bool MatchedGENtoRECO_MuonSeededTracksOutInDisplaced = FoundGENtoRECOmatch(muonSeededTracksOutInDisplaced_,genPart,deltaR_GENToRECO);

        eff1_EarlyDisplacedMuons_vs_dxy->Fill(MatchedGENtoRECO_EarlyDisplacedMuons,dxy);
        eff1_EarlyDisplacedMuons_vs_eta->Fill(MatchedGENtoRECO_EarlyDisplacedMuons,genPart->eta());

        eff1_MuonSeededTracksOutInDisplaced_vs_dxy->Fill(MatchedGENtoRECO_MuonSeededTracksOutInDisplaced,dxy);
        eff1_MuonSeededTracksOutInDisplaced_vs_eta->Fill(MatchedGENtoRECO_MuonSeededTracksOutInDisplaced,genPart->eta());
    
    }
    
    
    // Loop over earlyDisplacedMuons_
    for (reco::MuonCollection::const_iterator muo = earlyDisplacedMuons_->begin(); muo != earlyDisplacedMuons_->end(); ++muo){
        h1_EarlyDisplacedMuons_pt->Fill(muo->pt());
        h1_EarlyDisplacedMuons_eta->Fill(muo->eta());
        h1_EarlyDisplacedMuons_phi->Fill(muo->phi());
        h1_EarlyDisplacedMuons_num->Fill(0);
    }

    h1_MuonSeededTracksCandidatesOutInDisplaced_size->Fill(muonSeededTrackCandidatesOutInDisplaced_->size());
    // Loop over muonSeededTrackCandidatesOutInDisplaced_
    for (TrackCandidateCollection::const_iterator trkCand = muonSeededTrackCandidatesOutInDisplaced_->begin(); trkCand != muonSeededTrackCandidatesOutInDisplaced_->end(); ++trkCand){
        h1_MuonSeededTracksCandidatesOutInDisplaced_num->Fill(0);
        //~ h1_MuonSeededTracksCandidatesOutInDisplaced_eta->Fill(trkCand->eta())
        //~ h1_MuonSeededTracksCandidatesOutInDisplaced_pt->Fill(trkCand->pt())
        //~ h1_MuonSeededTracksCandidatesOutInDisplaced_phi->Fill(trkCand->phi())
        
    }
    // Loop over muonSeededTracksOutInDisplaced_
    for (reco::TrackCollection::const_iterator trk = muonSeededTracksOutInDisplaced_->begin(); trk != muonSeededTracksOutInDisplaced_->end(); ++trk){
        //~ printInfo(trk);
        h1_MuonSeededTracksOutInDisplaced_pt->Fill(trk->pt());
        h1_MuonSeededTracksOutInDisplaced_eta->Fill(trk->eta());
        h1_MuonSeededTracksOutInDisplaced_phi->Fill(trk->phi());
        h1_MuonSeededTracksOutInDisplaced_num->Fill(0);
        h1_MuonSeededTracksOutInDisplaced_HitPattern_numberOfValidStripLayersWithMonoAndStereo->Fill(trk->hitPattern().numberOfValidStripLayersWithMonoAndStereo());
        //~ h1_MuonSeededTracksOutInDisplaced_HitPattern_stripLayersWithMeasurement->Fill(trk->hitPattern().stripLayersWithMeasurement());
    }
    
    //~ std::vector<reco::Muon>::const_iterator mu1=muons->begin(); mu1!=muons->end(); ++mu1
    if (event.getByToken(MuonSeededSeedsOutInDisplacedNumSeedsToken_,MuonSeededSeedsOutInDisplacedNumSeeds_)){
        fillMuonSeedsDisplaced(MuonSeededSeedsOutInDisplacedNumSeeds_,event);
        for (unsigned int iter = 0; iter < MuonSeededSeedsOutInDisplacedNumSeeds_->size(); iter++){
            h1_Displaced_pt->Fill(tkDisplacedpt_->at(iter));
            h1_Displaced_eta->Fill(tkDisplacedeta_->at(iter));
            h1_Displaced_phi->Fill(tkDisplacedphi_->at(iter));
            h1_Displaced_d0->Fill(tkDisplacedd0_->at(iter));
            h1_Displaced_dxy->Fill(tkDisplaceddxy_->at(iter));
            
            h2_Displaced_NumSeeds_vs_pt->Fill(tkDisplacedpt_->at(iter),MuonSeededSeedsOutInDisplacedNumSeeds_->at(iter));
            h2_Displaced_NumSeeds_vs_eta->Fill(tkDisplacedeta_->at(iter),MuonSeededSeedsOutInDisplacedNumSeeds_->at(iter));
            h2_Displaced_NumSeeds_vs_abseta->Fill(std::abs(tkDisplacedeta_->at(iter)),MuonSeededSeedsOutInDisplacedNumSeeds_->at(iter));
            h2_Displaced_NumSeeds_vs_phi->Fill(tkDisplacedphi_->at(iter),MuonSeededSeedsOutInDisplacedNumSeeds_->at(iter));
            h2_Displaced_NumSeeds_vs_d0->Fill(tkDisplacedd0_->at(iter),MuonSeededSeedsOutInDisplacedNumSeeds_->at(iter));
            h2_Displaced_NumSeeds_vs_dxy->Fill(tkDisplaceddxy_->at(iter),MuonSeededSeedsOutInDisplacedNumSeeds_->at(iter));
        }
    }

    
    if (event.getByToken(MuonSeededSeedsOutInNumSeedsToken_,MuonSeededSeedsOutInNumSeeds_)){
        fillMuonSeeds(MuonSeededSeedsOutInNumSeeds_,event);
        for (unsigned int iter = 0; iter < MuonSeededSeedsOutInNumSeeds_->size(); iter++){
            h1_pt->Fill(tkpt_->at(iter));
            h1_eta->Fill(tketa_->at(iter));
            h1_phi->Fill(tkphi_->at(iter));
            h1_d0->Fill(tkd0_->at(iter));
            h1_dxy->Fill(tkdxy_->at(iter));
            
            
            h2_NumSeeds_vs_pt->Fill(tkpt_->at(iter),MuonSeededSeedsOutInNumSeeds_->at(iter));
            h2_NumSeeds_vs_eta->Fill(tketa_->at(iter),MuonSeededSeedsOutInNumSeeds_->at(iter));
            h2_NumSeeds_vs_abseta->Fill(std::abs(tketa_->at(iter)),MuonSeededSeedsOutInNumSeeds_->at(iter));
            h2_NumSeeds_vs_phi->Fill(tkphi_->at(iter),MuonSeededSeedsOutInNumSeeds_->at(iter));
            h2_NumSeeds_vs_d0->Fill(tkd0_->at(iter),MuonSeededSeedsOutInNumSeeds_->at(iter));
            h2_NumSeeds_vs_dxy->Fill(tkdxy_->at(iter),MuonSeededSeedsOutInNumSeeds_->at(iter));
        }
    }

    // endEvent();
    //~ tree_["muonTree"] -> Fill();
}

void MuonSeedsAnalyzer::fillMuonSeedsDisplaced(const edm::Handle<std::vector<int>>   &tnumseeds,
                                               const edm::Event                      &tevent
                                                )
{
    for (std::vector<int>::const_iterator num = tnumseeds->begin(); num != tnumseeds->end(); num++){
        h1_NumSeedsDisplaced->Fill(*num);
    }
}

void MuonSeedsAnalyzer::fillMuonSeeds(const edm::Handle<std::vector<int>>   &tnumseeds,
                                      const edm::Event                      &tevent
                                        )
{
    for (std::vector<int>::const_iterator num = tnumseeds->begin(); num != tnumseeds->end(); num++){
        h1_NumSeeds->Fill(*num);
    }
    
}

void MuonSeedsAnalyzer::fillEff(){
    
}

void MuonSeedsAnalyzer::printInfo(reco::TrackCollection::const_iterator & trk){
    std::cout << "Info about track: " << std::endl;
    std::cout << "pt = " << trk->pt() << " GeV" << std::endl;
    std::cout << "eta = " << trk->eta() << std::endl;
    std::cout << "phi = " << trk->phi() << std::endl;
    
}

void MuonSeedsAnalyzer::printGENInfo(reco::GenParticleCollection::const_iterator & genPart){
    std::cout << "Info about GEN track: " << std::endl;
    std::cout << "pt_GEN = " << genPart->pt() << " GeV" << std::endl;
    std::cout << "eta_GEN = " << genPart->eta() << std::endl;
    std::cout << "phi_GEN = " << genPart->phi() << std::endl;
    
}

bool MuonSeedsAnalyzer::FoundGENtoRECOmatch(const edm::Handle<reco::MuonCollection>       & muons ,
                                            const std::vector<reco::GenParticle>::const_iterator & genPart,
                                            const double                                  & deltaR ){
    
    bool FoundMatch = false;
        for (std::vector<reco::Muon>::const_iterator muo = muons->begin(); muo != muons->end(); ++muo){
            if (TMath::Sqrt((genPart->phi() - muo->phi())*(genPart->phi() - muo->phi()) + (genPart->eta() - muo->eta())*(genPart->eta() - muo->eta())) < deltaR){
                FoundMatch = true;
            }
            if (FoundMatch) break;
        }
    return FoundMatch;
    
}

bool MuonSeedsAnalyzer::FoundGENtoRECOmatch(const edm::Handle<reco::TrackCollection>      & tracks ,
                                            const std::vector<reco::GenParticle>::const_iterator & genPart,
                                            const double                                  & deltaR ){
    
    bool FoundMatch = false;
        for (reco::TrackCollection::const_iterator trk = tracks->begin(); trk != tracks->end(); ++trk){
            if (TMath::Sqrt((genPart->phi() - trk->phi())*(genPart->phi() - trk->phi()) + (genPart->eta() - trk->eta())*(genPart->eta() - trk->eta())) < deltaR){
                FoundMatch = true;
            }
            if (FoundMatch) break;
        }
    return FoundMatch;
    
}

//---------------------------------------------
void MuonSeedsAnalyzer::beginEvent()
{

  //~ event_.hlt.triggers.clear();
}


// define this as a plug-in
DEFINE_FWK_MODULE(MuonSeedsAnalyzer);
