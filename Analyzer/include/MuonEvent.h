#ifndef MuonEvent_h
#define MuonEvent_h
//FWK

class Muon {
    public:
        Int_t numberOfSeeds;
        Double_t pt;
        Double_t eta;
        Double_t phi;
        Double_t dxy;
        Double_t d0;
    
        explicit Muon() {};
        ~Muon();
};

Muon::~Muon(){}


class MuonEvent {
    public:
        Int_t runNumber;
        Int_t luminosityBlockNumber;
        Int_t eventNumber;
        Muon muon;
        
        
        explicit MuonEvent() {};
        ~MuonEvent();
        
};

MuonEvent::~MuonEvent(){}
#endif
