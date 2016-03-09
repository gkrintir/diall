#ifndef lwMuon_h
#define lwMuon_h

//
// light-weight muon candidate
//

#include "UserCode/diall/interface/particleBase.h"

class lwMuon : public particleBase {
 public:
  lwMuon();
  lwMuon(Double_t pt, Double_t eta, Double_t phi, Double_t m, Int_t id = -1);
  virtual ~lwMuon() {;}
  lwMuon(const lwMuon& obj); // copy constructor
  lwMuon& operator=(const lwMuon& other); // assignment
  
  void     SetClosestGen(Int_t i)      { fMatchId1 = i; }
  void     SetClosestPF(Int_t i)       { fMatchId2  = i; }
  
  void     SetTrkChi2(float val) {fTrkChi2 = val; }
  void     SetGlbChi2(float val) {fGlbChi2 = val; }
  void     SetNMuHits(float val) {fNMuHits = val; }
  void     SetMatchedStations(float val) {fMS = val; }
  void     SetDxy(float val) {fDxy = val; }
  void     SetDz(float val) {fDz = val; }
  void     SetTrkDxy(float val) {fTrkDxy = val; }
  void     SetTrkDz(float val) {fTrkDz = val; }
  void     SetNPixHits(float val) {fNPixHits = val; }
  void     SetTrkLWM(float val) {fTrkLWM = val; } 

  Int_t    GetClosestGen()   const     { return fMatchId1; }
  Int_t    GetClosestPF()    const     { return fMatchId2; }
  
  Float_t   GetTrkChi2() {return fTrkChi2 ; }
  Float_t   GetGlbChi2() {return fGlbChi2 ; }
  Int_t     GetNMuHits() {return fNMuHits ; }
  Int_t     GetMatchedStations() {return fMS ; }
  Float_t   GetDxy() {return fDxy ; }
  Float_t   GetDz() {return fDz ; }
  Float_t   GetTrkDxy() {return fTrkDxy ; }
  Float_t   GetTrkDz() {return fTrkDz ; }
  Int_t     GetNPixHits() {return fNPixHits ; }
  Int_t     GetTrkLWM() {return fTrkLWM ; }
  
 protected:
  
  Float_t                      fTrkChi2;     //  chi2
  Float_t                      fGlbChi2;     //  chi2
  Int_t                        fNMuHits;     //  muon hits
  Int_t                        fMS;          //  #matched stations
  Float_t                      fDxy;         //  dxy
  Float_t                      fDz;          //  dz
  Float_t                      fTrkDxy;      //  innerTrack dxy
  Float_t                      fTrkDz;       //  innerTrack dz
  Int_t                        fNPixHits;    //  pixel hits
  Int_t                        fTrkLWM;      //  tracker layer hits

  ClassDef(lwMuon,1)
};
#endif
