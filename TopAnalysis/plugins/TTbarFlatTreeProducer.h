#ifndef TTbarFlatTreeProducer_h
#define TTbarFlatTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"

class TTbarFlatTreeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit TTbarFlatTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~TTbarFlatTreeProducer();

  private:  
    virtual bool isGoodMuon(const pat::Muon &mu,const reco::Vertex &vtx,float rho);
    virtual bool isGoodElectron(const pat::Electron &el,const reco::Vertex &vtx,float rho);
    virtual bool isGoodJet(const pat::Jet &jet);
    float MuonRelIso(const reco::Candidate *cand,float rho);
    float ElectronRelIso(const reco::Candidate *cand,float rho);
    float LeptonRelIso(const reco::Candidate *cand,float rho){return cand->isElectron() ? ElectronRelIso(cand,rho) : MuonRelIso(cand,rho);}
    void initialize();
    void computeEventShapes(std::vector<const reco::Candidate *> myObj);
    int  findMatch(LorentzVector p4,std::vector<TLorentzVector> vP4);
    //---- configurable parameters --------   
    edm::InputTag srcJets_,srcMET_,srcMuons_,srcElectrons_,srcGenParticles_,srcRho_,srcVtx_,srcQGL_,triggerResults_,triggerPrescales_;
    std::string srcBtag_,srcPU_,kinfit_;
    int    nJetsMin_;
    int    nBJetsMin_;
    double ptMin_;
    double htMin_;
    double etaMax_;
    double btagMinThreshold_,btagMaxThreshold_;
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    TH1F *puHisto_,*cutFlowHisto_;
    //---- TRIGGER -------------------------
    std::vector<std::string> triggerNames_;
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int   run_,evt_,nVtx_,lumi_,nJets_,nBJets_,nLeptons_,status_;
    float rho_,met_,metSig_,ht_,htBtag_,prob_,chi2_,qglAve_,qglMin_,qglMedian_;
    std::vector<bool> *triggerBit_;
    std::vector<int>  *triggerPre_;
    //---- event-shape variables ------
    float sphericity_,aplanarity_,foxWolfram_[4];
    //---- top variables --------------
    int   idxQ_,idxQbar_,idxB_,idxP_,idxPbar_,idxBbar_;
    float mTop_[2],mW_[2],ptTop_[2],yTop_[2],dRbbTop_,mTTbar_,yTTbar_,ptTTbar_;
    float mTopReco_[2],mWReco_[2],dRbbTopReco_;
    //---- jet variables --------------
    std::vector<bool>  *isBtag_;
    std::vector<int>   *flavor_;
    std::vector<float> *pt_,*eta_,*phi_,*mass_,*energy_,*chf_,*nhf_,*phf_,*elf_,*muf_,*btag_,*qgl_,*puMva_;
    //---- lepton variables -----------
    std::vector<int>   *lId_;
    std::vector<float> *lPt_,*lEta_,*lPhi_,*lE_,*lIso_;
    //---- MC variables ---------------
    int npu_,decay_;
};





#endif
