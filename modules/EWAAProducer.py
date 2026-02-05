import ROOT
from ROOT import TLorentzVector
from itertools import combinations, permutations
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os,copy
import numpy as np
from numpy import sign
from numpy import argsort

def find_two_smallest(goodids, fakeids):
  allids=goodids+fakeids
  allids.sort()
  p1=allids[0]
  p2=allids[1]
  return (p1,p2)

class EWAAProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("met_user","F")
    self.out.branch("met_phi_user","F")
    self.out.branch("LooseElectron_id","I",lenVar="nLooseElectron")
    self.out.branch("LooseMuon_id","I",lenVar="nLooseMuon")
    self.out.branch("GoodPhoton_id","I",lenVar="nGoodPhoton")
    self.out.branch("FakePhoton_id","I",lenVar="nFakePhoton")
    self.out.branch("TightJet_id","I",lenVar="nTightJet")
    self.out.branch("TightJet_pt","F",lenVar="nTightJet")
    self.out.branch("TightJet_eta","F",lenVar="nTightJet")
    self.out.branch("TightJet_phi","F",lenVar="nTightJet")
    self.out.branch("TightJet_mass","F",lenVar="nTightJet")
    self.out.branch("j1_pt","F")
    self.out.branch("j1_eta","F")
    self.out.branch("j1_phi","F")
    self.out.branch("j1_mass","F")
    self.out.branch("j1_id","I")
    self.out.branch("j2_pt","F")
    self.out.branch("j2_eta","F")
    self.out.branch("j2_phi","F")
    self.out.branch("j2_mass","F")
    self.out.branch("j2_id","I")
    self.out.branch("j3_pt","F")
    self.out.branch("j3_eta","F")
    self.out.branch("j3_phi","F")
    self.out.branch("j3_mass","F")
    self.out.branch("j3_id","I")
    self.out.branch("dRjj","F")
    self.out.branch("dEtajj","F")
    self.out.branch("dPhijj","F")
    self.out.branch("mjj","F")
    self.out.branch("fake_flag","I")
    self.out.branch("pho1_pt_SR","F")
    self.out.branch("pho1_eta_SR","F")
    self.out.branch("pho1_phi_SR","F")
    self.out.branch("pho1_id_SR","I")
    self.out.branch("pho2_pt_SR","F")
    self.out.branch("pho2_eta_SR","F")
    self.out.branch("pho2_phi_SR","F")
    self.out.branch("pho2_id_SR","I")
    self.out.branch("dR_p1p2_SR","F")
    self.out.branch("dPhi_p1p2_SR","F")
    self.out.branch("dEta_p1p2_SR","F")
    self.out.branch("Maa","F")
    self.out.branch("Ptaa","F")
    self.out.branch("Etaaa","F")
    self.out.branch("Phiaa","F")
    self.out.branch("dR_p1j1_SR","F")
    self.out.branch("dR_p1j2_SR","F")
    self.out.branch("dR_p2j1_SR","F")
    self.out.branch("dR_p2j2_SR","F")
    self.out.branch("dPhi_p1j1_SR","F")
    self.out.branch("dPhi_p1j2_SR","F")
    self.out.branch("dPhi_p2j1_SR","F")
    self.out.branch("dPhi_p2j2_SR","F")
    self.out.branch("zepp_SR","F")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    if (event.PV_npvsGood<1): return False

    met_user=-99
    met_phi_user=-99
    if self.is_mc:
      met_user=event.MET_T1Smear_pt
      met_phi_user=event.MET_T1Smear_phi
    else:
      met_user=event.MET_T1_pt
      met_phi_user=event.MET_T1_phi

    self.out.fillBranch("met_user",met_user)
    self.out.fillBranch("met_phi_user",met_phi_user)

    # recover HLT
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
        for iobj in range(0,event.nTrigObj):
          if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
            HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    GoodPhoton_id = []
    FakePhoton_id = []

    photons = Collection(event, 'Photon')
    for ipho in range(0, event.nPhoton):
      if photons[ipho].pt<25:continue
      if abs(photons[ipho].eta)>1.4442 and abs(photons[ipho].eta)<1.57:continue
      if photons[ipho].pixelSeed:continue
      # barrel: -0.02, endcap: -0.26
      if photons[ipho].mvaID_WP90:
        GoodPhoton_id.append(ipho)

      elif (photons[ipho].mvaID>-0.9 and photons[ipho].mvaID<-0.7):
        FakePhoton_id.append(ipho)

    LooseElectron_id = []
    eles = Collection(event, 'Electron')
    for iele in range(0, event.nElectron):
      if eles[iele].convVeto:continue
      if eles[iele].pt>10 and eles[iele].cutBased>0:
        LooseElectron_id.append(iele)

    LooseMuon_id = []
    muons = Collection(event, 'Muon')
    for imu in range(0, event.nMuon):
      if muons[imu].looseId and muons[imu].pt>10:
        LooseMuon_id.append(imu)

    jets = Collection(event, 'Jet')
    TightJet_id = []
    TightJet_pt = []
    TightJet_eta = []
    TightJet_phi = []
    TightJet_mass = []
    TightJet_v4 = []
    jet_v4_temp=TLorentzVector()
    lep_v4_temp=TLorentzVector()
    pho_v4_temp=TLorentzVector()
    for ijet in range(0, event.nJet):
      if abs(jets[ijet].eta)>4.7 or jets[ijet].pt_nom<30: continue
      if jets[ijet].jetId<6:continue
      jet_v4_temp.SetPtEtaPhiM(jets[ijet].pt_nom,jets[ijet].eta,jets[ijet].phi,jets[ijet].mass_nom)
      pass_mu_dr=1
      pass_ele_dr=1
      pass_jet_dr=1
      pass_pho_dr=1

      for imu in range(0,len(LooseMuon_id)):
        if pass_mu_dr<1:continue
        mid_tmp=LooseMuon_id[imu]
        lep_v4_temp.SetPtEtaPhiM(muons[mid_tmp].pt, muons[mid_tmp].eta, muons[mid_tmp].phi, muons[mid_tmp].mass)
        if jet_v4_temp.DeltaR(lep_v4_temp)<0.4:pass_mu_dr=0
      
      for iele in range(0,len(LooseElectron_id)):
        if pass_ele_dr<1:continue
        eid_tmp=LooseElectron_id[iele]
        lep_v4_temp.SetPtEtaPhiM(eles[eid_tmp].pt, eles[eid_tmp].eta, eles[eid_tmp].phi, eles[eid_tmp].mass)
        if jet_v4_temp.DeltaR(lep_v4_temp)<0.4:pass_ele_dr=0

      if len(TightJet_id)>0:
        for ij in range(0,len(TightJet_id)):
          if jet_v4_temp.DeltaR(TightJet_v4[ij])<0.4:pass_jet_dr=0

      for ipho in range(0,len(GoodPhoton_id)):
        if pass_pho_dr<1:continue
        phoid_tmp=GoodPhoton_id[ipho]
        pho_v4_temp.SetPtEtaPhiM(photons[phoid_tmp].pt, photons[phoid_tmp].eta, photons[phoid_tmp].phi, photons[phoid_tmp].mass)
        if jet_v4_temp.DeltaR(pho_v4_temp)<0.4:pass_pho_dr=0
      
      if pass_mu_dr>0 and pass_ele_dr>0 and pass_jet_dr>0 and pass_pho_dr>0:
        TightJet_id.append(ijet)
        TightJet_v4.append(jet_v4_temp.Clone())
        TightJet_pt.append(jet_v4_temp.Clone().Pt())
        TightJet_eta.append(jet_v4_temp.Clone().Eta())
        TightJet_phi.append(jet_v4_temp.Clone().Phi())
        TightJet_mass.append(jet_v4_temp.Clone().M())

    self.out.fillBranch("LooseElectron_id",LooseElectron_id)
    self.out.fillBranch("LooseMuon_id",LooseMuon_id)
    self.out.fillBranch("GoodPhoton_id",GoodPhoton_id)
    self.out.fillBranch("FakePhoton_id",FakePhoton_id)
    self.out.fillBranch("TightJet_id", TightJet_id)
    self.out.fillBranch("TightJet_pt", TightJet_pt)
    self.out.fillBranch("TightJet_eta", TightJet_eta)
    self.out.fillBranch("TightJet_phi", TightJet_phi)
    self.out.fillBranch("TightJet_mass", TightJet_mass)

    if len(TightJet_id)<2:return False
    if len(GoodPhoton_id)+len(FakePhoton_id)<2:return False

    j1_pt=jets[TightJet_id[0]].pt_nom
    j1_eta=jets[TightJet_id[0]].eta
    j1_phi=jets[TightJet_id[0]].phi
    j1_mass=jets[TightJet_id[0]].mass_nom
    j1_id=TightJet_id[0]
    j2_pt=jets[TightJet_id[1]].pt_nom
    j2_eta=jets[TightJet_id[1]].eta
    j2_phi=jets[TightJet_id[1]].phi
    j2_mass=jets[TightJet_id[1]].mass_nom
    j2_id=TightJet_id[1]
    j1_p4=TLorentzVector()
    j2_p4=TLorentzVector()
    j3_p4=TLorentzVector()
    j1_p4.SetPtEtaPhiM(j1_pt,j1_eta,j1_phi,j1_mass)
    j2_p4.SetPtEtaPhiM(j2_pt,j2_eta,j2_phi,j2_mass)
    j3_pt=-99
    j3_eta=-99
    j3_phi=-99
    j3_mass=-99
    j3_id=-99
    if len(TightJet_id)>2:
      j3_pt=jets[TightJet_id[2]].pt_nom
      j3_eta=jets[TightJet_id[2]].eta
      j3_phi=jets[TightJet_id[2]].phi
      j3_mass=jets[TightJet_id[2]].mass_nom
      j3_id=TightJet_id[2]
      j3_p4.SetPtEtaPhiM(j3_pt,j3_eta,j3_phi,j3_mass)
    dRjj=j1_p4.DeltaR(j2_p4)
    dEtajj=abs(j1_eta - j2_eta)
    dPhijj=j1_p4.DeltaPhi(j2_p4)
    mjj=(j1_p4 + j2_p4).M()
    if mjj<150: return False
    
    self.out.fillBranch("j1_pt", j1_pt)
    self.out.fillBranch("j1_eta", j1_eta)
    self.out.fillBranch("j1_phi", j1_phi)
    self.out.fillBranch("j1_mass", j1_mass)
    self.out.fillBranch("j1_id", j1_id)
    self.out.fillBranch("j2_pt", j2_pt)
    self.out.fillBranch("j2_eta", j2_eta)
    self.out.fillBranch("j2_phi", j2_phi)
    self.out.fillBranch("j2_mass", j2_mass)
    self.out.fillBranch("j2_id", j2_id)
    self.out.fillBranch("j3_pt", j3_pt)
    self.out.fillBranch("j3_eta", j3_eta)
    self.out.fillBranch("j3_phi", j3_phi)
    self.out.fillBranch("j3_mass", j3_mass)
    self.out.fillBranch("j3_id", j3_id)
    self.out.fillBranch("dRjj", dRjj)
    self.out.fillBranch("dEtajj", dEtajj)
    self.out.fillBranch("dPhijj", dPhijj)
    self.out.fillBranch("mjj", mjj)

    # fake_flag, 00: both prompt, 10: leading photon is fake, 01: subleading photon is fake, 11: both are fake
    fake_flag=-99
    pho1_pt_SR=-99
    pho1_eta_SR=-99
    pho1_phi_SR=-99
    pho1_id_SR=-99
    pho2_pt_SR=-99
    pho2_eta_SR=-99
    pho2_phi_SR=-99
    pho2_id_SR=-99
    dR_p1p2_SR=-99
    dPhi_p1p2_SR=-99
    dEta_p1p2_SR=-99
    Maa=-99
    Ptaa=-99
    Etaaa=-99
    Phiaa=-99
    dR_p1j1_SR=-99
    dR_p1j2_SR=-99
    dR_p2j1_SR=-99
    dR_p2j2_SR=-99
    dPhi_p1j1_SR=-99
    dPhi_p1j2_SR=-99
    dPhi_p2j1_SR=-99
    dPhi_p2j2_SR=-99
    zepp_SR=-99

    photon1_p4=TLorentzVector()
    photon2_p4=TLorentzVector()
    # for fake photon estimation

    if len(GoodPhoton_id)+len(FakePhoton_id)==2:
      if len(GoodPhoton_id)==2:
        fake_flag=0
        photon1_p4.SetPtEtaPhiM(photons[GoodPhoton_id[0]].pt,photons[GoodPhoton_id[0]].eta,photons[GoodPhoton_id[0]].phi,0)
        photon2_p4.SetPtEtaPhiM(photons[GoodPhoton_id[1]].pt,photons[GoodPhoton_id[1]].eta,photons[GoodPhoton_id[1]].phi,0)
        pho1_id_SR=GoodPhoton_id[0]
        pho2_id_SR=GoodPhoton_id[1]
      elif len(GoodPhoton_id)==1:
        if photons[GoodPhoton_id[0]].pt<photons[FakePhoton_id[0]].pt:
          fake_flag=2
          photon1_p4.SetPtEtaPhiM(photons[FakePhoton_id[0]].pt,photons[FakePhoton_id[0]].eta,photons[FakePhoton_id[0]].phi,0)
          photon2_p4.SetPtEtaPhiM(photons[GoodPhoton_id[0]].pt,photons[GoodPhoton_id[0]].eta,photons[GoodPhoton_id[0]].phi,0)
          pho1_id_SR=FakePhoton_id[0]
          pho2_id_SR=GoodPhoton_id[0]
        else:
          fake_flag=1
          photon1_p4.SetPtEtaPhiM(photons[GoodPhoton_id[0]].pt,photons[GoodPhoton_id[0]].eta,photons[GoodPhoton_id[0]].phi,0)
          photon2_p4.SetPtEtaPhiM(photons[FakePhoton_id[0]].pt,photons[FakePhoton_id[0]].eta,photons[FakePhoton_id[0]].phi,0)
          pho1_id_SR=GoodPhoton_id[0]
          pho2_id_SR=FakePhoton_id[0]
      else:
        fake_flag=3
        photon1_p4.SetPtEtaPhiM(photons[FakePhoton_id[0]].pt,photons[FakePhoton_id[0]].eta,photons[FakePhoton_id[0]].phi,0)
        photon2_p4.SetPtEtaPhiM(photons[FakePhoton_id[1]].pt,photons[FakePhoton_id[1]].eta,photons[FakePhoton_id[1]].phi,0)
        pho1_id_SR=FakePhoton_id[0]
        pho2_id_SR=FakePhoton_id[1]

      pho1_pt_SR=photon1_p4.Pt()
      pho1_eta_SR=photon1_p4.Eta()
      pho1_phi_SR=photon1_p4.Phi()
      pho2_pt_SR=photon2_p4.Pt()
      pho2_eta_SR=photon2_p4.Eta()
      pho2_phi_SR=photon2_p4.Phi()
      dR_p1p2_SR=photon1_p4.DeltaR(photon2_p4)
      dPhi_p1p2_SR=photon1_p4.DeltaPhi(photon2_p4)
      dEta_p1p2_SR=abs(pho1_eta_SR - pho2_eta_SR)
      Maa=(photon1_p4+photon2_p4).M()
      Ptaa=(photon1_p4+photon2_p4).Pt()
      Etaaa=(photon1_p4+photon2_p4).Eta()
      Phiaa=(photon1_p4+photon2_p4).Phi()
      dR_p1j1_SR=photon1_p4.DeltaR(j1_p4)
      dR_p1j2_SR=photon1_p4.DeltaR(j2_p4)
      dR_p2j1_SR=photon2_p4.DeltaR(j1_p4)
      dR_p2j2_SR=photon2_p4.DeltaR(j2_p4)
      dPhi_p1j1_SR=photon1_p4.DeltaPhi(j1_p4)
      dPhi_p1j2_SR=photon1_p4.DeltaPhi(j2_p4)
      dPhi_p2j1_SR=photon2_p4.DeltaPhi(j1_p4)
      dPhi_p2j2_SR=photon2_p4.DeltaPhi(j2_p4)
      zepp_SR=abs((photon1_p4+photon2_p4).Eta() - 0.5*(j1_eta+j2_eta))

    elif len(GoodPhoton_id)+len(FakePhoton_id)>2:
      p1_id,p2_id=find_two_smallest(GoodPhoton_id, FakePhoton_id)
      photon1_p4.SetPtEtaPhiM(photons[p1_id].pt,photons[p1_id].eta,photons[p1_id].phi,0)
      photon2_p4.SetPtEtaPhiM(photons[p2_id].pt,photons[p2_id].eta,photons[p2_id].phi,0)
      pho1_id_SR=p1_id
      pho2_id_SR=p2_id

      if p1_id in GoodPhoton_id:
        if p2_id in GoodPhoton_id:
          fake_flag=0
        else:
          fake_flag=1
      else:
        if p2_id in GoodPhoton_id:
          fake_flag=2
        else:
          fake_flag=3

      pho1_pt_SR=photon1_p4.Pt()
      pho1_eta_SR=photon1_p4.Eta()
      pho1_phi_SR=photon1_p4.Phi()
      pho2_pt_SR=photon2_p4.Pt()
      pho2_eta_SR=photon2_p4.Eta()
      pho2_phi_SR=photon2_p4.Phi()
      dR_p1p2_SR=photon1_p4.DeltaR(photon2_p4)
      dPhi_p1p2_SR=photon1_p4.DeltaPhi(photon2_p4)
      dEta_p1p2_SR=abs(pho1_eta_SR - pho2_eta_SR)
      Maa=(photon1_p4+photon2_p4).M()
      Ptaa=(photon1_p4+photon2_p4).Pt()
      Etaaa=(photon1_p4+photon2_p4).Eta()
      Phiaa=(photon1_p4+photon2_p4).Phi()
      dR_p1j1_SR=photon1_p4.DeltaR(j1_p4)
      dR_p1j2_SR=photon1_p4.DeltaR(j2_p4)
      dR_p2j1_SR=photon2_p4.DeltaR(j1_p4)
      dR_p2j2_SR=photon2_p4.DeltaR(j2_p4)
      dPhi_p1j1_SR=photon1_p4.DeltaPhi(j1_p4)
      dPhi_p1j2_SR=photon1_p4.DeltaPhi(j2_p4)
      dPhi_p2j1_SR=photon2_p4.DeltaPhi(j1_p4)
      dPhi_p2j2_SR=photon2_p4.DeltaPhi(j2_p4)
      zepp_SR=abs((photon1_p4+photon2_p4).Eta() - 0.5*(j1_eta+j2_eta))

    self.out.fillBranch("fake_flag",fake_flag)
    self.out.fillBranch("pho1_pt_SR",pho1_pt_SR)
    self.out.fillBranch("pho1_eta_SR",pho1_eta_SR)
    self.out.fillBranch("pho1_phi_SR",pho1_phi_SR)
    self.out.fillBranch("pho1_id_SR",pho1_id_SR)
    self.out.fillBranch("pho2_pt_SR",pho2_pt_SR)
    self.out.fillBranch("pho2_eta_SR",pho2_eta_SR)
    self.out.fillBranch("pho2_phi_SR",pho2_phi_SR)
    self.out.fillBranch("pho2_id_SR",pho2_id_SR)
    self.out.fillBranch("dR_p1p2_SR",dR_p1p2_SR)
    self.out.fillBranch("dPhi_p1p2_SR",dPhi_p1p2_SR)
    self.out.fillBranch("dEta_p1p2_SR",dEta_p1p2_SR)
    self.out.fillBranch("Maa",Maa)
    self.out.fillBranch("Ptaa",Ptaa)
    self.out.fillBranch("Etaaa",Etaaa)
    self.out.fillBranch("Phiaa",Phiaa)
    self.out.fillBranch("dR_p1j1_SR",dR_p1j1_SR)
    self.out.fillBranch("dR_p1j2_SR",dR_p1j2_SR)
    self.out.fillBranch("dR_p2j1_SR",dR_p2j1_SR)
    self.out.fillBranch("dR_p2j2_SR",dR_p2j2_SR)
    self.out.fillBranch("dPhi_p1j1_SR",dPhi_p1j1_SR)
    self.out.fillBranch("dPhi_p1j2_SR",dPhi_p1j2_SR)
    self.out.fillBranch("dPhi_p2j1_SR",dPhi_p2j1_SR)
    self.out.fillBranch("dPhi_p2j2_SR",dPhi_p2j2_SR)
    self.out.fillBranch("zepp_SR",zepp_SR)

    return True

EWAA2016apv = lambda: EWAAProducer("2016apv")
EWAA2016 = lambda: EWAAProducer("2016")
EWAA2017 = lambda: EWAAProducer("2017")
EWAA2018 = lambda: EWAAProducer("2018")
