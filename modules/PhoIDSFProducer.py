import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os

class PhoIDSFProducer(Module):
  def __init__( self , year ):
    self.year = year
    self.pixel_veto = "Pho_PV.root"
    self.id_mva = "Pho_MVA90.root"
    self.SF_location_path = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/data/year%s/" %(os.environ['CMSSW_BASE'], self.year)
    print 'SF location:', self.SF_location_path

  def beginJob(self):
    print 'begin to set Photon ID SF --->>>'
    print 'start to open SF root file --->>>'
    # init the TH2F
    self.pixel_veto_mva_dir= ROOT.TDirectoryFile()
    self.pixel_veto_mva_th1f= ROOT.TH1F()
    self.id_mva_th2f= ROOT.TH2F()
    #Open the SF root file
    self.file_pixel_veto= ROOT.TFile.Open(self.SF_location_path+self.pixel_veto)
    self.file_id_mva= ROOT.TFile.Open(self.SF_location_path+self.id_mva)
    #access to the TH2F
    self.file_pixel_veto.GetObject('MVAID', self.pixel_veto_mva_dir)
    self.pixel_veto_mva_dir.GetObject('SF_HasPix_MVAID', self.pixel_veto_mva_th1f)
    self.file_id_mva.GetObject('EGamma_SF2D', self.id_mva_th2f)
    print 'open SF files successfully --->>>'

  def endJob(self):
    print 'close SF root file --->>>'
    self.file_pixel_veto.Close()
    self.file_id_mva.Close()
    print 'finish setting Photon ID SF --->>>'
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch('Photon_MVA90_ID_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_MVA90_ID_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_PixelVeto_MVAID_Inc_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_PixelVeto_MVAID_high_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_PixelVeto_MVAID_low_SF','F', lenVar='nPhoton')
    self.out.branch('Photon_PixelVeto_MVAID_Inc_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_PixelVeto_MVAID_high_SFerr','F', lenVar='nPhoton')
    self.out.branch('Photon_PixelVeto_MVAID_low_SFerr','F', lenVar='nPhoton')
  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    photons = Collection(event, "Photon")
    if not (len(photons)>0): pass
    Photon_MVA90_ID_SF = []
    Photon_MVA90_ID_SFerr = []
    Photon_PixelVeto_MVAID_Inc_SF = []
    Photon_PixelVeto_MVAID_Inc_SFerr = []
    Photon_PixelVeto_MVAID_high_SF = []
    Photon_PixelVeto_MVAID_high_SFerr = []
    Photon_PixelVeto_MVAID_low_SF = []
    Photon_PixelVeto_MVAID_low_SFerr = []
    
    for ipho in range(0, len(photons)):
      if photons[ipho].pt < 500: 
        Photon_MVA90_ID_SF.append(self.id_mva_th2f.GetBinContent(self.id_mva_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        Photon_MVA90_ID_SFerr.append(self.id_mva_th2f.GetBinError(self.id_mva_th2f.FindBin(photons[ipho].eta, photons[ipho].pt)))
        if abs(photons[ipho].eta)<1.566:
          Photon_PixelVeto_MVAID_Inc_SF.append(self.pixel_veto_mva_th1f.GetBinContent(1))
          Photon_PixelVeto_MVAID_Inc_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(1))
          if photons[ipho].r9>0.96:
            Photon_PixelVeto_MVAID_high_SF.append(self.pixel_veto_mva_th1f.GetBinContent(2))
            Photon_PixelVeto_MVAID_high_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(2))
          else:
            Photon_PixelVeto_MVAID_low_SF.append(self.pixel_veto_mva_th1f.GetBinContent(3))
            Photon_PixelVeto_MVAID_low_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(3))
        else:
          Photon_PixelVeto_MVAID_Inc_SF.append(self.pixel_veto_mva_th1f.GetBinContent(4))
          Photon_PixelVeto_MVAID_Inc_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(4))
          if photons[ipho].r9>0.96:
            Photon_PixelVeto_MVAID_high_SF.append(self.pixel_veto_mva_th1f.GetBinContent(5))
            Photon_PixelVeto_MVAID_high_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(5))
          else:
            Photon_PixelVeto_MVAID_low_SF.append(self.pixel_veto_mva_th1f.GetBinContent(6))
            Photon_PixelVeto_MVAID_low_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(6))

      else: 
        Photon_MVA90_ID_SF.append(self.id_mva_th2f.GetBinContent(self.id_mva_th2f.FindBin(photons[ipho].eta, 499)))
        Photon_MVA90_ID_SFerr.append(self.id_mva_th2f.GetBinError(self.id_mva_th2f.FindBin(photons[ipho].eta, 499)))
        if abs(photons[ipho].eta)<1.566:
          Photon_PixelVeto_MVAID_Inc_SF.append(self.pixel_veto_mva_th1f.GetBinContent(1))
          Photon_PixelVeto_MVAID_Inc_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(1))
          if photons[ipho].r9>0.96:
            Photon_PixelVeto_MVAID_high_SF.append(self.pixel_veto_mva_th1f.GetBinContent(2))
            Photon_PixelVeto_MVAID_high_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(2))
          else:
            Photon_PixelVeto_MVAID_low_SF.append(self.pixel_veto_mva_th1f.GetBinContent(3))
            Photon_PixelVeto_MVAID_low_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(3))
        else:
          Photon_PixelVeto_MVAID_Inc_SF.append(self.pixel_veto_mva_th1f.GetBinContent(4))
          Photon_PixelVeto_MVAID_Inc_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(4))
          if photons[ipho].r9>0.96:
            Photon_PixelVeto_MVAID_high_SF.append(self.pixel_veto_mva_th1f.GetBinContent(5))
            Photon_PixelVeto_MVAID_high_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(5))
          else:
            Photon_PixelVeto_MVAID_low_SF.append(self.pixel_veto_mva_th1f.GetBinContent(6))
            Photon_PixelVeto_MVAID_low_SFerr.append(self.pixel_veto_mva_th1f.GetBinError(6))

    self.out.fillBranch('Photon_MVA90_ID_SF', Photon_MVA90_ID_SF)
    self.out.fillBranch('Photon_MVA90_ID_SFerr', Photon_MVA90_ID_SFerr)
    self.out.fillBranch('Photon_PixelVeto_MVAID_Inc_SF', Photon_PixelVeto_MVAID_Inc_SF)
    self.out.fillBranch('Photon_PixelVeto_MVAID_Inc_SFerr', Photon_PixelVeto_MVAID_Inc_SFerr)
    self.out.fillBranch('Photon_PixelVeto_MVAID_high_SF', Photon_PixelVeto_MVAID_high_SF)
    self.out.fillBranch('Photon_PixelVeto_MVAID_high_SFerr', Photon_PixelVeto_MVAID_high_SFerr)
    self.out.fillBranch('Photon_PixelVeto_MVAID_low_SF', Photon_PixelVeto_MVAID_low_SF)
    self.out.fillBranch('Photon_PixelVeto_MVAID_low_SFerr', Photon_PixelVeto_MVAID_low_SFerr)

    return True

PhoIDSF2016apv = lambda: PhoIDSFProducer("2016apv")
PhoIDSF2016 = lambda: PhoIDSFProducer("2016")
PhoIDSF2017 = lambda: PhoIDSFProducer("2017")
PhoIDSF2018 = lambda: PhoIDSFProducer("2018")
