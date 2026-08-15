import ROOT
import os

ROOT.gROOT.SetBatch()

ROOT.gROOT.ProcessLine(".O 2") # Set optimization level to 2
ROOT.gROOT.ProcessLine(".O") # Show optimization level
ROOT.gSystem.SetBuildDir("./obj/",True);

ROOT.gROOT.ProcessLine(".L vetomap_utils.h+O")
ROOT.gROOT.ProcessLine(".L run_histograms.cc+O")

# ROOT.gROOT.ProcessLine(".L NtupleMaker_Functions_GenLevel.h+O")
## ROOT.gROOT.ProcessLine(".L LumiFilterV2.h+O")

