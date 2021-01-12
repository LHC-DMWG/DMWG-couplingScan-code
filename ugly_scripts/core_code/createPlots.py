from __future__ import division
import ROOT
from array import array
from exclCalc import *
import json
import sys

#==================== Load arguments ====================

# Needed configuration file on command line. Complain if you didn't get one.
if len(sys.argv) < 2 :
    print("You need to provide a configuration file to run with!")
    exit(1)

config_file = sys.argv[1]
print("Running with config",config_file)

#================ Load parameters and data ===============

with open(config_file) as json_file:
    config = json.load(json_file)

infile = ROOT.TFile.Open(config["input"]["file"],"READ")
graphname = str(config["input"]["graph"])
data = infile.Get(graphname)

# Use default range of mDM unless otherwise specified
mDM_upper = 2000
mDM_lower = 0
if "mDM_upper" in config["options"].keys() :
    mDM_upper = config["options"]["mDM_upper"]
if "mDM_lower" in config["options"].keys() :
    mDM_lower = config["options"]["mDM_lower"]

# Use range of masses available in input file
# or specified range, whichever is smaller
xvals = []
for index in range(data.GetN()) :
    xval = ROOT.Double()
    yval = ROOT.Double()    
    data.GetPoint(index,xval,yval)
    xvals.append(int(xval))
xvals = sorted(xvals)
mMed_upper = xvals[-1]
mMed_lower = xvals[0]
print("Before, lower x value was",mMed_lower)
if "mMed_upper" in config["options"].keys() :
    mMed_upper = min(mMed_upper,config["options"]["mMed_upper"])
if "mMed_lower" in config["options"].keys() :
    mMed_lower = max(mMed_lower,config["options"]["mMed_lower"])
    print("After, lower x value is",mMed_lower)

print("Using mMed range",mMed_lower,"-",mMed_upper)
print("And mDM range",mDM_lower,"-",mDM_upper)

scaleby = 1000.0 if "GeV" in str(config["input"]["units"]) else 1.0


#========= Create one output graph per scenario ==========


for scenario in config["scenarios"] :

    exclusion_list = []
    exclusion_hist = ROOT.TH2D("exclusion_hist","exclusion_hist", int(float(mMed_upper-mMed_lower)/10.0), mMed_lower, mMed_upper, int(float(mDM_upper-mDM_lower)/10.0), mDM_lower, mDM_upper) 
    exclusion_hist.SetDirectory(0)

    if str(scenario["model"]) == "AV" :
        print "Treating this as axial-vector."
    else :
        print "Treating this as vector."

    gq = scenario["gq"]
    gDM = scenario["gDM"]
    gl = scenario["gl"]

    # Row of empty points below and above contour, ensuring
    # that lower edge is drawn correctly. Do the same
    # at the top edge of the plot.
    for m in range(mDM_lower, mDM_upper+1, 10) :
      exclusion_list.append((mMed_lower-10,m,25))
      exclusion_list.append((mMed_upper+10,m,25))
    for M in range(mMed_lower, mMed_upper+1, 10) :
      exclusion_list.append((M, mDM_upper+10,25))

    # Look at individual points
    for M in range(mMed_lower, mMed_upper+1, 10):

        for m in range(mDM_lower, mDM_upper+1, 10):

            # Get corresponding g' to mDM
            if str(scenario["model"]) == "AV" :
                gPrime = gqPrimeAxial(M, m, gq, gDM, gl)
            else :
                gPrime = gqPrimeVector(M, m, gq, gDM, gl)

            # Exclusion depth
            analysis_limit = data.Eval(M) 
            depth = analysis_limit/gPrime
            exclusion_list.append((M, m, depth))
        
            # Exclusion histograms
            # Excluded iff larger than observed coupling limit
            if gPrime > analysis_limit :
                exclusion_hist.Fill(M,m)

    # Arrays to turn into a TGraph
    exclusion_list = sorted(exclusion_list)
    xPoints = array('f', [i[0]/scaleby for i in exclusion_list])
    yPoints = array('f', [i[1]/scaleby for i in exclusion_list])
    zPoints = array('f', [i[2] for i in exclusion_list])

    lowestExclMasses = ROOT.TGraph2D(len(xPoints),xPoints, yPoints, zPoints)
    savefile_name = str(config["output"]["file"])
    if ".root" in savefile_name :
        savefile_name = savefile_name.replace(".root",scenario["name"]+".root")
    else :
        savefile_name = savefile_name+"_"+scenario["name"]+".root"
    savefile = ROOT.TFile.Open(savefile_name, "RECREATE")
    lowestExclMasses.Write("Graph2D")
    exclusion_hist.Write()
    savefile.Close()

    c0 = ROOT.TCanvas("c0","c0",800,800)
    exclusion_hist.SetTitle("{0} ;mediator mass (TeV); DM mass (TeV);".format(scenario["name"]))
    exclusion_hist.Draw("colz")
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    ROOT.gStyle.SetOptStat(0)
    params = "gq="+str(gq)+"; gl="+str(gl)+"; gchi="+str(gDM)
    t = ROOT.TLatex()
    t.DrawLatexNDC(0.5,0.8, params)
    c0.SaveAs(savefile_name.replace(".root",".png"))
    del c0
