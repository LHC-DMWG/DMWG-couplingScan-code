import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import scipy as sp
import time

from couplingscan.scan import *
from couplingscan.rescaler import *

# Analysing results from 
# https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/EXOT-2016-32/

# Add this to any plots.
plot_tag = ""

analysis_tag = "EXOT-2016-32"
scenario_tag = "AV_gq0p25_gchi1p0"
paper_scenario = {"model" : "AV", "gq"  : 0.25, "gdm" : 1.0, "gl" : 0.0}

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

def get_official_contour(scenario_name) :
  # Now get official contour
  with open("monophoton_hepdata/hepdata_{0}_obscontour.json".format(scenario_name), "r") as read_file:
    contour_data = json.load(read_file)
  contour_values = contour_data["values"]
  contour_xlist = np.array([point["x"][0]["value"] for point in contour_values]).astype(float)
  contour_ylist = np.array([point["y"][0]["value"] for point in contour_values]).astype(float)

  # Now make a curve
  vertices = np.column_stack((contour_xlist, contour_ylist))
  string_path = mpath.Path(vertices)
  patch = mpatches.PathPatch(string_path, edgecolor="red", facecolor=None, fill=False, lw=2)

  return patch

def make_plot(xvals, yvals, zvals, this_tag, addText=None, addCurves=None, addPoints=False) :

  levels = range(26)  # Levels must be increasing.
  fig,ax=plt.subplots(1,1)
  plt.xlim(0, 1200)
  plt.ylim(0, 500)
  plt.rc('font',size=17)
  ratio = get_aspect_ratio(ax)
  ax.set_aspect(ratio)
  cp = ax.tricontourf(xvals, yvals, zvals, levels=levels, cmap='Blues_r')
  fig.colorbar(cp)

  # Want points under contour, if adding them.
  if addPoints :
    # Separate into two populations: excluded and non excluded.
    xexcl,yexcl = [],[]
    xnon,ynon = [],[]
    for x,y,z in zip(xvals,yvals,zvals) :
      if z < 1. : 
        xexcl.append(x)
        yexcl.append(y)
      else :
        xnon.append(x)
        ynon.append(y)
    #for i, j, k in zip(xvals,yvals,zvals) :
    ax.scatter(xnon,ynon,color='red', marker='o',facecolors='none',linewidths=2)
    ax.scatter(xexcl,yexcl,color='white', marker='o',facecolors='none',linewidths=2)

  ax.set_xlabel("m$_{ZA}$ [GeV]")
  ax.set_ylabel("m$_{\chi}$ [GeV]")   

  # Now add exclusion contour (if not doing official - harder to see with both)
  if not addCurves :
    ax.tricontour(xvals, yvals, zvals,levels=[1],colors=['w'],linewidths=[2])

  # Now add another for comparison if desired.
  if addCurves :
    for curve in addCurves :
      ax.add_patch(curve)

  # Add text
  if addText :
    plt.figtext(0.2,0.75,addText,size=14)
    #plt.figtext(0.2,0.75,addText,backgroundcolor="white",size=14)

  plt.savefig('plots/{0}_{1}.eps'.format(analysis_tag,this_tag),bbox_inches='tight')
  plt.savefig('plots/{0}_{1}.pdf'.format(analysis_tag,this_tag),bbox_inches='tight')

# Extract HEPData into useable format
with open("monophoton_hepdata/hepdata_{0}.json".format(scenario_tag), "r") as read_file:
  data = json.load(read_file)
values = data["values"]

# Extract as numpy arrays
# Note: to run on monojet, swap the [1] and [0] in first two rows.
xlist = np.array([val["x"][1]["value"] for val in values]).astype(float)
ylist = np.array([val["x"][0]["value"] for val in values]).astype(float)
zlist = np.array([val["y"][0]["value"] for val in values]).astype(float)

# Make a contour plot matching the one from the analysis
text_original="Original\nAxial-vector\ng$_{q}$=0.25, g$_{\chi}$=1"
make_plot(xlist,ylist,zlist,scenario_tag+"_original",addText=text_original,addPoints=True)

# And one with formal line
contour_A1 = get_official_contour("A1")
make_plot(xlist,ylist,zlist,scenario_tag+"_original_withcontour",addText=text_original,addCurves=[contour_A1],addPoints=True)

# Make a scan containing this grid
scan_A1 = DMAxialModelScan(
mmed=xlist,
mdm=ylist,
gq=paper_scenario["gq"],
gdm=paper_scenario["gdm"],
gl=paper_scenario["gl"],
)
rescaler_fromA1 = Rescaler(scan_A1)

# Now convert to each of our three other scenarios and see how it looks.
# Do each in a different way.
new_scenarios = {
  "AV_gq0p1_gl0p1_gchi1p0" : {"model" : "AV", "gq" : 0.1, "gdm" : 1.0, "gl" : 0.1, "scenario_tag" : "A2"},
  "V_gq0p25_gchi1p0" : {"model" : "V", "gq"  : 0.25, "gdm" : 1.0, "gl" : 0.0, "scenario_tag" : "V1"},
  "V_gq0p1_gl0p01_gchi1p0" : {"model" : "V", "gq" : 0.1, "gdm" : 1.0, "gl" : 0.01, "scenario_tag" : "V2"}
}

# First: rescale to another AV model.
# For this let's use the propagator method.
scalefactors_A1toA2 = rescaler_fromA1.rescale_by_propagator(new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["gq"], new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["gdm"], new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["gl"])
# Extract the relevant entry (only one here)
scalefactors_A1toA2 = scalefactors_A1toA2[(new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["gq"], new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["gdm"], new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["gl"])]
# Z values are weird, they are actually like 1/z value. So here I want to divide by my scale factors,
# not multiply by them.
limits_A2 = zlist/scalefactors_A1toA2
new_scenarios["AV_gq0p1_gl0p1_gchi1p0"]["limits"] = limits_A2
print("Finished A2 limits")

# Second: scale to a vector model.
# Here we need to use the hadron level cross section.
scalefactors_A1toV1 = rescaler_fromA1.rescale_by_hadronic_xsec_monox(new_scenarios["V_gq0p25_gchi1p0"]["gq"],new_scenarios["V_gq0p25_gchi1p0"]["gdm"],new_scenarios["V_gq0p25_gchi1p0"]["gl"],'vector')
scalefactors_A1toV1 = scalefactors_A1toV1[(new_scenarios["V_gq0p25_gchi1p0"]["gq"],new_scenarios["V_gq0p25_gchi1p0"]["gdm"],new_scenarios["V_gq0p25_gchi1p0"]["gl"])]
limits_V1 = zlist/scalefactors_A1toV1
new_scenarios["V_gq0p25_gchi1p0"]["limits"] = limits_V1
print("Finished V1 limits")

# Third: to get to V2, we should use
# propagator scaling from V1. For this let's
# create a new scan and a new rescaler.
scan_V1 = DMVectorModelScan(
mmed=xlist,
mdm=ylist,
gq=new_scenarios["V_gq0p25_gchi1p0"]["gq"],
gdm=new_scenarios["V_gq0p25_gchi1p0"]["gdm"],
gl=new_scenarios["V_gq0p25_gchi1p0"]["gl"],
)
rescaler_fromV1 = Rescaler(scan_V1)
scalefactors_V1toV2 = rescaler_fromV1.rescale_by_propagator(new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["gq"],new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["gdm"],new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["gl"],'vector')
scalefactors_V1toV2 = scalefactors_V1toV2[(new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["gq"],new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["gdm"],new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["gl"])]
limits_V2 = limits_V1/scalefactors_V1toV2
new_scenarios["V_gq0p1_gl0p01_gchi1p0"]["limits"] = limits_V2

# Now ready to make plots.
from matplotlib.path import Path
for name, this_scenario in new_scenarios.items() :

  # Baseline plot with homemade contour
  text_here = "Rescaled\n{0}\ng$_{4}$={1}, g$_{5}$={2}, g$_{6}$={3}".format(("Axial-vector" if this_scenario["model"]=="AV" else "Vector"),this_scenario["gq"],this_scenario["gdm"],this_scenario["gl"],"q","\chi","l")
  make_plot(xlist,ylist,this_scenario["limits"],name+plot_tag,addText=text_here,addPoints=True)

  contour = get_official_contour(this_scenario["scenario_tag"])

  text_here = "Rescaled\n{0}\ng$_{4}$={1}, g$_{5}$={2}, g$_{6}$={3}".format(("Axial-vector" if this_scenario["model"]=="AV" else "Vector"),this_scenario["gq"],this_scenario["gdm"],this_scenario["gl"],"q","\chi","l")
  make_plot(xlist,ylist, this_scenario["limits"],name+"_compare"+plot_tag, addText=text_here, addCurves=[contour],addPoints=True)
