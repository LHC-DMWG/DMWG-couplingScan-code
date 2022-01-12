import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import scipy as sp
import time

# Analysing results from 
# https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/EXOT-2016-27/

# Add this to any plots.
plot_tag = ""

analysis_tag = "EXOT-2016-27"
scenario_tag = "AV_gq0p25_gchi1p0"
paper_scenario = {"model" : "AV", "gq"  : 0.25, "gDM" : 1.0, "gl" : 0.0}

point_limits = [3000,1500]

# This is for benchmarking. If I just want to test the time taken
# to do a few points and check everything runs, turn this on.
doTest = False

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

def make_plot(xvals, yvals, zvals, this_tag, addText=None, addCurves=None, addPoints=False) :

  levels = range(26)  # Levels must be increasing.
  fig,ax=plt.subplots(1,1)
  plt.xlim(0, 2600)
  plt.ylim(0, 1100)
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
    ax.scatter(xnon,ynon,color='red', marker='o',facecolors='none')
    ax.scatter(xexcl,yexcl,color='white', marker='o',facecolors='none')

  # Now add exclusion contour
  ax.tricontour(xvals, yvals, zvals,levels=[1],colors=['w'],linewidths=[2])
  ax.set_xlabel("m$_{ZA}$ [GeV]")
  ax.set_ylabel("m$_{\chi}$ [GeV]")

  # Now add another for comparison if desired.
  if addCurves :
    for curve in addCurves :
      ax.add_patch(curve)

  # Add text
  if addText :
    plt.figtext(0.2,0.7,addText,backgroundcolor="white")

  plt.savefig('plots/{0}_{1}.eps'.format(analysis_tag,this_tag),bbox_inches='tight')
  plt.savefig('plots/{0}_{1}.pdf'.format(analysis_tag,this_tag),bbox_inches='tight')


# Now begin extracting data.

# HEPData file is a pretty weird json.
with open("hepdata_{0}.json".format(scenario_tag), "r") as read_file:
  data = json.load(read_file)

values = data["values"]

# Extract as numpy arrays
xlist = np.array([val["x"][0]["value"] for val in values]).astype(float)
ylist = np.array([val["x"][1]["value"] for val in values]).astype(float)
zlist = np.array([val["y"][0]["value"] for val in values]).astype(float)

# Make a contour plot matching the one from the analysis
text_original="Original\nAxial-vector\ng$_{q}$=0.25, g$_{\chi}$=1"
make_plot(xlist,ylist,zlist,scenario_tag+"_original",addText=text_original,addPoints=True)

# Now convert to each of our other scenarios and let's see how plausible it looks
from package.simple_functions_pybind import *
scenarios = {
#  "AV_gq0p25_gchi1p0" : {"model" : "AV", "gq"  : 0.25, "gDM" : 1.0, "gl" : 0.0},
#  "AV_gq0p1_gl0p1_gchi1p0" : {"model" : "AV", "gq" : 0.1, "gDM" : 1.0, "gl" : 0.1},
  "V_gq0p25_gchi1p0" : {"model" : "V", "gq"  : 0.25, "gDM" : 1.0, "gl" : 0.0},
#  "V_gq0p1_gl0p01_gchi1p0" : {"model" : "V", "gq" : 0.1, "gDM" : 1.0, "gl" : 0.01}
}

# Conversion - scale by xsec(old)/xsec(new) bc theory_new is on the bottom
for scenario in scenarios.keys() :

  print("Beginning scenario",scenario)
  start = time.perf_counter()

  this_scenario = scenarios[scenario]
  converted_x = []
  converted_y = []
  converted_z = []

  # paper_val = sigma_obs/sigma_theory
  # Want to convert to sigma_obs/sigma_new_theory
  # So multiply by (sigma_theory/sigma_new_theory)
  nPointsKept = 0
  for (mMed, mDM, paper_val) in zip(xlist,ylist,zlist) :

    # Various escape conditions
    if mMed > point_limits[0] or mDM > point_limits[1] : 
      continue
    nPointsKept = nPointsKept+1
    if doTest and nPointsKept > 5 :
      continue
    print("Beginning point",mMed, mDM)

    # Propagator only version: appropriate within scenarios?
    # paper_xsec = relative_monox_propagator_integral_axial(mMed, mDM, paper_scenario["gq"], paper_scenario["gDM"], paper_scenario["gl"])
    # if this_scenario["model"] == "AV" :
    #   this_xsec = relative_monox_propagator_integral_axial(mMed, mDM, this_scenario["gq"], this_scenario["gDM"], this_scenario["gl"])
    # else :
    #   this_xsec = relative_monox_propagator_integral_vector(mMed, mDM, this_scenario["gq"], this_scenario["gDM"], this_scenario["gl"])

    # Most full version
    time_one = time.perf_counter()
    paper_xsec = relative_monox_xsec_hadron_axial(mMed, mDM, paper_scenario["gq"], paper_scenario["gDM"], paper_scenario["gl"])
    stop_one = time.perf_counter()
    print("One integral took",round(stop_one-time_one),"seconds.")
    print(paper_xsec)
    if this_scenario["model"] == "AV" :
      this_xsec = relative_monox_xsec_hadron_axial(mMed, mDM, this_scenario["gq"], this_scenario["gDM"], this_scenario["gl"])
    else :
      time_two = time.perf_counter()
      this_xsec = relative_monox_xsec_hadron_vector(mMed, mDM, this_scenario["gq"], this_scenario["gDM"], this_scenario["gl"])  
      stop_two = time.perf_counter()
      print("Next integral took",round(stop_two-time_two),"seconds.")
    print(this_xsec)    

    # We get zero cross section off-shell with this setup
    if this_xsec == 0 :
      continue
    else :
      converted_val = paper_val*(paper_xsec/this_xsec)
      converted_z.append(converted_val)
      converted_x.append(mMed)
      converted_y.append(mDM)

  stop = time.perf_counter()
  print("Finished in {0} seconds.".format(round(stop-start)))

  # Now ready to make a plot.
  array_x = np.array(converted_x)
  array_y = np.array(converted_y)
  array_z = np.array(converted_z)
  text_here = "{0}\ng$_{4}$={1}, g$_{5}$={2}, g$_{6}$={3}".format(("Axial-vector" if this_scenario["model"]=="AV" else "Vector"),this_scenario["gq"],this_scenario["gDM"],this_scenario["gl"],"q","\chi","l")
  make_plot(array_x,array_y,array_z,scenario+plot_tag,addText=text_here,addPoints=True)

  # and save.
  scenarios[scenario]["xvals"] = array_x
  scenarios[scenario]["yvals"] = array_y
  scenarios[scenario]["zvals"] = array_z

