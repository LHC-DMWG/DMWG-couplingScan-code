import json
import numpy as np
import matplotlib.pyplot as plt
import os

from couplingscan.scan import *
from couplingscan.rescaler import *
from couplingscan.limitparsers import *
from common_functions import *

# Analysing results from http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-16-056/

plot_tag = ""

analysis_tag = "ATLAS-dilepton-internal"

# Create scan beginning from multiple observed limits of different widths
# and a cross section limit (approximate here).
# Extract HEPData into useable format
with open("dilepton_data/hepdata_observed_xseclimits_atlas139ifb.json", "r") as read_file:
  data = json.load(read_file)
values = data["values"]
widths = data["qualifiers"]["RELATIVE WIDTH"]
# And convert to x-y numpy arrays
xlist = np.array([val["x"][0]["value"] for val in values]).astype(float)
ylists = {}
for group in range(6) :
  ylist = np.array([val["y"][group]["value"] for val in values]).astype(float)
  width = widths[group]["value"]
  float_width = float(width.replace(" %",""))
  frac_width = float_width/100
  ylists[frac_width] = ylist

# Read in our theory curve - very approximate but will work for this test
xvals = []
yvals = []
with open("dilepton_data/approximate_theorycurve.txt", "r") as read_file:
  lines = read_file.readlines()
  for line in lines :
    tokens = line.split(", ")
    xvals.append(1000*float(tokens[0])) # this was in TeV
    yvals.append(float(tokens[1]))
x_theory = np.array(xvals)
y_theory = np.array(yvals)

# Now create a visible limit object with this, and extract our 2d limits from it.
# We will give it a full set of observed limits.
# When we do this, we imply that larger intrinsic widths than those passed to it are
# not valid to exclude with this analysis.
# If only one observed limit is given, we treat it as valid everywhere and leave it
# to the user to decide when to cut it off.
dilepton_limit = CrossSectionLimit_Dilepton(
    mmed_limit=xlist,
    xsec_limit=ylists,
    mmed_theory=x_theory,
    xsec_theory=y_theory,
    mdm=2.5,
    gq=0.1,
    gdm=1.0,
    gl=0.01,
    coupling='vector'
)

# Plot to validate.
make_xsec_plot(xlist,ylists[0.03],x_theory,y_theory,analysis_tag,"check_1dinputs")

# A1 and V1 have no lepton coupling so we don't need to worry about them. 
# We'll extract directly to A2 and V2 from our scan.
target_xvals = np.linspace(300,2000,171)
target_yvals = np.linspace(0,1700,35)
target_xgrid, target_ygrid = np.meshgrid(target_xvals,target_yvals)

# Let's start with going straight to V2 since it's literally the same thing.
# The input should be one line through the V2 curve.
scan_V2 = DMVectorModelScan(
mmed=target_xgrid.flatten(),
mdm=target_ygrid.flatten(),
gq=0.1,
gdm=1.0,
gl=0.01,
)
all_depths_V2 = dilepton_limit.extract_exclusion_depths(scan_V2)
values_V2 = dilepton_limit.select_depths(scan_V2,all_depths_V2)
x, y, z = clean_grid(scan_V2.mmed, scan_V2.mdm, values_V2,)
make_plot(x, y, z, [0, 3500], [0, 1700], analysis_tag, "V2_direct", addText=None, addCurves=None, addPoints=True)

# Make an A2 scan
scan_A2 = DMAxialModelScan(
mmed=target_xgrid.flatten(),
mdm=target_ygrid.flatten(),
gq=0.1,
gdm=1.0,
gl=0.1,
)

all_depths_A2 = dilepton_limit.extract_exclusion_depths(scan_A2)
values_A2 = dilepton_limit.select_depths(scan_A2,all_depths_A2)
x, y, z = clean_grid(scan_A2.mmed, scan_A2.mdm, values_A2)

# Draw it!
make_plot(x, y, z, [0, 3500], [0, 1700], analysis_tag, "A2", addText=None, addCurves=None, addPoints=True)

# And get another V2 by rescaling from this scan
print("Inputting values:")
print(all_depths_A2)
rescaleA2 = Rescaler(scan_A2,all_depths_A2)
V2_depths = rescaleA2.rescale_by_br_leptons(target_gq=0.1,target_gdm=1,target_gl=0.01,model='vector')[(0.1,1.0,0.01)]
x, y, z = clean_grid(scan_A2.mmed, scan_A2.mdm, V2_depths)
make_plot(x, y, z, [0, 3500],[0, 1700], analysis_tag, "V2_rescaled", addPoints = True)
