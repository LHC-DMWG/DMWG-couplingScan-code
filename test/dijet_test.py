import json
import numpy as np

from couplingscan.scan import *
from couplingscan.rescaler import *
from couplingscan.limit_1d import *

# Analysing results from http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-16-056/

plot_tag = ""

analysis_tag = "CMS-EXO-16-056"

#scenario_tag = "AV_gq0p25_gchi1p0"
#paper_scenario = {"model" : "AV", "gq"  : 0.25, "gdm" : 1.0, "gl" : 0.0}

# Create scan beginning from gq limit.
# Extract HEPData into useable format
with open("dijet_hepdata/hepdata_gqplot_cms36ifb.json", "r") as read_file:
  data = json.load(read_file)
values = data["values"]
# And convert to x-y numpy arrays
xlist = np.array([val["x"][0]["value"] for val in values]).astype(float)
ylist = np.array([val["y"][0]["value"] for val in values]).astype(float)

# Now create a 1d visible limit object with this, and extract our 2d limits from it.
# This model is unclear re vector/axial-vector; let's try both and see if they differ.
# It HAS TO BE leptophobic and decoupled from dark matter; otherwise the cancellations
# in the cross section don't work out, or else in the same coupling scenario as what we're
# trying to go to. Those are really the only options.
# Need to tell the code that so it knows what to do.
gq_limit = DMLimit1D_AxialDijet(
    mmed=xlist,
    gq_limits=ylist,
    mdm=10000,
    gdm=0.0,
    gl=0.0,
)

# This is what we want to get our limits in: a scan in A1 scenario with plenty of points.
target_xvals = np.linspace(200,3600,35)
target_yvals = np.linspace(0,1700,171)
target_xgrid, target_ygrid = np.meshgrid(target_xvals,target_yvals)
scan_A1 = DMAxialModelScan(
mmed=target_xgrid.flatten(),
mdm=target_ygrid.flatten(),
gq=0.25,
gdm=1.0,
gl=0.0,
)

values = gq_limit.extract_exclusion_depths(scan_A1)

# Alternatively, create scan beginning from cross-section limit by just making gq limit (one line)