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

analysis_tag = "CMS-EXO-16-056"

# Extract HEPData into useable format
with open("dijet_hepdata/hepdata_gqplot_cms36ifb.json", "r") as read_file:
  data = json.load(read_file)
invalues = data["values"]
# And convert to x-y numpy arrays
xlist = np.array([val["x"][0]["value"] for val in invalues]).astype(float)
ylist = np.array([val["y"][0]["value"] for val in invalues]).astype(float)

# Create scan beginning from gq limit.
# Create a 1d visible limit object with this, and extract our 2d limits from it.
# This is based on what the settings are for the world in which the gq plot is made.
gq_limit = CouplingLimit_Dijet(
    mmed=xlist,
    gq_limits=ylist,
    mdm=10000,
    gdm=0.0,
    gl=0.0,
    coupling='vector'
)

# This is what we want to get our limits in: a scan in A1 scenario with plenty of points.
target_xvals = np.linspace(200,3600,35)
target_yvals = np.linspace(0,1700,35)
print(target_yvals)
target_xgrid, target_ygrid = np.meshgrid(target_xvals,target_yvals)
scan_A1 = DMAxialModelScan(
mmed=target_xgrid.flatten(),
mdm=target_ygrid.flatten(),
gq=0.25,
gdm=1.0,
gl=0.0,
)

values = gq_limit.extract_exclusion_depths(scan_A1)
print(values)

# Make a plot
make_plot(scan_A1.mmed, scan_A1.mdm, values, [500, 3500], [0, 1700], analysis_tag, "A1", addText=None, addCurves=None, addPoints=True)

# Get the other three scenarios by rescaling from our scan
rescaleA1 = Rescaler(scan_A1, values)
A2_depths = rescaleA1.rescale_by_br_quarks(target_gq=0.1,target_gdm=1,target_gl=0.1,model='axial')[(0.1,1,0.1)]
V_depths = rescaleA1.rescale_by_br_quarks(target_gq=[0.1,0.25],target_gdm=1,target_gl=[0.0,0.01],model='vector')
V1_depths = V_depths[(0.25,1.0,0.0)]
V2_depths = V_depths[(0.1,1.0,0.01)]
make_plot(scan_A1.mmed, scan_A1.mdm, A2_depths, [500, 3500], [0, 1700], analysis_tag, "A2_rescaled", addPoints = True)
make_plot(scan_A1.mmed, scan_A1.mdm, V1_depths, [500, 3500], [0, 1700], analysis_tag, "V1_rescaled", addPoints = True)
make_plot(scan_A1.mmed, scan_A1.mdm, V2_depths, [500, 3500], [0, 1700], analysis_tag, "V2_rescaled", addPoints = True)

# Now get the other three scenarios via the gq plot directly.
# Confirm they are identical.
scan_A2 = DMAxialModelScan(mmed=target_xgrid.flatten(), mdm=target_ygrid.flatten(), gq=0.1, gdm=1.0, gl=0.1)
A2_direct = gq_limit.extract_exclusion_depths(scan_A2)
scan_V1 = DMVectorModelScan(mmed=target_xgrid.flatten(), mdm=target_ygrid.flatten(), gq=0.25, gdm=1.0, gl=0.0)
V1_direct = gq_limit.extract_exclusion_depths(scan_V1)
scan_V2 = DMVectorModelScan(mmed=target_xgrid.flatten(), mdm=target_ygrid.flatten(), gq=0.1, gdm=1.0, gl=0.01)
V2_direct = gq_limit.extract_exclusion_depths(scan_V2)
make_plot(scan_A2.mmed, scan_A2.mdm, A2_direct, [500, 3500], [0, 1700], analysis_tag, "A2_from1dlimit", addPoints = True)
make_plot(scan_V1.mmed, scan_V1.mdm, V1_direct, [500, 3500], [0, 1700], analysis_tag, "V1_from1dlimit", addPoints = True)
make_plot(scan_V2.mmed, scan_V2.mdm, V2_direct, [500, 3500], [0, 1700], analysis_tag, "V2_from1dlimit", addPoints = True)

# Alternatively, we can create scan beginning from cross-section limit by just making gq limit (one line).
# This is what we'll show in the paper example.

# We don't have public data that includes both an observed and theory curve though, 
# so theory curve was extracted from fig 11a using WebPlotDigitizer
xvals = []
yvals = []
with open("dijet_hepdata/approximate_theorycurve.txt","r") as read_file :
  lines = read_file.readlines()
  for line in lines :
    tokens = line.split(", ")
    xvals.append(float(tokens[0]))
    yvals.append(float(tokens[1]))
x_theory = np.array(xvals)
y_theory = np.array(yvals)

# Now get the 1D observed cross section limit for 11a.
with open("dijet_hepdata/hepdata_crosssectionlimit_cms36ifb.json", "r") as read_file:
  data = json.load(read_file)
xsecvalues = data["values"]
xlist_xsec = np.array([val["x"][0]["value"] for val in xsecvalues]).astype(float)
ylist_xsec = np.array([val["y"][2]["value"] for val in xsecvalues]).astype(float)

# Plot to check it ...
make_xsec_plot(xlist_xsec,ylist_xsec,x_theory,y_theory,analysis_tag,"check_1dinputs")
# And it looks fine. So we can proceed and use these two lines to create a dijet scan.

# And set this as input for the new scan.
xsec_limit = CrossSectionLimit_Dijet(
  mmed_limit=xlist_xsec,
  xsec_limit=ylist_xsec,
  mmed_theory=x_theory,
  xsec_theory=y_theory,
  mdm=1,
  gq=0.25,
  gdm=1.0,
  gl=0.0,
  coupling='vector',
  max_intrinsic_width=0.15
)

# Get paper exclusion limits and make them into contours.
for scanname in ["a1","v1"] :
  draw_contours = []
  for contour_i in [1,2] :
    contour_x = []
    contour_y = []
    with open("dijet_hepdata/{0}_contour_{1}.txt".format(scanname,contour_i),"r") as read_file :
      lines = read_file.readlines()
      for line in lines :
        tokens = line.split(", ")
        contour_x.append(float(tokens[0]))
        contour_y.append(float(tokens[1]))
    thisline = plt.Line2D(contour_x, contour_y,lw = 2, color ='red')
    draw_contours.append(thisline)
  
  if "a1" in scanname : 
    text = "Axial-vector\ng$_q$=0.25, g$_\chi$=1.0, g$_l$=0.0"
    values_new_A1 = xsec_limit.extract_exclusion_depths(scan_A1)
    x, y, z = clean_grid(scan_A1.mmed, scan_A1.mdm, values_new_A1)
    make_plot(x, y, z, [500, 3500], [0, 1700], analysis_tag, "A1_from_approx_xsec", addText=text, addCurves=draw_contours, addPoints=True)
  else : 
    text = "Vector\ng$_q$=0.25, g$_\chi$=1.0, g$_l$=0.0"
    values_new_V1 = xsec_limit.extract_exclusion_depths(scan_V1)
    x, y, z = clean_grid(scan_V1.mmed, scan_V1.mdm, values_new_V1)
    make_plot(x, y, z, [500, 3500],[0, 1700], analysis_tag, "V1_from_approx_xsec", addText=text, addCurves=draw_contours, addPoints=True)

# These don't match quite as exactly, because this method of extracting the cross section is a little wonky,
# but you can see that the two results are equivalent so long as the inputs are equivalent.
# Use whichever is easiest given your analysis results format.