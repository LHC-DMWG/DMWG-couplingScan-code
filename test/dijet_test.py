import json
import numpy as np
import matplotlib.pyplot as plt
import os

from couplingscan.scan import *
from couplingscan.rescaler import *
from couplingscan.limitparsers import *

# Analysing results from http://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-16-056/

plot_tag = ""

analysis_tag = "CMS-EXO-16-056"

def clean_grid(xvals, yvals, zvals) :
  if zvals.size == 0 :
    return np.array([]), np.array([]), np.array([])
  xclean = xvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  yclean = yvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  zclean = zvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  return xclean, yclean, zclean

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

def make_plot(xvals, yvals, zvals, this_tag, addText=None, addCurves=None, addPoints=False) :

  levels = range(26)  # Levels must be increasing.
  fig,ax=plt.subplots(1,1)
  plt.xlim(500, 3500)
  plt.ylim(0, 1700)
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
    ax.scatter(xnon,ynon,color='red', marker='o',facecolors='none',linewidths=1,s=1)
    ax.scatter(xexcl,yexcl,color='white', marker='o',facecolors='none',linewidths=1,s=1)

  ax.set_xlabel("m$_{ZA}$ [GeV]")
  ax.set_ylabel("m$_{\chi}$ [GeV]")   

  # Now add exclusion contour
  ax.tricontour(xvals, yvals, zvals,levels=[1],colors=['w'],linewidths=[2])

  # Now add another for comparison if desired.
  if addCurves :
    for curve in addCurves :
      ax.add_line(curve)

  # Add text
  if addText :
    plt.figtext(0.2,0.75,addText,size=14,bbox=dict(facecolor='white', alpha=0.75))
    #plt.figtext(0.2,0.75,addText,backgroundcolor="white",size=14)

  plt.savefig('plots/{0}_{1}.eps'.format(analysis_tag,this_tag),bbox_inches='tight')
  plt.savefig('plots/{0}_{1}.pdf'.format(analysis_tag,this_tag),bbox_inches='tight')

# Create scan beginning from gq limit.
# Extract HEPData into useable format
with open("dijet_hepdata/hepdata_gqplot_cms36ifb.json", "r") as read_file:
  data = json.load(read_file)
invalues = data["values"]
# And convert to x-y numpy arrays
xlist = np.array([val["x"][0]["value"] for val in invalues]).astype(float)
ylist = np.array([val["y"][0]["value"] for val in invalues]).astype(float)

# Now create a 1d visible limit object with this, and extract our 2d limits from it.
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
make_plot(scan_A1.mmed, scan_A1.mdm, values, "A1", addText=None, addCurves=None, addPoints=True)

# Get the other three scenarios by rescaling from our scan
rescaleA1 = Rescaler(scan_A1, values)
A2_depths = rescaleA1.rescale_by_br_quarks(target_gq=0.1,target_gdm=1,target_gl=0.1,model='axial')[(0.1,1,0.1)]
V_depths = rescaleA1.rescale_by_br_quarks(target_gq=[0.1,0.25],target_gdm=1,target_gl=[0.0,0.01],model='vector')
V1_depths = V_depths[(0.25,1.0,0.0)]
V2_depths = V_depths[(0.1,1.0,0.01)]
make_plot(scan_A1.mmed, scan_A1.mdm, A2_depths, "A2_rescaled", addPoints = True)
make_plot(scan_A1.mmed, scan_A1.mdm, V1_depths, "V1_rescaled", addPoints = True)
make_plot(scan_A1.mmed, scan_A1.mdm, V2_depths, "V2_rescaled", addPoints = True)

# Now get the other three scenarios via the gq plot directly.
# Confirm they are identical.
scan_A2 = DMAxialModelScan(mmed=target_xgrid.flatten(), mdm=target_ygrid.flatten(), gq=0.1, gdm=1.0, gl=0.1)
A2_direct = gq_limit.extract_exclusion_depths(scan_A2)
scan_V1 = DMVectorModelScan(mmed=target_xgrid.flatten(), mdm=target_ygrid.flatten(), gq=0.25, gdm=1.0, gl=0.0)
V1_direct = gq_limit.extract_exclusion_depths(scan_V1)
scan_V2 = DMVectorModelScan(mmed=target_xgrid.flatten(), mdm=target_ygrid.flatten(), gq=0.1, gdm=1.0, gl=0.01)
V2_direct = gq_limit.extract_exclusion_depths(scan_V2)
make_plot(scan_A2.mmed, scan_A2.mdm, A2_direct, "A2_from1dlimit", addPoints = True)
make_plot(scan_V1.mmed, scan_V1.mdm, V1_direct, "V1_from1dlimit", addPoints = True)
make_plot(scan_V2.mmed, scan_V2.mdm, V2_direct, "V2_from1dlimit", addPoints = True)

# Alternatively, we can create scan beginning from cross-section limit by just making gq limit (one line).
# This is what we'll show in the paper example.

# We don't have public data that includes both an observed and theory curve though, so theory curve was extracted
# from fig 11a using WebPlotDigitizer
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
plt.clf()
plt.plot(xlist_xsec,ylist_xsec)
plt.plot(x_theory, y_theory)
plt.xlim(500, 8000)
plt.yscale('log')
plt.ylim(1e-5, 1e3)
plt.savefig('plots/validation/dijet_check_1dinputs.pdf',bbox_inches='tight')
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
    make_plot(x, y, z, "A1_fromxsec", addText=text, addCurves=draw_contours, addPoints=True)
  else : 
    text = "Vector\ng$_q$=0.25, g$_\chi$=1.0, g$_l$=0.0"
    values_new_V1 = xsec_limit.extract_exclusion_depths(scan_V1)
    x, y, z = clean_grid(scan_V1.mmed, scan_V1.mdm, values_new_V1)
    make_plot(x, y, z, "V1_fromxsec", addText=text, addCurves=draw_contours, addPoints=True)