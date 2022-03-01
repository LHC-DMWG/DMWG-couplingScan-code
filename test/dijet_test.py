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

def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

def make_plot(xvals, yvals, zvals, this_tag, addText=None, addCurves=None, addPoints=False) :

  levels = range(26)  # Levels must be increasing.
  fig,ax=plt.subplots(1,1)
  plt.xlim(0, 3500)
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

# Make a plot
make_plot(scan_A1.mmed, scan_A1.mdm, values, "A1", addText=None, addCurves=None, addPoints=True)

# Get the other three scenarios by rescaling from our scan
rescaleA1 = Rescaler(scan_A1)
A2_sfs = rescaleA1.rescale_by_br_quarks(target_gq=0.1,target_gdm=1,target_gl=0.1,model='axial')[(0.1,1,0.1)]
V_sfs = rescaleA1.rescale_by_br_quarks(target_gq=[0.1,0.25],target_gdm=1,target_gl=[0.0,0.01],model='vector')
V1_sfs = V_sfs[(0.25,1.0,0.0)]
V2_sfs = V_sfs[(0.1,1.0,0.01)]
make_plot(scan_A1.mmed, scan_A1.mdm, values/A2_sfs, "A2_rescaled", addPoints = True)
make_plot(scan_A1.mmed, scan_A1.mdm, values/V1_sfs, "V1_rescaled", addPoints = True)
make_plot(scan_A1.mmed, scan_A1.mdm, values/V2_sfs, "V2_rescaled", addPoints = True)

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
# We don't have public data that includes both an observed and theory curve though, so to demonstrate this
# we have to make a theory line to use.

# First let's convert the g_q plot to something like this so we have something to start from.
# We know what the 1D observed limit is, so all we need to do is extract the theory curve.
with open("dijet_hepdata/hepdata_couplinglimit_cms36ifb.json", "r") as read_file:
  data = json.load(read_file)
xsecvalues = data["values"]
xlist_xsec = np.array([val["x"][0]["value"] for val in xsecvalues]).astype(float)
ylist_xsec = np.array([val["y"][2]["value"] for val in xsecvalues]).astype(float)
# Now: theory curve goes like gq2 so theory curve is just
# gqth^2/gqlim^2*observed with all else the same and some fixed gqth.
# Target model will be gq = 0.25, 0, 0 for a mediator mass of 10 TeV (decoupled mediator), vector mediator.
gqlimit_interp = np.interp(xlist_xsec,xlist,ylist,left=10.,right=10.)
theory_curve = ylist_xsec*(0.25**2/gqlimit_interp**2)
# Plot to check it ...
plt.clf()
plt.plot(xlist_xsec,ylist_xsec)
plt.plot(xlist_xsec,theory_curve)
plt.xlim(500, 8000)
plt.yscale('log')
plt.ylim(1e-5, 1e5)
plt.savefig('plots/validation/dijet_extract_theoryxsec.pdf',bbox_inches='tight')
# And it looks fine. So we can proceed and use these two lines to create a dijet scan.

# And set this as input for the new scan.
xsec_limit = CrossSectionLimit_Dijet(
  mmed_limit=xlist_xsec,
  xsec_limit=ylist_xsec,
  mmed_theory=xlist_xsec,
  xsec_theory=theory_curve,
  mdm=10000,
  gq=0.25,
  gdm=0.0,
  gl=0.0,
  coupling='vector',
  max_intrinsic_width=0.15
)

# Now extract some limit plots and compare them to the ones using the original method.
values_new = xsec_limit.extract_exclusion_depths(scan_A1)
# And make a plot
make_plot(scan_A1.mmed, scan_A1.mdm, values_new, "A1_fromxsec", addText=None, addCurves=None, addPoints=True)