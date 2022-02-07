import json
import numpy as np
import matplotlib.pyplot as plt

from couplingscan.scan import *
from couplingscan.rescaler import *
from couplingscan.limit_1d import *

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
gq_limit = DMLimit1D_Dijet(
    mmed=xlist,
    gq_limits=ylist,
    mdm=10000,
    gdm=0.0,
    gl=0.0,
    coupling='axial'
)

# This is what we want to get our limits in: a scan in A1 scenario with plenty of points.
target_xvals = np.linspace(200,3600,35)
target_yvals = np.linspace(0,1700,18)
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
make_plot(scan_A1.mmed, scan_A1.mdm, values, analysis_tag, addText=None, addCurves=None, addPoints=True)

# Get the other four scenarios by rescaling from our scan?

# Alternatively, create scan beginning from cross-section limit by just making gq limit (one line)