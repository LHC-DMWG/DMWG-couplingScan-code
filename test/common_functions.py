# Functions unrelated to couplingscan package usage but helpful for other components of demo scripts.

# Supports plot making
def get_aspect_ratio(ax) :
  ratio = 1.0
  xleft, xright = ax.get_xlim()
  ybottom, ytop = ax.get_ylim()
  return abs((xright-xleft)/(ybottom-ytop))*ratio

# Just a plotter.
import matplotlib.pyplot as plt
import os
def make_plot(xvals, yvals, zvals, xaxisLims, yaxisLims, analysis_tag, this_tag, addText=None, addCurves=None, addPoints=False) :
  
  # Check for output dir; make it if missing
  if not os.path.exists("plots/validation") :
    os.makedirs("plots/validation")

  levels = range(26)  # Levels must be increasing.
  fig,ax=plt.subplots(1,1)
  plt.xlim(xaxisLims[0], xaxisLims[1])
  plt.ylim(yaxisLims[0], yaxisLims[1])
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
      ax.add_line(curve)

  # Add text
  if addText :
    plt.figtext(0.2,0.75,addText,size=14,bbox=dict(facecolor='white', alpha=0.75))
    #plt.figtext(0.2,0.75,addText,backgroundcolor="white",size=14)

  print("Making plot",'plots/{0}_{1}.pdf'.format(analysis_tag,this_tag))
  plt.savefig('plots/{0}_{1}.pdf'.format(analysis_tag,this_tag),bbox_inches='tight')

# A cross section limit plot.
def make_xsec_plot(xvals_obs, yvals_obs, xvals_th, yvals_th, analysis_tag, this_tag) :
  # Check for output dir; make it if missing
  if not os.path.exists("plots/validation") :
    os.makedirs("plots/validation")

  plt.clf()
  plt.plot(xvals_obs,yvals_obs,label="observed")
  plt.plot(xvals_th,yvals_th,label="theory")
  plt.legend()
  saveas = 'plots/validation/{0}_{1}.pdf'.format(analysis_tag,this_tag)
  print("Making plot",saveas)
  plt.savefig(saveas,bbox_inches='tight')

# Cleans up numpy arrays to get rid of nans and infs
import numpy as np
def clean_grid(xvals, yvals, zvals) :
  if zvals.size == 0 :
    return np.array([]), np.array([]), np.array([])
  xclean = xvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  yclean = yvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  zclean = zvals[np.logical_and(np.logical_and(np.isfinite(xvals),np.isfinite(yvals)),np.isfinite(zvals))]
  return xclean, yclean, zclean 