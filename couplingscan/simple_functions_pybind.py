import lhapdfwrap as pdfwrap

from exclCalc import *
import numpy as np
import scipy as sc
import scipy.integrate as integrate
import math

## Believe that full cross section is not really analytically integrable.
## New attempt is to do it numerically per point.
ECM = 13000.**2

wrapper = pdfwrap.IntegrandHandler("NNPDF30_nlo_as_0118", ECM)

## Integrate and multiply by couplings to obtain full cross sections
def relative_monox_xsec_parton_vector(M, mDM, gq, gDM, gl) :
  gamma = totalWidthVector(M, mDM, gq, gDM, gl)
  # Integrate dSig/ds over kinematic limits: s = 4mDM^2 to CME
  # Give integration some points of interest to check, to help
  # it discover the resonance peak.
  intpoints = [M,M**2-gamma,M**2,M**2+gamma]
  integral = integrate.quad(wrapper.integrand_parton_vector,4.*mDM**2,ECM,args=(gamma,M,mDM),points=intpoints,limit=500)
  # Use when checking integral precision
  #print("Est.",integral[0],"pm",integral[1])  
  xsec = gq**2 * gDM**2 * integral[0]
  return xsec

## Integrate and multiply by couplings to obtain full cross sections
def relative_monox_xsec_parton_axial(M, mDM, gq, gDM, gl) :
  gamma = totalWidthAxial(M, mDM, gq, gDM, gl)
  # Integrate dSig/ds over kinematic limits: s = 4mDM^2 to CME
  # Give integration some points of interest to check, to help
  # it discover the resonance peak.
  intpoints = [M,M**2-gamma,M**2,M**2+gamma]
  integral = integrate.quad(wrapper.integrand_parton_axialvector,4.*mDM**2,ECM,args=(gamma,M,mDM),points=intpoints,limit=500)
  # Use when checking integral precision
  #print("Est.",integral[0],"pm",integral[1])
  xsec = gq**2 * gDM**2 * integral[0]
  return xsec

# Lower limit is actually a curve
# Upper limit is 1
def limit_x1(x2,pid,gamma,M,mDM) :
  lower_lim = (4.*mDM**2)/(x2 * ECM)
  return [lower_lim, 1]

# This is from its smallest value when x1 is largest
# and goes up to 1.
def limit_x2(pid,gamma,M,mDM) :
  lower_lim = (4.*mDM**2)/ECM
  return [lower_lim, 1]

# Here is where I can make it check the
# various points in S that allow the integral
# to converge correctly.
# Add the point at which the integrand goes to zero
# since that's a discontinuity now.
def opts_x1(x2,pid,gamma,M,mDM) :
  points_list = [M**2/(x2 * ECM)]  
  # Don't both returning ridge location if we're off-shell
  if (M < 2.*mDM) :
    return {}
  else :
    return {'points' : points_list}

# Nothing special.
def opts_x2(pid,gamma,M,mDM) :
  return {}

# Only doing first 2 generations since everything else is negligible
def relative_monox_xsec_hadron_vector(M, mDM, gq, gDM, gl) :
  gamma = totalWidthVector(M, mDM, gq, gDM, gl)
  xsec = 0
  # Note: only doing up and down PDFs because I tested
  # with all 4 and saw no discernable difference.
  for q_pid in range(1,3) :
    # Variables of integration are x1 and x2
    # Limits and options are functions.
    # The limits are similarly functions.   
    integral = integrate.nquad(wrapper.integrand_hadronic_vector,[limit_x1,limit_x2],args=(q_pid,gamma,M,mDM),opts=[opts_x1,opts_x2])    
    xsec = xsec + integral[0]
  return gq**2 * gDM**2 * xsec

def relative_monox_xsec_hadron_axial(M, mDM, gq, gDM, gl) :
  gamma = totalWidthAxial(M, mDM, gq, gDM, gl)
  xsec = 0
  # Note: only doing up and down PDFs because I tested
  # with all 4 and saw no discernable difference.
  for q_pid in range(1,3) :
    # Variables of integration are x1 and x2
    # Limits and options are functions.
    # The limits are similarly functions.    
    integral = integrate.nquad(wrapper.integrand_hadronic_axialvector,[limit_x1,limit_x2],args=(q_pid,gamma,M,mDM),opts=[opts_x1,opts_x2],full_output=True) 
    xsec = xsec + integral[0]
  return gq**2 * gDM**2 * xsec

# Functions which use only integral of propagator.
# Fails in V <-> AV translation; too simple.
# However, seems quite successful within a scenario.
# Use these after conversion from one model to another.
def relative_monox_propagator_integral_vector(M, mDM, gq, gDM, gl) :  

  gamma = totalWidthVector(M, mDM, gq, gDM, gl)
  arctan_factor = math.pi/2.0 + np.arctan((M**2 - 4.*mDM**2)/(M*gamma))
  sigma = gq**2 * gDM**2 * arctan_factor/(M*gamma)
  return sigma

def relative_monox_propagator_integral_axial(M, mDM, gq, gDM, gl) :  

  gamma = totalWidthAxial(M, mDM, gq, gDM, gl)
  arctan_factor = math.pi/2.0 + np.arctan((M**2 - 4.*mDM**2)/(M*gamma))
  sigma = gq**2 * gDM**2 * arctan_factor/(M*gamma)  
  return sigma
