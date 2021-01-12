from __future__ import division
import math
import numpy as np
from array import array


    
# ======================== Partial widths ============================

# Quark and lepton masses, PDG 2018
quarkMasses = [0.0022, 0.0047, 0.095, 1.275, 4.180, 173.0]
leptonMasses = [0.000511, 0.105658, 1.77686]

def z(M,m): 
    return m**2 / M**2

#DM particles
def width_axial_DM (m_med, m_DM, g_DM):
    
    width = 0
    if(m_med >= 2*m_DM):
        z_DM = z(m_med, m_DM)
        width = (g_DM**2 * m_med)/(12*math.pi) * (1 - 4*z_DM)**(3/2)
    return width

def width_vector_DM (m_med, m_DM, g_DM):
    width = 0
    if(m_med >= 2*m_DM):
        z_DM = z(m_med, m_DM)
        width = (g_DM**2 * m_med)/(12*math.pi) * (1 - 4*z_DM)**(1/2) * (1 + 2*z_DM)
    return width

#quarks
def width_axial_q (m_med, g_q):
    width = 0
    for m_q in quarkMasses:
        if(m_med >= 2*m_q):
            z_q = z(m_med, m_q)
            width += 3*(g_q**2 * m_med)/(12*math.pi) * (1 - 4*z_q)**(3/2)
    return width
def width_vector_q (m_med, g_q):
    width = 0
    for m_q in quarkMasses:
        if(m_med >= 2*m_q):
            z_q = z(m_med, m_q)
            width += 3*(g_q**2 * m_med)/(12*math.pi) * (1 - 4*z_q)**(1/2) * (1 + 2*z_q)
    return width

#charged leptons
def width_axial_l (m_med, g_l):   
    width = 0
    for m_l in leptonMasses:
        if(m_med >= 2*m_l):
            z_l = z(m_med, m_l)
            width += (g_l**2 * m_med)/(12*math.pi) * (1 - 4*z_l)**(3/2)
    return width
def width_vector_l (m_med, g_l):   
    width = 0
    for m_l in leptonMasses:
        if(m_med >= 2*m_l):
            z_l = z(m_med, m_l)
            width += (g_l**2 * m_med)/(12*math.pi) * (1 - 4*z_l)**(1/2) * (1 + 2*z_l)
    return width

#neutrinos
def width_axial_nu (m_med, g_l):
    width = 3*g_l**2 * m_med / (24*math.pi)
    return width
def width_vector_nu (m_med, g_l):
    width = 3*g_l**2 * m_med / (24*math.pi)
    return width

# ======================== Total widths ============================

def totalWidthAxial(M, mDM, gq, gDM, gl) :
    return width_axial_q(M, gq) + width_axial_DM(M, mDM, gDM) + width_axial_l(M, gl) + width_axial_nu(M, gl)

def totalWidthVector(M, mDM, gq, gDM, gl) :
    return width_vector_q(M, gq) + width_vector_DM(M, mDM, gDM) + width_vector_l(M, gl) + width_vector_nu(M, gl)

# ======================== Dijet family conversions ============================

#calculate the quark coupling in the quark-only model (analysis limits) in terms of the quark coupling g' in the full model (all decay channels) 
def gqPrimeAxial (M, mDM, gq, gDM, gl):
    totalWidth = totalWidthAxial(M, mDM, gq, gDM, gl)
    gprime = math.sqrt( 1/(width_axial_q(M, 1)) * width_axial_q(M, gq)**2 / totalWidth )
    return gprime  

def gqPrimeVector(M, mDM, gq, gDM, gl):
    totalWidth = totalWidthVector(M, mDM, gq, gDM, gl)
    gprime = math.sqrt( 1/(width_vector_q(M, 1)) * width_axial_q(M, gq)**2 / totalWidth )
    return gprime
