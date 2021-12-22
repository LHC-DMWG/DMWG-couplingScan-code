from dataclasses import dataclass
from enum import Enum
import numpy as np
import abc
import imp
import scipy.integrate as integrate

# Check if lhapdf was available at compile time. 
try:
    imp.find_module('package')
    hasLHAPDF = True
except:
    hasLHAPDF = False

PI = np.pi

class Quarks(Enum):
    up=0.0024
    down=0.0048
    strange=0.104
    charm=1.27
    bottom=4.2
    top=171.2

class Leptons(Enum):
    electron=0.000511
    muon=0.105658
    tau=1.77682


def alpha(x, y):
    """
    Convenience function that implements part of the width formulae.
    """
    return 1 + 2 * x**2 / y**2

def beta(x, y):
    """
    Convenience function that implements part of the width formulae.
    """
    # Not sure why these are here?
    return np.sqrt(1 - 4 * x**2 / y**2)

@dataclass
class DMModelScan(abc.ABC):
    '''
    Abstract parent class for scans of different model types.
    '''
    mmed: float
    mdm: float
    gq: float
    gdm: float
    gl: float
    _coupling: str

    ECM: float = 13000.**2

    # Note: only doing up and down PDFs because I tested
    # with all 4 light quarks and saw no discernable difference.
    # If needed, increase here.
    _nquarks_pdf: int = 2

    # Initialise handler for lhapdfwrap if compiled
    # with lhapdf available.
    # To determine reliably from here, check for generated file.
    _wrapper = None
    if (hasLHAPDF) :
        import lhapdfwrap as pdfwrap
        _wrapper = pdfwrap.IntegrandHandler("NNPDF30_nlo_as_0118", ECM)

    def __post_init__(self):

        # Various safety controls:
        # If any starting parameter is just a float, make it into a 1-item array.
        # For all the others, make sure they have type float.
        for attr in ["mmed", "mdm", "gq", "gdm", "gl"] :
            attrval = getattr(self,attr)
            if type(attrval) is not np.ndarray :
                setattr(self,attr,np.array([attrval],dtype=float))
            else :
                setattr(self,attr,attrval.astype(float))

        # Check that the arrays we have been given match in shape where necessary.
        if self.mmed.shape != self.mdm.shape :
            print("Error: mass points have mismatching shapes!")
            print("These are meant to be matching x and y values. Please fix.")
            exit(1)
        if not (self.gq.shape == self.gdm.shape == self.gl.shape) :
            print("Error: coupling points have mismatching shapes!")
            print("Each point in your scan must have exactly one gq, gdm, and gl.")
            exit(1)

    @abc.abstractmethod
    def mediator_total_width(self):
        pass

    @abc.abstractmethod
    def mediator_partial_width_quarks(self):
        pass

    @abc.abstractmethod
    def mediator_partial_width_dm(self):
        pass

    # "Relative" in function names from here on
    # indicates that these are not full cross sections
    # but do correctly define ratios of cross sections
    # calculated using the same functions.

    # May want to move these since some are signature dependent
    # I.e. are there equivalents to these for scalar/pseudoscalar?
    # @abc.abstractmethod
    # def propagator_relative(self) :
    #     pass

    # @abc.abstractmethod
    # def hadron_level_xsec_monox_relative(self) :
    #     pass

    # @abc.abstractmethod
    # def parton_level_xsec_monox_relative(self) :
    #     pass

    # TODO:
    # Add for resonances!

    # Next four functions used by both vector and axial-vector, so put them here

    # Make it check the various points in S 
    # that allow the integral to converge correctly.
    # Add the point at which the integrand goes to zero
    # since that's a discontinuity now.
    def opts_x1(self,x2,pid,gamma,M,mDM) :
        points_list = [M**2/(x2 * self.ECM)]  
        # Don't both returning ridge location if we're off-shell
        if (M < 2.*mDM) :
            return {}
        else :
            return {'points' : points_list}

    # Nothing special.
    def opts_x2(self,pid,gamma,M,mDM) :
        return {}

    def limit_x1(self,x2,pid,gamma,M,mDM) :
        """
        Integration limits for x1, x2 space
        """
        # Lower limit is actually a curve
        # Upper limit is 1
        lower_lim = (4.*mDM**2)/(x2 * self.ECM)
        return [lower_lim, 1]

    def limit_x2(self,pid,gamma,M,mDM) :
        """
        Integration limits for x1, x2 space
        """    
        # This is from its smallest value when x1 is largest
        # and goes up to 1.
        lower_lim = (4.*mDM**2)/self.ECM
        return [lower_lim, 1]        


@dataclass
class DMScalarModelScan(DMModelScan):
    '''
    Specific implementation of a parameter scan
    for a scalar mediator.
    '''
    _coupling: str = 'scalar'

    def mediator_total_width(self):
        return self.mediator_partial_width_quarks() + self.mediator_partial_width_dm() + self.mediator_partial_width_gluon()
    
    def mediator_partial_width_quarks(self):
        width = 0
        v = 246
        for mq in Quarks:
            yq = np.sqrt(2) * mq.value / v
            iwidth = np.piecewise(
                self.mmed,
                mq.value < self.mmed * 0.5,
                [lambda x : 3 * self.gq**2 * yq**2 * x / (16 * PI) * beta(mq.value, x)**3, 0]
            )

            width += iwidth

        return width

    def mediator_partial_width_dm(self):
        width = np.piecewise(
            self.mmed,
            self.mdm < self.mmed * 0.5,
            [lambda x : self.gdm **2 * x / (8 * PI) * beta(self.mdm, x) ** 3, 0]

        )
        return width
    
    def mediator_partial_width_gluon(self):
        alphas = 0.130
        v = 246
        width = alphas ** 2 * self.gq**2 * self.mmed**3 / (32 * PI**3 * v**2)
        width = width * np.abs(self.fs(4 * (Quarks.top.value / self.mmed)**2))**2
        return width

    def fs(self,tau):
        tau = np.complex(tau,0)
        return tau * (1 + (1 - tau) * (np.arctan(1. / np.sqrt(tau - 1)))**2)
        
@dataclass
class DMPseudoModelScan(DMModelScan):
    '''
    Specific implementation of a parameter scan
    for a pseudoscalar mediator.
    '''
    _coupling: str = 'pseudo'

    def mediator_total_width(self):
        return self.mediator_partial_width_quarks() + self.mediator_partial_width_dm() + self.mediator_partial_width_gluon()
    
    def mediator_partial_width_quarks(self):
        width = 0
        v = 246
        for mq in Quarks:
            yq = np.sqrt(2) * mq.value / v
            iwidth = np.piecewise(
                self.mmed,
                mq.value < self.mmed * 0.5,
                [lambda x : 3 * self.gq**2 * yq**2 * x / (16 * PI) * beta(mq.value, x), 0]
            )
            width += iwidth

        return width

    def mediator_partial_width_dm(self):
        width = np.piecewise(
            self.mmed,
            self.mdm < self.mmed * 0.5,
            [lambda x : self.gdm **2 * x / (8 * PI) * beta(self.mdm, x), 0]
        )
        
        return width
    
    def mediator_partial_width_gluon(self):
        alphas = 0.130
        v = 246
        width = alphas ** 2 * self.gq**2 * self.mmed**3 / (32 * PI**3 * v**2)
        width = width * np.abs(self.fps(4 * (Quarks.top.value / self.mmed)**2))**2
        return width

    def fps(self,tau):
        tau = np.complex(tau,0)
        return tau * (np.arctan(1. / np.sqrt(tau - 1)))**2

@dataclass
class DMVectorModelScan(DMModelScan):
    '''
    Specific implementation of a parameter scan
    for a vector mediator.
    '''
    _coupling: str = 'vector'

    def mediator_total_width(self):
        return self.mediator_partial_width_quarks() + self.mediator_partial_width_dm() + self.mediator_partial_width_leptons()

    def mediator_partial_width_quarks(self):
        '''
        On-shell width for mediator -> q q.
        '''

        width = 0
        for mq in Quarks:
            # Only calculate when mq < 0.5 mmed, or we get errors.
            iwidth = np.piecewise(
                self.mmed,
                [mq.value < self.mmed * 0.5],
                [lambda x : 3 * self.gq**2 * x / (12 * PI) * alpha(mq.value, x) * beta(mq.value, x), 0]
            )
            width += iwidth

        return width

    def mediator_partial_width_dm(self):
        '''
        On-shell width for mediator -> DM DM.
        '''
        width = np.piecewise(
            self.mmed,
            [self.mdm < self.mmed * 0.5],
            [lambda x : self.gdm**2 * x / (12 * PI) * alpha(self.mdm, x) * beta(self.mdm, x), 0]
        )
        return width

    def mediator_partial_width_leptons(self):
        '''
        On-shell width for mediator -> l l, where l is a charged or neutral lepton.
        '''
        # Neutrinos
        width = self.gl**2 / (24*PI) * self.mmed

        # Charged leptons
        for ml in Leptons:
            # Only add width for ml < mmed
            iwidth = np.piecewise(
                self.mmed,
                [ml.value < self.mmed * 0.5],
                [lambda x : self.gl**2 * x / (12*PI) * alpha(ml.value, x) * beta(ml.value, x),0]
            )
            
            width += iwidth

        return width

    def propagator_relative(self) :
        '''
        Integral of full propagator expression for vector mediator
        '''        
        gamma = self.mediator_total_width()
        arctan_factor = PI/2.0 + np.arctan((self.mmed**2 - 4.*self.mdm**2)/(self.mmed*gamma))
        sigma = self.gq**2 * self.gdm**2 * arctan_factor/(self.mmed*gamma)
        return sigma

    def hadron_level_xsec_monox_relative(self) :
        '''
        (Relative) hadron-level cross section for vector mediator to DM
        '''        
        gamma = self.mediator_total_width()
        xsecs = []
        for mmed_i, mdm_i, gamma_i in zip(self.mmed, self.mdm, gamma) :
            xsec = 0
            for q_pid in range(1,self._nquarks_pdf) :  
                integral = integrate.nquad(self._wrapper.integrand_hadronic_vector,[self.limit_x1,self.limit_x2],args=(q_pid,gamma_i,mmed_i,mdm_i),opts=[self.opts_x1,self.opts_x2])
                xsec = xsec + integral[0]
            xsecs.append(xsec)
        # For properly broadcasting gq and gdm dependence
        xsecs = self.gq**2 * self.gdm**2 * xsecs
        return xsecs

    # In case of future relevance: parton level relative xsec
    def parton_level_xsec_monox_relative(self) :
        '''
        (Relative) parton-level cross section for vector mediator to DM
        ''' 
        gamma = self.mediator_total_width()

        # Integrate is adaptive and fundamentally
        # doesn't work with broadcasting.
        # So for this function we are going to have to 
        # actually do the values one at a time.
        # if type(self.mmed) is np.ndarray or type(self.mdm) is np.ndarray :
        xsecs = []
        for mmed_i, mdm_i, gamma_i in zip(self.mmed, self.mdm, gamma) :
            intpoints = [mmed_i,mmed_i**2-gamma_i,mmed_i**2,mmed_i**2+gamma_i]
            integral = integrate.quad(self._wrapper.integrand_parton_vector,4.*mdm_i**2,self.ECM,args=(gamma_i,mmed_i,mdm_i),points=intpoints,limit=500)
            xsecs.append(integral[0])
        xsecs = self.gq**2 * self.gdm**2 * xsecs
        return xsecs 

            
@dataclass
class DMAxialModelScan(DMModelScan):
    '''
    Specific implementation of a parameter scan
    for an axial vector mediator.
    '''
    _coupling: str = 'axial'

    def mediator_total_width(self):
        return self.mediator_partial_width_quarks() + self.mediator_partial_width_dm() + self.mediator_partial_width_leptons()

    def mediator_partial_width_quarks(self):
        '''
        On-shell width for mediator -> q q.
        '''
        width = 0
        for mq in Quarks:
            # Only add width for mq < mmed
            iwidth = np.piecewise(
                self.mmed,
                [mq.value < self.mmed * 0.5],     
                [lambda x : 3 * self.gq**2 * x / (12 * PI) * beta(mq.value, x)**3, 0]
            )

            width += iwidth

        return width

    def mediator_partial_width_dm(self):
        '''
        On-shell width for mediator -> DM DM.
        '''
        width = np.piecewise(
            self.mmed,
            [self.mdm < self.mmed * 0.5],        
            [lambda x : self.gdm**2 * x / (12 * PI) * beta(self.mdm, x)**3, 0]
        )

        return width

    def mediator_partial_width_leptons(self):
        '''
        On-shell width for mediator -> l l, where l is a charged or neutral lepton.
        '''
        # Neutrinos
        width = self.gl**2 / (24*PI) * self.mmed

        # Charged leptons
        for ml in Leptons:
           # Only add width for ml < mmed
            iwidth = np.piecewise(
                self.mmed,
                [ml.value < self.mmed * 0.5],     
                [lambda x : self.gl**2 * x / (12*PI) * beta(ml.value, x)**3, 0]
            )
            
            width += iwidth

        return width

    def propagator_relative(self) :
        '''
        Integral of full propagator expression for axial-vector mediator
        '''        
        gamma = self.mediator_total_width()
        arctan_factor = PI/2.0 + np.arctan((self.mmed**2 - 4.*self.mdm**2)/(self.mmed*gamma))
        sigma = self.gq**2 * self.gdm**2 * arctan_factor/(self.mmed*gamma)
        return sigma       

    def hadron_level_xsec_monox_relative(self) :
        '''
        (Relative) hadron-level cross section for axial-vector mediator to DM
        '''        
        gamma = self.mediator_total_width()
        xsecs = []
        for mmed_i, mdm_i, gamma_i in zip(self.mmed, self.mdm, gamma) :
            xsec = 0
            for q_pid in range(1,self._nquarks_pdf) :  
                integral = integrate.nquad(self._wrapper.integrand_hadronic_axialvector,[self.limit_x1,self.limit_x2],args=(q_pid,gamma_i,mmed_i,mdm_i),opts=[self.opts_x1,self.opts_x2],full_output=True) 
                xsec = xsec + integral[0]
            xsecs.append(xsec)
        xsecs = self.gq**2 * self.gdm**2 * xsecs
        return xsecs

    # In case of future relevance: parton level relative xsec
    def parton_level_xsec_monox_relative(self) :
        '''
        (Relative) parton-level cross section for axial-vector mediator to DM
        ''' 
        gamma = self.mediator_total_width()

        # Integrate is adaptive and fundamentally
        # doesn't work with broadcasting.
        # So for this function we are going to have to 
        # actually do the values one at a time.
        xsecs = []
        for mmed_i, mdm_i, gamma_i in zip(self.mmed, self.mdm, gamma) :
            intpoints = [mmed_i,mmed_i**2-gamma_i,mmed_i**2,mmed_i**2+gamma_i]
            integral = integrate.quad(self._wrapper.integrand_parton_axialvector,4.*mdm_i**2,self.ECM,args=(gamma_i,mmed_i,mdm_i),points=intpoints,limit=500)
            xsecs.append(integral[0])
        xsecs = self.gq**2 * self.gdm**2 * xsecs
        return xsecs  