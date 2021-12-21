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
    #x = np.complex(x,0)
    #y = np.complex(y,0)
    print(x, y)
    print(4 * x**2/y**2)
    print(1 - 4 * x**2 / y**2)
    #return np.sqrt(1 - 4 * x**2 / y**2)
    val = np.sqrt(1 - 4 * x**2 / y**2)
    print(val)
    return val

def limit_x1(x2,pid,gamma,M,mDM,ECM) :
    """
    Integration limits for x1, x2 space
    """
    # Lower limit is actually a curve
    # Upper limit is 1
    lower_lim = (4.*mDM**2)/(x2 * ECM)
    return [lower_lim, 1]


def limit_x2(pid,gamma,M,mDM,ECM) :
    """
    Integration limits for x1, x2 space
    """    
    # This is from its smallest value when x1 is largest
    # and goes up to 1.
    lower_lim = (4.*mDM**2)/ECM
    return [lower_lim, 1]

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

    # Initialise handler for lhapdfwrap if compiled
    # with lhapdf available.
    # To determine reliably from here, check for generated file.
    _wrapper = None
    if (hasLHAPDF) :
        import lhapdfwrap as pdfwrap
        _wrapper = pdfwrap.IntegrandHandler("NNPDF30_nlo_as_0118", ECM)

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
    # def propagator_monox_relative(self) :
    #     pass

    # @abc.abstractmethod
    # def hadron_level_xsec_monox_relative(self) :
    #     pass

    # @abc.abstractmethod
    # def parton_level_xsec_monox_relative(self) :
    #     pass

    # TODO:
    # Add for resonances!


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
            iwidth = 3 * self.gq**2 * yq**2 * self.mmed / \
                    (16 * PI) * beta(mq.value, self.mmed)**3

            # Only add width for mq < mmed
            width = np.where(
                mq.value < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width

    def mediator_partial_width_dm(self):
        width = self.gdm **2 * self.mmed / (8 * PI) * beta(self.mdm, self.mmed) ** 3
        return np.where(
            self.mdm < self.mmed * 0.5,
            width,
            0
        )
    
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
            iwidth = 3 * self.gq**2 * yq**2 * self.mmed / \
                    (16 * PI) * beta(mq.value, self.mmed)

            # Only add width for mq < mmed
            width = np.where(
                mq.value < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width

    def mediator_partial_width_dm(self):
        width = self.gdm **2 * self.mmed / (8 * PI) * beta(self.mdm, self.mmed)
        return np.where(
            self.mdm < self.mmed * 0.5,
            width,
            0
        )
    
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
            iwidth = self.gl**2 * self.mmed / (12*PI) * alpha(ml.value, self.mmed) * beta(ml.value, self.mmed)

            # Only add width for ml < mmed
            width = np.where(
                ml.value < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width

    # Vector
    def propagator_monox_relative(self) :
        print("Start")
        print("Mediator masses:",self.mmed)
        print("DM masses:",self.mdm
        )
        gamma = self.mediator_total_width()
        print("Got gammas")
        arctan_factor = PI/2.0 + np.arctan((self.mmed**2 - 4.*self.mdm**2)/(self.mmed*gamma))
        sigma = self.gq**2 * self.gdm**2 * arctan_factor/(self.mmed*gamma)
        return sigma

    def hadron_level_xsec_monox_relative(self) :
        return 1  


    # In case of future relevance: parton level relative xsec
    def parton_level_xsec_monox_relative(self) :

        gamma = self.mediator_total_width()
        intpoints = [self.mmed,self.mmed**2-gamma,self.mmed**2,self.mmed**2+gamma]

        # Integrate is adaptive and fundamentally
        # doesn't work with broadcasting.
        # So for this function we are going to have to 
        # actually do the values one at a time.
        xsecs = []
        if type(self.mmed) is np.ndarray or type(self.mdm) is np.ndarray :
          for mmed_i, mdm_i, gamma_i, points_i in zip(self.mmed, self.mdm, gamma, intpoints) :
            integral = integrate.quad(self._wrapper.integrand_parton_vector,4.*mdm_i**2,self.ECM,args=(gamma_i,mmed_i,mdm_i),points=points_i,limit=500)
            xsecs.append(self.gq**2 * self.gdm**2 * integral[0])
        else :
          integral = integrate.quad(self._wrapper.integrand_parton_vector,4.*self.mdm**2,self.ECM,args=(gamma,self.mmed,self.mdm),points=intpoints,limit=500)
          xsecs.append(self.gq**2 * self.gdm**2 * integral[0])
        
        return np.array(xsecs)
    
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
            iwidth = 3 * self.gq**2 * self.mmed / \
                    (12 * PI) * beta(mq.value, self.mmed)**3

            # Only add width for mq < mmed
            width = np.where(
                mq.value < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width

    def mediator_partial_width_dm(self):
        '''
        On-shell width for mediator -> DM DM.
        '''
        width = self.gdm**2 * self.mmed / (12 * PI) * beta(self.mdm, self.mmed)**3
        return np.where(
            self.mdm < self.mmed * 0.5,
            width,
            0
        )

    def mediator_partial_width_leptons(self):
        '''
        On-shell width for mediator -> l l, where l is a charged or neutral lepton.
        '''
        # Neutrinos
        width = self.gl**2 / (24*PI) * self.mmed

        # Charged leptons
        for ml in Leptons:
            iwidth = self.gl**2 * self.mmed / (12*PI) * beta(ml.value, self.mmed)**3

            # Only add width for ml < mmed
            width = np.where(
                ml.value < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width