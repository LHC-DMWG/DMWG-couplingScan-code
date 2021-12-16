from dataclasses import dataclass
from enum import Enum
import numpy as np
import abc

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
    return np.sqrt(1 - 4 * x**2 / y**2)

def limit_x1(x2,pid,gamma,M,mDM) :
    """
    Integration limits for x1, x2 space
    """
    # Lower limit is actually a curve
    # Upper limit is 1
    lower_lim = (4.*mDM**2)/(x2 * ECM)
    return [lower_lim, 1]


def limit_x2(pid,gamma,M,mDM) :
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

    # Default CME = LHC current
    self.ECM = 13000.**2
    self.PDF = "NNPDF30_nlo_as_0118"
    self.wrapper = pdfwrap.IntegrandHandler(self.PDF, ECM)

    self.import lhapdfwrap as pdfwrap

    @abc.abstractmethod
    def mediator_total_width(self):
        pass

    @abc.abstractmethod
    def mediator_partial_width_quarks(self):
        pass

    @abc.abstractmethod
    def mediator_partial_width_dm(self):
        pass

    # and so on for more partial width functions,
    # propagator terms, whatever ...

    def propagator(self):
        pass

@dataclass
class DMScalarModelScan(DMModelScan):
    '''
    Specific implementation of a parameter scan
    for a scalar mediator.
    '''
    _coupling: str = 'scalar'

    def mediator_total_width(self):
        # TODO
        pass

    def mediator_partial_width_quarks(self):
        # TODO
        pass

    def mediator_partial_width_dm(self):
        gamma = self.gdm **2 * self.mmed / (8 * PI) * beta(self.mdm, self.mmed) ** 3
        return np.where(
            self.mdm < self.mmed * 0.5,
            gamma,
            0
        )


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
            iwidth = 3 * self.gq**2 * self.mmed / \
                    (12 * PI) * alpha(mq, self.mmed) * beta(mq, self.mmed)

            # Only add width for mq < mmed
            width = np.where(
                mq < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width

    def mediator_partial_width_dm(self):
        '''
        On-shell width for mediator -> DM DM.
        '''
        width = self.gdm**2 * self.mmed / (12 * PI) * alpha(self.mdm, self.mmed) * beta(self.mdm, self.mmed)
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
            iwidth = self.gl**2 * self.mmed / (12*PI) * alpha(ml, self.mmed) / beta(ml, self.mmed)

            # Only add width for ml < mmed
            width = np.where(
                ml < self.mmed * 0.5,
                width + iwidth,
                width
            )


        return width


@dataclass
class DMAxialVectorModelScan(DMModelScan):
    '''
    Specific implementation of a parameter scan
    for an axial vector mediator.
    '''
    _coupling: str = 'axialvector'

    def mediator_total_width(self):
        return self.mediator_partial_width_quarks() + self.mediator_partial_width_dm() + self.mediator_partial_width_leptons()

    def mediator_partial_width_quarks(self):
        '''
        On-shell width for mediator -> q q.
        '''
        width = 0
        for mq in Quarks:
            iwidth = 3 * self.gq**2 * self.mmed / (12 * PI) * beta(mq, self.mmed)**3
            # Only add width for mq < mmed
            width = np.where(
                mq < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width        

    def mediator_partial_width_dm(self):
        '''
        On-shell width for mediator -> DM DM.
        '''    
        width = (self.gdm**2 * self.mmed)/(12.*PI) * beta(self.mdm, self.mmed)**3
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
            iwidth = self.gl**2 * self.mmed / (12*PI) * alpha(ml, self.mmed) / beta(ml, self.mmed)**3

            # Only add width for ml < mmed
            width = np.where(
                ml < self.mmed * 0.5,
                width + iwidth,
                width
            )


        return width        