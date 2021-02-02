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
    x = np.complex(x,0)
    y = np.complex(y,0)
    return np.sqrt(1 - 4 * x**2 / y**2)

@dataclass
class DMScalarModelScan:
    '''
    Specific implementation of a parameter scan
    for a scalar mediator.
    '''
    _coupling: str = 'scalar'
    mmed: int = 0
    mdm: int = 0
    gq: float = 0.0
    gdm: float = 0.0

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
class DMPseudoModelScan:
    '''
    Specific implementation of a parameter scan
    for a scalar mediator.
    '''
    _coupling: str = 'pseudo'
    mmed: int = 0
    mdm: int = 0
    gq: float = 0.0
    gdm: float = 0.0

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
class DMVectorModelScan:
    '''
    Specific implementation of a parameter scan
    for a vector mediator.
    '''
    _coupling: str = 'vector'
    mmed: int = 0
    mdm: int = 0
    gq: float = 0.0
    gdm: float = 0.0
    gl: float = 0.0

    def mediator_total_width(self):
        return self.mediator_partial_width_quarks() + self.mediator_partial_width_dm() + self.mediator_partial_width_leptons()

    def mediator_partial_width_quarks(self):
        '''
        On-shell width for mediator -> q q.
        '''
        width = 0
        for mq in Quarks:
            iwidth = 3 * self.gq**2 * self.mmed / \
                    (12 * PI) * alpha(mq.value, self.mmed) * beta(mq.value, self.mmed)

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
            iwidth = self.gl**2 * self.mmed / (12*PI) * alpha(ml.value, self.mmed) * beta(ml.value, self.mmed)

            # Only add width for ml < mmed
            width = np.where(
                ml.value < self.mmed * 0.5,
                width + iwidth,
                width
            )
        return width
    
@dataclass
class DMAxialModelScan:
    '''
    Specific implementation of a parameter scan
    for a vector mediator.
    '''
    _coupling: str = 'axial'
    mmed: int = 0
    mdm: int = 0
    gq: float = 0.0
    gdm: float = 0.0
    gl: float = 0.0
    

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
