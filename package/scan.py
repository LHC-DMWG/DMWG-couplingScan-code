from dataclasses import dataclass
import numpy as np
import abc

PI = np.pi


def beta(x, y):
    """
    Convenience function that implements part of the width formulae.
    """
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