import numpy as np
from couplingscan.scan import *

@dataclass
class CouplingLimiter(abc.ABC) :

    def __init__(self) :
        pass

@dataclass
class Dijet_CouplingLimiter(CouplingLimiter) :
    '''
    Extracts limits on couplings from a scan within a dijet scenario
    '''
    
    def extract_limits_gq(self, scan, exclusiondepths) :

        