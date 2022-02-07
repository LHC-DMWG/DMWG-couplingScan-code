from dataclasses import dataclass
import numpy as np
import abc
import math

from couplingscan.scan import *

@dataclass
class DMVisibleLimit1D(abc.ABC):
    mmed: float
    gq_limits: float
    mdm: float
    gdm: float
    gl: float
    #gq: float

    def __post_init__(self):

        # Various safety controls:
        # If any starting parameter is just a float, make it into a 1-item array.
        # For all the others, make sure they have type float.
        for attr in ["mmed", "mdm", "gq", "gdm", "gl"] :
            attrval = getattr(self,attr)
            if isinstance(attrval,list) :
                setattr(self,attr,np.array(attrval,dtype=float))
            elif type(attrval) is not np.ndarray :
                setattr(self,attr,np.array([attrval],dtype=float))
            else :
                setattr(self,attr,attrval.astype(float))

        # Check that the arrays we have been given have the shapes we expect.
        # For the 1d visible limit like this, we expect mmed and one of the visible couplings to match.
        # The other has to have a single value, and there can be only one value of mdm.
        if (len(self.mdm) > 1) :
            print("Error: there should be a fixed DM mass for this type of limit!")
            print("If you are treating DM as decoupled, you can just set that value very high.")
            exit(1)
        if not ((self.mmed.shape == self.gq_limit.shape and len(self.gl)==1) or
                (self.mmed.shape == self.gl.shape and len(self.gq)==1)) :
            print("""Error: you must have an equal number of mediator mass and visible limit (coupling) values,
                and the other coupling to SM must be a single fixed value.""")
            print("The mediator and coupling limit points are meant to be matching x and y values. Please fix.")
            exit(1)

    @abc.abstractmethod
    def extract_exclusion_depths(self) :
        pass


@dataclass
class DMLimit1D_AxialDijet(DMVisibleLimit1D) :

    def __post_init__(self):
        gq: -1

    # This is dijet at a hadron collider: quarks in, quarks out.
    def extract_exclusion_depths(self,scan) :

        # Limit scenario is the one in which our input limit (and this class) is defined.
        # Scan scenario is the one we're going towards.

        # Extract full cross sections for scan scenario.
        xsec_scan = scan.mediator_partial_width_quarks()**2/scan.mediator_total_width()

        # These are limit gq values in the plot for each mediator mass of interest
        interpolated_limit_gq = np.interp(scan.mmed,self.mmed, self.gq_limits)

        # These are gq': the equivalent gq in the world of the plot
        plot_world = DMAxialModelScan(
            mmed=scan.mmed,
            mdm=self.mdm,
            gq=1.0,
            gdm=0.0,
            gl=0.0,
        )
        plot_world_widths = plot_world.mediator_partial_width_quarks()

        exclusion_depth = interpolated_limit_gq**2 * plot_world_widths / xsec_scan
        return exclusion_depth


# For dilepton, can't use the same math because it doesn't cancel out!
# Should we do everything Etienne's way instead?
# Compare dijet done with Phil's method to dijet done with Etienne's method. Different?