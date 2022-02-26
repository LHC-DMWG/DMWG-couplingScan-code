from dataclasses import dataclass
import numpy as np
import abc
import math

from couplingscan.scan import *

# TODO TEMP
def interpolated_intercepts(x, y1, y2):
    """Find the intercepts of two curves, given by the same x data"""

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """    

        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0]*p2[1] - p2[0]*p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D  = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x,y

        L1 = line([point1[0],point1[1]], [point2[0],point2[1]])
        L2 = line([point3[0],point3[1]], [point4[0],point4[1]])

        R = intersection(L1, L2)

        return R

    idxs = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)

    xcs = []
    ycs = []

    for idx in idxs:
        xc, yc = intercept((x[idx], y1[idx]),((x[idx+1], y1[idx+1])), ((x[idx], y2[idx])), ((x[idx+1], y2[idx+1])))
        xcs.append(xc[0])
        ycs.append(yc[0])
    return np.array([xcs,ycs])
# TODO END TEMP

@dataclass
class DMVisibleLimit1D(abc.ABC):
    mmed: float
    gq_limits: float
    mdm: float
    gdm: float
    gl: float
    coupling : str

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
        
        if ('axial' not in self.coupling and 'vector' not in self.coupling) :
            print("This is only defined for axial or vector couplings!")
            exit(1)            

    @abc.abstractmethod
    def extract_exclusion_depths(self) :
        pass


@dataclass
class CouplingLimit_Dijet(DMVisibleLimit1D) :

    def __post_init__(self):
        gq: -1

    # This is dijet at a hadron collider: quarks in, quarks out.
    def extract_exclusion_depths(self,scan) :

        # Limit scenario is the one in which our input limit (and this class) is defined.
        # Scan scenario is the one we're going towards.

        # Extract full cross sections for scan scenario.
        xsec_scan = scan.mediator_partial_width_quarks()**2/scan.mediator_total_width()

        # These are limit gq values in the plot for each mediator mass of interest
        # Any points in grid that are actually above or below analysis mmed
        # values should never be excluded, so we give them a very large value.
        interpolated_limit_gq = np.interp(scan.mmed,self.mmed, self.gq_limits,left=10.,right=10.)

        # These are gq': the equivalent gq in the world of the plot
        if self.coupling == 'axial' :
            plot_world = DMAxialModelScan(
                mmed=scan.mmed,
                mdm=self.mdm,
                gq=1.0,
                gdm=0.0,
                gl=0.0,
            )
        else :
            plot_world = DMVectorModelScan(
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
# Instead, scaling theory cross section and looking
# at ratio with respect to observed limit line.
@dataclass
class CrossSectionLimit_Dilepton(abc.ABC) :

    mmed_limit : float
    xsec_limit : float
    mmed_theory : float
    xsec_theory : float
    mdm: float
    gq: float
    gdm: float
    gl: float
    coupling : str

    def __post_init__(self):
        
        # Add some checks: if we have only one limit,
        # do we assign it a default width or do we toss it 
        # back and say "don't do this with dilepton"?

        # TODO set this up so that everything is turned into the format we expect
        # if it isn't already.
        # If only a list given and not a dict for widths, print comprehensive error and quit.

        self.widths = list(self.xsec_limit.keys())
        self.xsec_limits = np.array([self.xsec_limit[i] for i in self.xsec_limit.keys()])

    def pick_appropriate_limit(self, test_widths, granular_limits) :
        # Linear interpolate between observed limits at points of interest.
        # If smaller width than smallest provided, use smallest provided.
        # If larger than largest provided, no limit can be set.
        appropriate_limits = []
        for width, limits in zip(test_widths,granular_limits.transpose()) :
            appropriate_observed = np.interp([width], self.widths, limits,left=limits[0],right=math.nan)
            appropriate_limits.append(appropriate_observed[0])
        return np.array(appropriate_limits)

    # This is dilepton at a hadron collider: quarks in, leptons out.
    def extract_exclusion_depths(self,scan) :

        # Here we have a cross section and a (set of) observed limit(s).
        # Want to scale the cross section to the target scenarios and take
        # the ratio w.r.t. the observed limit.

        # Where more than one observed limit given, interpolate to the 
        # best value. Linear interpolation will give the most reproducible result.

        # Extract full cross sections for scan scenario.
        xsec_scan = scan.mediator_partial_width_quarks()*scan.mediator_partial_width_leptons()/scan.mediator_total_width()

        # Calculate full widths also: we'll need this.
        widths_scan = scan.mediator_total_width()/scan.mmed

        # Interpolate theory curve and observed limit curves at the requested granularity.
        interp_xsec_theory = np.interp(scan.mmed, self.mmed_theory, self.xsec_theory,left=0,right=0)
        # The limits may have multiple width curves.
        interp_limit = lambda mylist : np.interp(scan.mmed, self.mmed_limit, mylist,left=1e3,right=1e3)
        interp_xsec_limit = np.array([interp_limit(i) for i in self.xsec_limits])

        # Get equivalent theory cross sections at the desired mass points in the world of the input xsec limit plot.
        if self.coupling == 'axial' :
            plot_world = DMAxialModelScan(
                mmed=scan.mmed,
                mdm=self.mdm,
                gq=self.gq,
                gdm=self.gdm,
                gl=self.gl,
            )
        else :
            plot_world = DMVectorModelScan(
                mmed=scan.mmed,
                mdm=self.mdm,
                gq=self.gq,
                gdm=self.gdm,
                gl=self.gl,
            )
        xsec_plot_world = plot_world.mediator_partial_width_quarks()*plot_world.mediator_partial_width_leptons()/plot_world.mediator_total_width()

        # Scale theory curve to the equivalent values for the scenario of interest.
        scaled_theory = interp_xsec_theory*(xsec_scan/xsec_plot_world)

        # Exclusion depth in world of plot is observed over theory
        # Which observed line? This depends on the intrinsic width.
        # Obtain and calculate the ratio.
        relevant_limits = self.pick_appropriate_limit(widths_scan,interp_xsec_limit)

        # TODO: Clear before publishing
        # Want to take slices through this and check curves. Why so wobbly?
        # Pattern: mmed goes up while mdm stays fixed, then starts again, etc
        # test_strips = np.reshape(relevant_limits,(-1,171))
        # theory_lines = np.reshape(scaled_theory,(-1,171))
        # dmvals = np.reshape(scan.mdm,(-1,171))
        # these_xvals = scan.mmed[:171]
        # import matplotlib.pyplot as plt
        # myarray = np.array([[0,0]])
        # for strip, line, val in zip(test_strips, theory_lines, dmvals[:,0]) :
        #     plt.clf()
        #     plt.plot(these_xvals,strip,label="observed")
        #     plt.plot(these_xvals,line,label="theory")
        #     plt.xlim(400, 1900)
        #     plt.legend()
        #     plt.savefig('plots/validation/dilepton_mdm{0}.pdf'.format(int(val)),bbox_inches='tight')

        #     # Second validation: where is curve if I do it exactly Etienne's way,
        #     # with the crossing points?
        #     edge_points = interpolated_intercepts(these_xvals, strip, line)
        #     for xval, yval in zip(edge_points[0],edge_points[1]) :
        #         myarray = np.vstack([myarray,[xval,val]])
        # sorted_edges = myarray[myarray[:,0].argsort()]
        # plt.clf()
        # plt.plot(sorted_edges[:,0],sorted_edges[:,1])
        # plt.legend()
        # plt.savefig('plots/validation/etienne_method.pdf',bbox_inches='tight')

        # Don't want to divide by zero in case we have impossible theory points,
        # so let's do this piecewise and set any points with zero cross section to an
        # arbitrary very-unexcluded value
        exclusion_depth = np.select([scaled_theory == 0.0, scaled_theory > 0],
                [100, relevant_limits/scaled_theory],
                default=np.nan)

        return exclusion_depth    