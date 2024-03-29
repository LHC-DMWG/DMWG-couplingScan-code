from dataclasses import dataclass
import numpy as np
import abc
import math

from couplingscan.scan import *

@dataclass
class CouplingLimit_Dijet(abc.ABC) :
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
        for attr in ["mmed", "mdm", "gdm", "gl"] :
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
        if not ((self.mmed.shape == self.gq_limits.shape and len(self.gl)==1) or
                (self.mmed.shape == self.gl.shape and len(self.gq)==1)) :
            print("""Error: you must have an equal number of mediator mass and visible limit (coupling) values,
                and the other coupling to SM must be a single fixed value.""")
            print("The mediator and coupling limit points are meant to be matching x and y values. Please fix.")
            exit(1)
        
        if ('axial' not in self.coupling and 'vector' not in self.coupling) :
            print("This is only defined for axial or vector couplings!")
            exit(1)

    # This is dijet at a hadron collider: quarks in, quarks out.
    def extract_exclusion_depths(self,scan) :

        # Limit scenario is the one in which our input limit (and this class) is defined.
        # Scan scenario is the one we're going towards.

        # Extract full cross sections for target scan scenario.
        xsec_scan = scan.mediator_partial_width_quarks()**2/scan.mediator_total_width()

        # Create scan in world of input plot.
        # Need a placeholder gq around which we interpret: pick 1.
        if self.coupling == 'axial' :
            plot_world = DMAxialModelScan(
                mmed=scan.mmed,
                mdm=self.mdm,
                gq=1.0,
                gdm=self.gdm,
                gl=self.gl,
            )
        else :
            plot_world = DMVectorModelScan(
                mmed=scan.mmed,
                mdm=self.mdm,
                gq=1.0,
                gdm=self.gdm,
                gl=self.gl,
            )

        # Interpolate input gq limit curve to get all the mass points we need
        # Any points in grid that are actually above or below analysis mmed
        # values should never be excluded, so we give them a very large value        
        interpolated_limit_gq = np.interp(scan.mmed, self.mmed, self.gq_limits,left=10.,right=10.)

        # This math comes from the CMS original versions of the calculation, and works well, 
        # but is limited to cases where gdm=0 and gl=0.
        # plot_world_widths = plot_world.mediator_partial_width_quarks()
        # exclusion_depth = interpolated_limit_gq**2 * plot_world_widths / xsec_scan

        # This should be more general.
        xsec_plot_world = plot_world.mediator_partial_width_quarks()**2/plot_world.mediator_total_width()
        exclusion_depth = (xsec_plot_world/xsec_scan) * interpolated_limit_gq**2

        return exclusion_depth

# Class for scaling theory cross section and looking
# at ratio with respect to observed limit line.
@dataclass
class CrossSectionLimit1D(abc.ABC):
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
        # Check that the arrays we have been given match in shape where necessary.
        # Observed limits:
        if type(self.xsec_limit) is dict :
            for width, limit in self.xsec_limit.items() :
                if width < 0 or width > 1.0 :
                    print("Widths are interpreted as intrinsic width to mass ratio.")
                    print("The width value",width,"does not make sense in this context.")
                    print("If you have entered a percentage, please divide by 100 to convert to a fraction.")
                    exit(1)
                if (limit.shape != self.mmed_limit.shape) :
                    print("Error: limit masses and cross section values have different shapes!")
                    print("These are meant to be matching x and y values. Please fix.")
                    exit(1)
        elif (self.xsec_limit.shape != self.mmed_limit.shape) :
            print("Error: limit masses and cross section values have different shapes!")
            print("These are meant to be matching x and y values. Please fix.")
            exit(1)

        # Theory curves:
        if (self.mmed_theory.shape != self.xsec_theory.shape) :
            print("Error: theory mass points and cross section values have mismatching shapes!")
            print("These are meant to be matching x and y values. Please fix.")
            exit(1)

        # Couplings and dark matter masses: since this is a single input plot,
        # they had better all be single-valued. Make sure they're floats.
        for attr in ["mdm", "gq", "gdm", "gl"] :
            attrval = getattr(self,attr)
            if isinstance(attrval,list) or type(attrval) is np.ndarray :
                if len(attrval) > 1 :
                    print("This needs to match an input plot that is one-dimensional and has just one theory line.")
                    print("You have not provided single values for one of mdm, gq, gdm, or gl.")
                    print("This is therefore ambiguious. Please provide exactly one value for each.")
                    exit(1)
                setattr(self,attr,float(attrval[0]))
            else :
                setattr(self,attr,float(attrval))

        pass

    @abc.abstractmethod
    def get_approx_xsec(self, scan) :
        pass

    # limit_sets are the ones to select from
    # test_widths are the widths for which you need limit values
    def pick_appropriate_limit(self, test_widths, limit_widths, limit_sets) :
        # Linear interpolate between observed limits at points of interest.
        # If smaller width than smallest provided, use smallest provided.
        # If larger than largest provided, no limit can be set.
        appropriate_limits = []
        for width, limits in zip(test_widths,limit_sets.transpose()) :
            appropriate_observed = np.interp([width], limit_widths, limits,left=limits[0],right=math.nan)
            appropriate_limits.append(appropriate_observed[0])
        return np.array(appropriate_limits)        

    # This will call the inheriting methods where the cross sections differ.
    def extract_exclusion_depths(self, scan) :

        # Here we have a cross section and a (set of) observed limit(s).
        # Want to scale the cross section to the target scenarios and take
        # the ratio w.r.t. the observed limit.

        # Where more than one observed limit given, interpolate to the 
        # best value. Linear interpolation will give the most reproducible result.

        # Extract full cross sections for scan scenario.
        xsec_scan = self.get_approx_xsec(scan)

        # Interpolate theory curve and observed limit curves at the requested granularity.
        # Want interpolation to be log, not linear (that is, linear in y log axis plot).
        interp_xsec_theory = np.exp(np.interp(scan.mmed, self.mmed_theory, np.log(self.xsec_theory),left=np.nan,right=np.nan))
        # The limits may have multiple width curves.
        interp_limit = lambda mylist : np.exp(np.interp(scan.mmed, self.mmed_limit, np.log(mylist),left=np.nan,right=np.nan))
        interp_xsec_limits = np.array([interp_limit(i) for i in self.xsec_limits])

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
        xsec_plot_world = self.get_approx_xsec(plot_world)

        # Scale theory curve to the equivalent values for the scenario of interest.
        scaled_theory = interp_xsec_theory*(xsec_scan/xsec_plot_world)

        # Exclusion depth in world of plot is observed over theory
        # Want to return full set of exclusion depths for the widths given,
        # if multiple, so they can be used for other rescalings.
        compute_depth = lambda mylist : np.select([scaled_theory == 0.0, scaled_theory > 0],
                [np.nan, mylist/scaled_theory],
                default=np.nan)
        with np.errstate(divide='ignore'):
            exclusion_depths = [compute_depth(i) for i in interp_xsec_limits]

        # If we were given a width dict, return one. Otherwise just return the depths.
        if type(self.xsec_limit) is dict :
            return {k : v for k, v in zip(self.widths,exclusion_depths)}
        else : return exclusion_depths[0]

    # If we are not planning on using a rescaler
    # and instead really want the appropriate depths per width for a given scan:
    def select_depths(self,scan,depths) :
        
        widths_scan = scan.mediator_total_width()/scan.mmed
        # Was: one per mediator value being tested
        use_limits = self.pick_appropriate_limit(widths_scan,list(depths.keys()),np.array([depths[i] for i in depths.keys()]))
        return use_limits

    def get_rhs_couplinglimit(self) :
        if self.coupling == 'axial' :
            plot_world = DMAxialModelScan(
                mmed=self.mmed,
                mdm=self.mdm,
                gq=self.gq,
                gdm=self.gdm,
                gl=self.gl,
            )
        else :
            plot_world = DMVectorModelScan(
                mmed=self.mmed,
                mdm=self.mdm,
                gq=self.gq,
                gdm=self.gdm,
                gl=self.gl,
            )
        xsec_plot_world = self.get_approx_xsec(plot_world)        
        rhs_couplinglimit = xsec_plot_world*(self.xsec_limit/self.xsec_theory)
        return rhs_couplinglimit      


# For dijet, unless a different request is made of us in future,
# let's assume just one width and take everything from there.
# Make an optional intrinsic width cutoff above which the limit doesn't apply
# so user can automatically not go outside acceptable regions.
# 10% by default (very loose)
@dataclass
class CrossSectionLimit_Dijet(CrossSectionLimit1D) :

    max_intrinsic_width : float = 0.1

    def __post_init__(self) :
        super().__post_init__()

        # Now add formatting for width, any other
        # dijet specific checks
        if type(self.xsec_limit) is dict :
            print("""You've supplied a dictionary for the limits. The appropriate limit to use
            for each point will be selected based on width. For intrinsic width to mass ratios larger 
            than the largest dictionary key given, a NaN will be returned.""")
            self.widths = list(self.xsec_limit.keys())
            self.xsec_limits = np.array([self.xsec_limit[i] for i in self.widths])
        else :
            print("""You have supplied a single limit curve. This will be considered the appropriate
            limit for all signal points up to an intrinsic width to mass ratio of {0}.
            To adjust the maximum intrinsic width, please set the value of max_intrinsic_width at initialisation
            or supply a dictionary instead. For intrinsic width to mass ratios larger than this value,
            a NaN will be returned.""".format(self.max_intrinsic_width))
            self.widths = [self.max_intrinsic_width]
            self.xsec_limits = np.array([self.xsec_limit])

    # This is dijet at hadron colliders: quarks in, quarks out.
    def get_approx_xsec(self, scan) :
        return scan.mediator_partial_width_quarks()**2/scan.mediator_total_width()      

# For dilepton, different visible final state
# but also include explicit support for varying widths
# since the resolutions are so different.
@dataclass
class CrossSectionLimit_Dilepton(CrossSectionLimit1D) :

    def __post_init__(self):
        super().__post_init__()
        
        # If only a list given and not a dict for widths, print comprehensive error and quit.
        if type(self.xsec_limit) is not dict :
            print("""For dilepton limits, you need to provide input xsec limits as a dict
            with corresponding intrinsic widths as keys. This is because the limits change
            noticeably with width for dilepton signatures. If you want to ignore this issue
            and use just one limit, you still need to pick a maximum intrinsic width for
            which to consider it valid. Please try again with a dictionary.""")
            exit(1)

        self.widths = list(self.xsec_limit.keys())
        self.xsec_limits = np.array([self.xsec_limit[i] for i in self.xsec_limit.keys()])

    # This is dilepton at hadron colliders: quarks in, leptons out.
    def get_approx_xsec(self, scan) :
        return scan.mediator_partial_width_quarks()*scan.mediator_partial_width_leptons()/scan.mediator_total_width()  