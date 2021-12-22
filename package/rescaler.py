import numpy as np

# Each rescaler has a reference scan against which the others are scaled.
class Rescaler():
    # Class that houses rescaling operations
    # All methods return a multiplicative factor used to rescale signal cross sections
    # TODO: decide if we should include the observed limits in here, thereby
    # returning the final scaled values,
    # or not include them and just return scale factors.

    def __init__(self, reference_scan) :
        self.reference_scan = reference_scan

        self.check_ref_scan()

    def check_ref_scan(self) :
        '''Need to confirm the reference scan makes sense.
        Key items: only one value of each coupling.'''
        if (self.reference_scan.gq is np.ndarray and len(self.reference_scan.gq) > 1) or \
            (self.reference_scan.gdm is np.ndarray and len(self.reference_scan.gdm) > 1) or \
            (self.reference_scan.gl is np.ndarray and len(self.reference_scan.gl) > 1) :
            print("You can only have one unique value of each coupling in your reference scan!")
            exit(1)

    def check_targets(self,target_gq, target_gdm, target_gl) :
        pass

    def rescale_by_br(self,target_gq, target_gdm, target_gl):
        '''Rescale according to gq^2 * BR(med->DM DM)'''
        
        # check that the two scans actually make sense
        # i.e. same variables being scanned, etc
        self.consistency_check(reference_scan, target_scan)

        # Create a target scan that has the enormous dimensionality required
        # to broadcast across the full set of scanned values


        # Squish output down to a manageable format


        # Return masses (?) and scale factors for data points

        reference_factor = reference_scan.br_mediator_dm() * reference_scan.gq ** 2
        target_factor = target_scan.br_mediator_dm() * target_scan.gq ** 2

        return target_factor / reference_factor

    def rescale_by_another_method(self,reference_scan, target_scan):
        pass
