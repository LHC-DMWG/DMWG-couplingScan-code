import numpy as np
from couplingscan.scan import *

# Each rescaler has a reference scan against which the others are scaled.
class Rescaler():
    # Class that houses rescaling operations
    # All methods return a multiplicative factor used to rescale signal cross sections

    def __init__(self, reference_scan) :
        self.reference_scan = reference_scan

        self.check_ref_scan()

    def check_ref_scan(self) :
        '''Need to confirm the reference scan makes sense.
        Key items: only one value of each coupling.'''
        if (self.reference_scan.gq == np.ndarray and len(self.reference_scan.gq) > 1) or \
            (self.reference_scan.gdm == np.ndarray and len(self.reference_scan.gdm) > 1) or \
            (self.reference_scan.gl == np.ndarray and len(self.reference_scan.gl) > 1) :
            print("You can only have one unique value of each coupling in your reference scan!")
            exit(1)
    
    def check_models_methods(self, method, target_model) :
        # If more models or methods introduced, modify this list and the list of
        # permitted rescalings.
        # Each method has a list of permitted groups to rescale between.
        # E.g. [['vector','axial'],['scalar','pseudoscalar']] means
        # that you can convert axial-vector to vector and vice versa,
        # and scalar to pseudoscalar and vice versa,
        # but you can't convert vector to pseudoscalar with this approach.
        methods = {
            'BR' : [['vector','axial'],['scalar','pseudoscalar']], # Pretty okay I think. Test scalar and pseudoscalar.
            # is BR same as what we were doing before for dijet? i think so, but validate.
            'propagator': [], # For now, consider this only valid within a model. Best recommended method in that case.
            'parton-level': [['vector','axial']], # Functions not yet implemented for scalar and pseudoscalar
            'hadron-level': [['vector','axial']] # Functions not yet implemented for scalar and pseudoscalar
        }

        # Commonest scenario: within a model. All good.
        if self.reference_scan._coupling == target_model : return True
        
        # Otherwise, check they exist in the same sub-list.
        options = methods[method]
        valid = False
        available = False
        for option in options :
            if target_model in option : available = True
            if self.reference_scan._coupling in option and target_model in option :
                valid = True

        if not valid :
            if not available : 
                print("Error: this rescaling method is not available for the target model!")
                print("Please choose a different method.")
                exit(1)
            else :
                print("Error: you cannot use this method to convert between",self.reference_scan._coupling,"and",target_model,"models!")
                print("Please choose a different method.")
                exit(1)                

        return

    def create_target_arrays(self,target_gq, target_gdm, target_gl) :
        '''Creates target arrays with full grid of requested couplings'''

        # Gets all combinations in 3 rows of values
        # And sets everything to floats.
        target_grid = np.array(np.meshgrid(target_gq,target_gdm,target_gl),dtype=float).reshape(3,-1)

        return target_grid

    def create_target_scan(self, target_ID, target_arrays) :

        # We already have the desired mass points from the reference scan.
        # But to make broadcasting work here we will need to increase the dimensionality
        # to include all of the target couplings as well.
        # This time we don't want meshgrid - we need to keep mass points paired up correctly

        # Repeat full set of mass points by number of tested couplings,
        # and individually repeat couplings by number of mass points.
        n_couplings = np.size(target_arrays,1)
        n_masspoints = np.size(self.reference_scan.mmed)
        target_mmed = np.tile(self.reference_scan.mmed, n_couplings)
        target_mdm = np.tile(self.reference_scan.mdm, n_couplings)
        target_couplings = np.repeat(target_arrays,n_masspoints,axis=1)

        # Now create the appropriate scan.
        if target_ID == 'axial' : 
            target_scan = DMAxialModelScan(mmed=target_mmed, mdm=target_mdm, gq=target_couplings[0],
                gdm=target_couplings[1], gl=target_couplings[2])
        elif target_ID == 'vector' :
            target_scan = DMVectorModelScan(mmed=target_mmed, mdm=target_mdm, gq=target_couplings[0],
                gdm=target_couplings[1], gl=target_couplings[2])
        elif target_ID == 'scalar' :
            target_scan = DMScalarModelScan(mmed=target_mmed, mdm=target_mdm, gq=target_couplings[0],
                gdm=target_couplings[1], gl=target_couplings[2])
        elif target_ID == 'pseudoscalar' :
            target_scan = DMPseudoModelScan(mmed=target_mmed, mdm=target_mdm, gq=target_couplings[0],
                gdm=target_couplings[1], gl=target_couplings[2])
        else :
            print("Unrecognized target model!")
            exit(1)

        return target_scan

    def format_output(self, scale_factors, target_arrays) :

        # Squish output down to a manageable format.
        # Return: {tuple of couplings : [scale factor per mass point]}
        output_dict = {}
        for i in range(np.shape(target_arrays)[1]) :
            gq, gdm, gl = target_arrays[:,i]
            output_dict[(gq, gdm, gl)] = scale_factors[i]

        return output_dict

    def rescale_by_br_quarks(self,target_gq, target_gdm, target_gl, model=None) :
        '''Rescale according to gq^2 * BR(med->DM DM). All possible
        combinations of specified couplings will be tested and
        results will be returned along with the coupling values they correspond to.'''
        
        # Check that this method of rescaling makes sense for the
        # target and reference scan types:
        if not model : model = self.reference_scan._coupling
        self.check_models_methods("BR",model)

        # Create a target scan that has the enormous dimensionality required
        # to broadcast across the full set of scanned values
        target_arrays = self.create_target_arrays(target_gq, target_gdm, target_gl)
        target_scan = self.create_target_scan(model, target_arrays)

        # Calculate scale factor at each point.
        reference_factor = self.reference_scan.mediator_partial_width_quarks() ** 2 / self.reference_scan.mediator_total_width()
        target_factors_1d = target_scan.mediator_partial_width_quarks() ** 2 / target_scan.mediator_total_width()

        # Reshape to have one row per coupling
        target_factors = np.reshape(target_factors_1d,(np.size(target_arrays,1),-1))       
        # Now this is also broadcastable
        scale_factors = target_factors / reference_factor

        return self.format_output(scale_factors,target_arrays)

    def rescale_by_br_leptons(self, target_gq, target_gdm, target_gl,model=None):
        '''Rescale according to gq^2 * BR(med->DM DM). All possible
        combinations of specified couplings will be tested and
        results will be returned along with the coupling values they correspond to.'''
        
        # Check that this method of rescaling makes sense for the
        # target and reference scan types:
        if not model : model = self.reference_scan._coupling
        self.check_models_methods("BR",model)

        # Create a target scan that has the enormous dimensionality required
        # to broadcast across the full set of scanned values
        target_arrays = self.create_target_arrays(target_gq, target_gdm, target_gl)
        target_scan = self.create_target_scan(model, target_arrays)

        # Calculate scale factor at each point.
        reference_factor = self.reference_scan.mediator_partial_width_quarks() * self.reference_scan.mediator_partial_width_leptons() / self.reference_scan.mediator_total_width()
        target_factors_1d = target_scan.mediator_partial_width_quarks() * target_scan.mediator_partial_width_leptons() / target_scan.mediator_total_width()

        # Reshape to have one row per coupling
        target_factors = np.reshape(target_factors_1d,(np.size(target_arrays,1),-1))       
        # Now this is also broadcastable        
        scale_factors = target_factors / reference_factor

        # Return nicely formatted results
        return self.format_output(scale_factors,target_arrays)

    def rescale_by_propagator(self,target_gq, target_gdm, target_gl, model=None):
        # Check that this method of rescaling makes sense for the
        # target and reference scan types:
        if not model : model = self.reference_scan._coupling
        self.check_models_methods("propagator",model)

        # Create a target scan that has the dimensionality required
        # to broadcast across the full set of scanned values
        target_arrays = self.create_target_arrays(target_gq, target_gdm, target_gl)
        target_scan = self.create_target_scan(model, target_arrays)        
        
        # Calculate scale factor at each point
        reference_factor = self.reference_scan.propagator_relative()
        target_factors_1d = target_scan.propagator_relative()

        # Reshape to have one row per coupling
        target_factors = np.reshape(target_factors_1d,(np.size(target_arrays,1),-1))       
        # Now this is also broadcastable        
        scale_factors = target_factors / reference_factor        
        
        # Return nicely formatted results
        return self.format_output(scale_factors,target_arrays)

    def rescale_by_hadronic_xsec_monox(self,target_gq, target_gdm, target_gl, model=None):
        '''Rescale using hadronic-level cross sections.'''

        # Check that this method of rescaling makes sense for the
        # target and reference scan types:
        if not model : model = self.reference_scan._coupling
        self.check_models_methods("hadron-level",model)

        # Create a target scan that has the dimensionality required
        # to broadcast across the full set of scanned values
        target_arrays = self.create_target_arrays(target_gq, target_gdm, target_gl)
        target_scan = self.create_target_scan(model, target_arrays)        
        
        for this_array in target_arrays :
            if (this_array == np.ndarray and len(this_array) > 1) :
                print("""Warning: the hadronic rescaling method takes a long time!
                We don't recommend that you use it for more than one target coupling scenario.
                Instead, try rescaling to a single target and then using the propagator scaling method
                to arrive at additional scenarios.""")

        # Calculate scale factor at each point
        reference_factor = self.reference_scan.hadron_level_xsec_monox_relative()
        target_factors_1d = target_scan.hadron_level_xsec_monox_relative()

        # Reshape to have one row per coupling
        target_factors = np.reshape(target_factors_1d,(np.size(target_arrays,1),-1))       
        # Now this is also broadcastable        
        scale_factors = target_factors / reference_factor        
        
        # Return nicely formatted results
        return self.format_output(scale_factors,target_arrays)

    def rescale_by_parton_level_xsec_monox(self,target_gq, target_gdm, target_gl, model=None):
        '''Rescale using parton-level cross sections.'''

        print('''Warning: the parton-level cross section is not the best-performing rescaling method
        in any hadron collider scenario. Consider using something else!''')

        # Check that this method of rescaling makes sense for the
        # target and reference scan types:
        if not model : model = self.reference_scan._coupling
        self.check_models_methods("parton-level",model)

        # Create a target scan that has the dimensionality required
        # to broadcast across the full set of scanned values
        target_arrays = self.create_target_arrays(target_gq, target_gdm, target_gl)
        target_scan = self.create_target_scan(model, target_arrays)        
        
        # Calculate scale factor at each point
        reference_factor = self.reference_scan.hadron_level_xsec_monox_relative()
        target_factors_1d = target_scan.hadron_level_xsec_monox_relative()

        # Reshape to have one row per coupling
        target_factors = np.reshape(target_factors_1d,(np.size(target_arrays,1),-1))       
        # Now this is also broadcastable        
        scale_factors = target_factors / reference_factor        
        
        # Return nicely formatted results
        return self.format_output(scale_factors,target_arrays)