from couplingscan.scan import *
from couplingscan.rescaler import *

# Example one: mostly constant parameters, scan in 1d
# numpy broadcasting does its magic

scan1 = DMScalarModelScan(
mmed=1000,
mdm=np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
gq=0.25,
gdm=1.0,
gl=0.0,
)
print("Scan 1 test:")
print(scan1.mediator_partial_width_dm() / scan1.mediator_total_width())

# Example two: multiple variables scanned
# numpy broadcasting still handles this
scan2 = DMScalarModelScan(
mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
mdm=np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
gq=0.25,
gdm=1.0,
gl=0.0,
)
print("Scan 2 test:")
print(scan2.mediator_partial_width_dm() / scan2.mediator_total_width())

# Example 2.5: make sure the pseudoscalar
# model isn't broken
scan2 = DMPseudoModelScan(
mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
mdm=np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
gq=0.25,
gdm=1.0,
gl=0.0,
)
print("Scan 2.5 test:")
print(scan2.mediator_partial_width_dm() / scan2.mediator_total_width())

# Example three: propagators, arrays, vector model
scan3 = DMVectorModelScan(mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
mdm=np.array([1,10,50,100,150,200,250,300,350,400,450], dtype=float),
gq=0.25,
gdm=1.0,
gl=0.0,
)
print("Scan 3, propagator:")
print(scan3.propagator_relative())

#Example 3.2: parton-level
print("Parton level:")
print(scan3.parton_level_xsec_monox_relative())

#Example 3.3: hadron-level
print("Hadron level:")
print(scan3.hadron_level_xsec_monox_relative())

# Example four: same tests but for axial-vector
scan4 = DMAxialModelScan(mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450]),
mdm=np.array([1,10,50,100,150,200,250,300,350,400,450]),
gq=0.25,
gdm=1.0,
gl=0.0,
)
print("Scan 4, propagator:")
print(scan4.propagator_relative())

# Example 4.2: parton-level
print("Parton level:")
print(scan4.parton_level_xsec_monox_relative())

# Example 4.3: hadron-level
# Slow, so comment out if you want faster tests.
print("Hadron level:")
print(scan4.hadron_level_xsec_monox_relative())

# Now let's try a rescaler.
rescaleA1 = Rescaler(scan4)
# Put some target couplings and target model,
# and pick method to use to get scale factors.
scalefactors_A2 = rescaleA1.rescale_by_br_quarks(target_gq=0.1,target_gdm=1,target_gl=0.01,model='axial')
print("Simple scale factors:")

# And let's try a more complicated one.
scalefactors_several = rescaleA1.rescale_by_br_quarks(target_gq=[0.25, 0.2],target_gdm=1,target_gl=[0.0, 0.05, 0.1],model='axial')
print("Mulitdimensional scale factors:")
print(scalefactors_several)
