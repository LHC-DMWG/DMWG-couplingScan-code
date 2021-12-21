from package.scan import *

# Example one: mostly constant parameters, scan in 1d
# numpy broadcasting does its magic
# scan1 = DMScalarModelScan(
# mmed=1000,
# mdm=np.array([1,10,50,100,150,200,250,300,350,400,450]),
# gq=0.25,
# gdm=1.0,
# gl=0.0,
# )
# print("Scan 1:")
# print(scan1.mediator_partial_width_dm() / scan1.mmed)

# # Example two: multiple variables scanned
# # numpy broadcasting still handles this
# scan2 = DMScalarModelScan(
# mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450]),
# mdm=np.array([1,10,50,100,150,200,250,300,350,400,450]),
# gq=0.25,
# gdm=1.0,
# gl=0.0,
# )
# print("Scan 2:")
# print(scan2.mediator_partial_width_dm() / scan2.mmed)

# Example three: propagators, arrays
scan3 = DMVectorModelScan(mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450]),
mdm=np.array([1,10,50,100,150,200,250,300,350,400,450]),
gq=0.25,
gdm=1.0,
gl=0.0,
)
print("Scan 3:")
print(scan3.mediator_partial_width_quarks())
#print(scan3.propagator_monox_relative())

# Example four: parton-level, no arrays
# scan4 = DMVectorModelScan(mmed=1000,
# mdm=400,
# gq=0.25,
# gdm=1.0,
# gl=0.0,
# )
# print("Scan 4:")
# print(scan4.parton_level_xsec_monox_relative())
# print("done")

# # Example four and a half: propagators, arrays
# scan4p5 = DMVectorModelScan(mmed=3*np.array([1,10,50,100,150,200,250,300,350,400,450]),
# mdm=np.array([1,10,50,100,150,200,250,300,350,400,450]),
# gq=0.25,
# gdm=1.0,
# gl=0.0,
# )
# print("Scan 4.5:")
# print(scan4p5.parton_level_xsec_monox_relative())

# Example four: parton level xsec


# Example five: hadron level xsec


# Example three: Inconsistent array sizes
# broadcasting fails
# -> Write checks to prevent users from doing this
# Kate thinking: should user access this at all? 
# Can we set it up to do grid over these dimensions instead of broadcasting?
# Would be more natural for user's expectation.
# Or else make it very clear what they're getting.
# I guess this is better because it allows uneven grids.
# scan3 = DMScalarModelScan(
# mmed=3*np.array([1,10]),
# mdm=np.array([1,10,50,100,150,200,250,300,350,400,450]),
# gq=0.25,
# gdm=1.0,
# gl=0.0,
# )
# print(scan3.mediator_partial_width_dm() / scan3.mmed)
