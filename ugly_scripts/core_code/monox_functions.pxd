# This is sort of a header file.
# Tell it to use c++:
# distutils: language=c++
# cython: language_level=3, boundscheck=False

# Four arguments are: S, Gamma, M, mDM
# where S is variable of integration
cdef double integrand_parton_vector(int, double[4])

# Six arguments are: x1, x2, pid, Gamma, M, mDM
# where x1 and x2 will both be integrated over.
cdef double integrand_hadronic_vector(int, double[6])

# Four arguments are: S, Gamma, M, mDM
# where S is variable of integration
cdef double integrand_parton_axialvector(int, double[4])

# Six arguments are: x1, x2, pid, Gamma, M, mDM
# where x1 and x2 will both be integrated over.
cdef double integrand_hadronic_axialvector(int, double[6])