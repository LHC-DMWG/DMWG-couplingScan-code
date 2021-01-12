from libc.math cimport sqrt
from libcpp.string cimport string

## Center of mass energy won't change here.
ECM = 13000.**2

# Importing from LHAPDF and wrapping functions and classes
#cdef extern from "LHAPDF/PDF.h" namespace "LHAPDF" :
cdef extern from "LHAPDF/LHAPDF.h" namespace "LHAPDF" :

     # Class
     cdef cppclass PDF :
       PDF()
       double xfxQ2(int id, double x, double q2)

     # Functions
     PDF* mkPDF(const string& setname, int member)
     inline void pathsPrepend(const string& p)

# Making a pdf to use
cdef char* pdfpath = "/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/"
cdef char* lhapdfpath = "/cvmfs/sft.cern.ch/lcg/releases/LCG_97python3/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc9-opt/share/LHAPDF/"
pathsPrepend(pdfpath)
pathsPrepend(lhapdfpath)
# cdef char* pdfset = "NNPDF30_lo_as_0118_nf_4"
cdef char* pdfset = "NNPDF30_nlo_as_0118"
cdef PDF* myPDFSet = mkPDF(pdfset,0)

# Parton-level cross section integrand
cdef double integrand_parton_vector(int n, double[4] args):

     # Four arguments are: S, Gamma, M, mDM
     # where S is variable of integration.
     S = args[0]
     Gamma = args[1]
     M = args[2]
     mDM = args[3]

     # Zero for S < 4 mDM^2, so just return that.
     if S < 4.*mDM**2 : return 0
     numerator = sqrt(S - 4.*mDM**2) * (S + 2.*mDM**2)
     denominator = sqrt(S)*(Gamma**2 * M**2 + (M**2 - S)**2)
     return numerator/denominator

# Hadron-level cross section integrand
cdef double integrand_hadronic_vector(int n, double[6] args):

     # Six arguments, parse here:
     x1 = args[0] 
     x2 = args[1]
     pid = int(args[2])
     Gamma = args[3] 
     M = args[4] 
     mDM = args[5]

     sHat = ECM*x1*x2
     integrand_basic = integrand_parton_vector(4, [sHat, Gamma, M, mDM] )
     total_integrand = myPDFSet.xfxQ2(pid,x1,sHat) * myPDFSet.xfxQ2(-pid,x2,sHat) * integrand_basic
     return 1e8*total_integrand

# Differential cross section to integrate: axial-vector
cdef double integrand_parton_axialvector(int n, double[4] args) :

     # Four arguments are: S, Gamma, M, mDM
     # where S is variable of integration.
     S = args[0]
     Gamma = args[1]
     M = args[2]
     mDM = args[3]

     if S < 4.*mDM**2 : return 0  
     numerator = (S - 4.*mDM**2)**(3./2.)
     denominator = sqrt(S)*(Gamma**2 * M**2 + (M**2 - S)**2)
     return numerator/denominator

cdef double integrand_hadronic_axialvector(int n, double[6] args) :

     # Six arguments, parse here:
     x1 = args[0] 
     x2 = args[1]
     pid = int(args[2])
     Gamma = args[3] 
     M = args[4] 
     mDM = args[5]

     sHat = ECM*x1*x2
     integrand_basic = integrand_parton_axialvector(4, [sHat, Gamma, M, mDM])
     total_integrand = myPDFSet.xfxQ2(pid,x1,sHat) * myPDFSet.xfxQ2(-pid,x2,sHat) * integrand_basic
     return 1e8*total_integrand 
