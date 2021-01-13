#include <pybind11/pybind11.h>
//#include "LHAPDF/LHAPDF.h"
#include <math.h>

#include "lhapdf_integrands.hpp"

IntegrandHandler::IntegrandHandler(const std::string& setname, double ECM) {

  m_PDFSet = LHAPDF::mkPDF(setname,0);
  m_ECM = ECM;

}

// Parton-level cross section integrand: vector
double IntegrandHandler::integrand_parton_vector(int n, double S, double Gamma, double M, double mDM) {

     // Zero for S < 4 mDM^2, so just return that.
     if (S < 4.*pow(mDM,2)) return 0;
     double numerator = sqrt(S - 4.*pow(mDM,2)) * (S + 2.*pow(mDM,2));
     double denominator = sqrt(S)*(pow(Gamma,2) * pow(M,2) + pow((pow(M,2) - S),2));
     return numerator/denominator;
}

// Hadron-level cross section integrand: vector
double IntegrandHandler::integrand_hadronic_vector(int n, double x1, double x2, double pid, double Gamma, double M, double mDM) {

     double sHat = m_ECM*x1*x2;
     double integrand_basic = integrand_parton_vector(4, sHat, Gamma, M, mDM );
     double total_integrand = m_PDFSet->xfxQ2(pid,x1,sHat) * m_PDFSet->xfxQ2(-pid,x2,sHat) * integrand_basic;
     // These numbers are super tiny so scale them up to make this calculable
     // Overall scale doesn't matter, only relative scales
     return 1e8*total_integrand;
 }

// Parton-level cross section integrand: axial-vector
double IntegrandHandler::integrand_parton_axialvector(int n, double S, double Gamma, double M, double mDM) {

     if (S < 4.*pow(mDM,2)) return 0;
     double numerator = pow((S - 4.*pow(mDM,2)),(3./2.));
     double denominator = sqrt(S)*(pow(Gamma,2) * pow(M,2) + pow((pow(M,2) - S),2));
     return numerator/denominator;
 }

// Hadron level cross section integrand: axial-vector
double IntegrandHandler::integrand_hadronic_axialvector(int n, double x1, double x2, double pid, double Gamma, double M, double mDM) {

     double sHat = m_ECM*x1*x2;
     double integrand_basic = integrand_parton_axialvector(4, sHat, Gamma, M, mDM);
     double total_integrand = m_PDFSet->xfxQ2(pid,x1,sHat) * m_PDFSet->xfxQ2(-pid,x2,sHat) * integrand_basic;
     // Again, scale up
     return 1e8*total_integrand;
}


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(lhapdfwrap, m) {
    m.doc() = R"pbdoc(
        Pybind11 for wrapping lhapdf IntegrandHandler
        -----------------------
        .. currentmodule:: lhapdfwrap
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

    py::class_<IntegrandHandler>(m, "IntegrandHandler")
        .def("integrand_parton_vector", &IntegrandHandler::integrand_parton_vector, R"pbdoc(
        Parton-level cross section integrand for vector mediators.)pbdoc")
       .def("integrand_hadronic_vector", &IntegrandHandler::integrand_hadronic_vector, R"pbdoc(
        Hadron-level cross section integrand for vector mediators.)pbdoc")
       .def("integrand_parton_axialvector", &IntegrandHandler::integrand_parton_axialvector, R"pbdoc(
        Parton-level cross section integrand for axial-vector mediators.)pbdoc")
       .def("integrand_hadronic_axialvector", &IntegrandHandler::integrand_hadronic_axialvector, R"pbdoc(
        Hadron-level cross section integrand for axial-vector mediators.)pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
