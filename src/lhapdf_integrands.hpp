#include "LHAPDF/LHAPDF.h"
#include <math.h>
#include <iostream>

class IntegrandHandler {

    public :

        IntegrandHandler(const std::string& setname, double ECM);

        double integrand_parton_vector(double S, double Gamma, double M, double mDM);

        double integrand_hadronic_vector(double x1, double x2, double pid, double Gamma, double M, double mDM);

        double integrand_parton_axialvector(double S, double Gamma, double M, double mDM);

        double integrand_hadronic_axialvector(double x1, double x2, double pid, double Gamma, double M, double mDM);

    private :

        LHAPDF::PDF * m_PDFSet;

        double m_ECM;
};