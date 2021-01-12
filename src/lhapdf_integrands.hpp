#include "LHAPDF/LHAPDF.h"
#include <math.h>
#include <iostream>

class IntegrandHandler {

    public :

        IntegrandHandler(const std::string& setname, double ECM);

        double integrand_parton_vector(int n, double S, double Gamma, double M, double mDM);

        double integrand_hadronic_vector(int n, double x1, double x2, double pid, double Gamma, double M, double mDM);

        double integrand_parton_axialvector(int n, double S, double Gamma, double M, double mDM);

        double integrand_hadronic_axialvector(int n, double x1, double x2, double pid, double Gamma, double M, double mDM);

    private :

        LHAPDF::PDF * m_PDFSet;

        double m_ECM;
};