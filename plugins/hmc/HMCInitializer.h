/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _8dfbf825_3e8d_410a_b192_93f36f8329d5
#define _8dfbf825_3e8d_410a_b192_93f36f8329d5

#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "InputParameters.h"

class HMCInitializer
{
public:
    typedef HMCInitializer Self;

    HMCInitializer();
    ~HMCInitializer();
    
    void operator()(vnl_matrix<double> const & Chain, 
                    InputParameters const & params);
    
    vnl_vector<double> const & GetPIi() const;
    vnl_matrix<double> const & Getaij() const;
    vnl_matrix<double> const & GetMu() const;
    std::vector<vnl_matrix<double> > const & GetSigma() const;

private:
    static double gauss1(double Data, double Mu, double Sigma);
    
    static void Histo_EM(vnl_vector<double> & histo, int quantification, int NbClasses,
                         vnl_vector<double> & Pi, vnl_vector<double> & Mu,
                         vnl_vector<double> & Sigma);
    
    static void Histogramme(vnl_vector<double> const & Chaine, int Taille, 
                            int NbClasses, vnl_vector<double> & Mu, 
                            vnl_vector<double> & Sigma);
    
    static void KMeans(vnl_matrix<double> const & Chain, int TailleEch, 
                       vnl_matrix<double> & Mu, vnl_matrix<double> & Sig,
                       InputParameters const & params);
    
    vnl_vector<double> PIi;
    vnl_matrix<double> aij;
    vnl_matrix<double> Mu;
    std::vector<vnl_matrix<double> > Sigma;
    
    void IniPI(InputParameters const & params);
    void Iniaij(InputParameters const & params);
    void IniKMeans(vnl_matrix<double> const & Chain, 
                   InputParameters const & params);
};

#endif // _8dfbf825_3e8d_410a_b192_93f36f8329d5
