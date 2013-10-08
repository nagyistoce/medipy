/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef _e7e3635a_542e_4b88_859d_ab2da9eaf14b
#define _e7e3635a_542e_4b88_859d_ab2da9eaf14b

#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

class EM
{
public:
    typedef EM Self;
    
    EM(vnl_vector<double> const & PIi, vnl_matrix<double> const & aij, 
       vnl_matrix<double> const & Mu, std::vector<vnl_matrix<double> > const & Sigma);
    
    ~EM();
    
    void HMCRobusteAtlasFlair(vnl_matrix<double> const & Chain, 
                              vnl_matrix<double> const & ChainAtlas,
                              unsigned int iterations,
                              int FlairImage, int RobusteType, float threshold);
    
    void SegMPMRobusteAtlas(vnl_matrix<double> const & Chain, 
                            vnl_matrix<double> const & AtlasChain,
                            vnl_vector<int> & ChainSeg,
                            vnl_vector<int> & ChainLesions,
                            bool AfficheOutliers);

private:
    vnl_matrix<double> Forward;
    vnl_matrix<double> Backward;
    vnl_matrix<double> ForwardBackward;
    vnl_matrix<double> SommePsi;
    vnl_vector<double> PIi;
    vnl_matrix<double> aij;
    vnl_matrix<double> Mu;
    std::vector<vnl_matrix<double> > Sigma;
    double Changement;
    vnl_matrix<double> MatrixProp;
    std::vector<double> Residus;
    std::vector<double> Probas;
    
    vnl_vector<int> PositionOutliers;
    
    vnl_vector<double> PIi2;
    vnl_matrix<double> aij2;
    vnl_matrix<double> Mu2;
    std::vector<vnl_matrix<double> > Sigma2;
    
    void CalculForwardRobusteAtlas(vnl_matrix<double> const & Chain, 
                                   vnl_matrix<double> const & AtlasChain);
    void CalculBackwardRobusteAtlas(vnl_matrix<double> const & Chain, 
                                    vnl_matrix<double> const & AtlasChain);
    void CalculForwardBackward();
    
    void CalculSommePsiRobusteAtlas(vnl_matrix<double> const & Chain, 
                                    vnl_matrix<double> const & AtlasChain);
    
    void EstimatePIi();
    void EstimateAij();
    
    void EstimateMuRobuste(vnl_matrix<double> const & Chain);
    void EstimateSigmaRobuste(vnl_matrix<double> const & Chain);
    
    void CalculChangement(int iter);
    
    static void MPM(vnl_matrix<double> const & ForwardBackward, 
             vnl_vector<int> & ChainSeg);
    
    static double gaussienne(const double Data, double Mu, double Sigma);
    static double gaussienne(vnl_vector<double> const & Data, 
                      vnl_vector<double> const & Mu, 
                      vnl_matrix<double> const & Sigma);
    static double gaussienne2D(vnl_vector<double> const & Data, 
                        vnl_vector<double> const & Mu, 
                        vnl_matrix<double> const & Sigma);
    static double gaussienne3D(vnl_vector<double> const & Data, 
                        vnl_vector<double> const & Mu, 
                        vnl_matrix<double> const & Sigma);
    
    void CalculResidusAtlasFlair(vnl_matrix<double> const & Chain, 
                                 vnl_matrix<double> const & AtlasChain,
                                 int FlairImage);
    
    void FindPositionOutliers(vnl_matrix<double> const & Chain, int FlairImage,
                              float PourcentageResidus);
    void FindPositionOutliersAtlasSeuil(vnl_matrix<double> const & Chain, 
                                        vnl_matrix<double> const & AtlasChain,
                                        int FlairImage, float Seuil);
};
    
#endif // _e7e3635a_542e_4b88_859d_ab2da9eaf14b
    
