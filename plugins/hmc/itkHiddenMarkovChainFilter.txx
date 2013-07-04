#ifndef _itkHiddenMarkovChainFilter_txx
#define _itkHiddenMarkovChainFilter_txx

#define math_Pi 3.141592654

#include "itkImageToImageFilter.h"


#include <vector>
#include <math.h>
#include <limits>


#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkPathIterator.h>
#include <itkPathConstIterator.h>




namespace itk{


template<typename TInputImage, typename TOutputImage>
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::HiddenMarkovChainFilter()
:m_Number_images(0), m_Atlas_bool(true), m_Number_iter(5), m_Number_classes(3), m_Criterion_outliers(1), m_Criterion_outliers_value(0.05), m_Flair_bool(true), m_Position_Flair(0)
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(0);
 
    this->SetNthOutput( 0, this->MakeOutput(0) );
    this->SetNthOutput( 1, this->MakeOutput(1) );
}

template<typename TInputImage, typename TOutputImage>
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::~HiddenMarkovChainFilter()
{
}

template<typename TInputImage, typename TOutputImage>
void HiddenMarkovChainFilter<TInputImage, TOutputImage>
::SetInputImage(unsigned int idx, TInputImage* image)
{
    /*typename InputImageType::IndexType pixelIndex;
    pixelIndex[0] = 128;
    pixelIndex[1] = 128;
    pixelIndex[2] = 128;*/
    SetNthInput(idx, image);
    //std::cout<<idx<<", "<<image->GetPixel(pixelIndex)<<std::endl;
}

template<typename TInputImage, typename TOutputImage>
typename TInputImage::Pointer HiddenMarkovChainFilter<TInputImage, TOutputImage>
::GetInputImage(unsigned int idx)
{
    return static_cast<TInputImage*>(this->ProcessObject::GetInput(idx));
}

template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Number of input images : " << this->m_Number_images << "\n";
    os << indent << "Boolean for atlas : " << this->m_Atlas_bool << "\n";
    os << indent << "Number of iterations : " << this->m_Number_iter << "\n";
    os << indent << "Number of classes : " << this->m_Number_classes << "\n";
    os << indent << "Type of criterion for outliers (0 for percentage and 1 for threshold) : " << this->m_Criterion_outliers << "\n";
    os << indent << "Value for the criterion for outliers : " << this->m_Criterion_outliers_value << "\n";
}


template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::GenerateData()
{
    /*typename InputImageType::Pointer input0 = this->GetInputImage(0);
    ImageRegionConstIterator<InputImageType> input0_iterator(input0, input0->GetLargestPossibleRegion());

    typename InputImageType::Pointer input1 = this->GetInputImage(1);
    ImageRegionConstIterator<InputImageType> input1_iterator(input1, input1->GetLargestPossibleRegion());

    typename OutputImageType::Pointer output_seg_image = this->GetOuputSegImage();
    output_seg_image->SetRegions(input0->GetLargestPossibleRegion());
    output_seg_image->Allocate();
    ImageRegionIterator<OutputImageType> output_seg_iterator(output_seg_image, output_seg_image->GetLargestPossibleRegion());

    typename OutputImageType::Pointer output_outliers_image = this->GetOuputOutliersImage();
    output_outliers_image->SetRegions(input1->GetLargestPossibleRegion());
    output_outliers_image->Allocate();
    ImageRegionIterator<OutputImageType> output_outliers_iterator(output_outliers_image, output_outliers_image->GetLargestPossibleRegion());

    typename PathType::Pointer hilbert_path = PathType::New();
    hilbert_path->SetHilbertOrder(m_Taille_cube);
    hilbert_path->Initialize();
    PathConstIterator<OutputImageType, PathType> hilbert_iterator(input0, hilbert_path);

    while(!output_seg_iterator.IsAtEnd()){
        output_seg_iterator.Set(input0_iterator.Get());
        ++output_seg_iterator;
        ++input0_iterator;
    }

    while(!output_outliers_iterator.IsAtEnd()){
        output_outliers_iterator.Set(input1_iterator.Get());
        ++output_outliers_iterator;
        ++input1_iterator;
    }

    unsigned int i=1;
    while(i<10){
        ++hilbert_iterator;
        if(hilbert_iterator.Get()>0){
            std::cout<<hilbert_iterator.Get()<<std::endl;
            i++;
        }
    }*/

    std::cout<<"Start !!!"<<std::endl;

    //Input iterators construction and input images reading 
    typename PathType::Pointer hilbert_path = PathType::New(); //Creation of hilbert-peano iterator
    hilbert_path->SetHilbertOrder(m_Taille_cube);
    hilbert_path->Initialize();
    
    std::vector<ConstIteratorType *> inputs_its;
    std::vector<ConstIteratorType *> atlas_its;

    ConstIteratorType * mask_it;

    for(unsigned int i=0; i < m_Number_images; i++){
        typename InputImageType::Pointer input = this->GetInputImage(i);
        ConstIteratorType *it = new ConstIteratorType(input, hilbert_path);
        it->GoToBegin();
        inputs_its.push_back(it);
    }

    if(m_Atlas_bool){
        for(unsigned int i=m_Number_images; i < m_Number_images+m_Number_classes; i++){
            typename InputImageType::Pointer input = this->GetInputImage(i);
            ConstIteratorType *it = new ConstIteratorType(input, hilbert_path);
            it->GoToBegin();
            atlas_its.push_back(it);
        }
        mask_it = new ConstIteratorType(this->GetInputImage(m_Number_images+m_Number_classes), hilbert_path);
        mask_it->GoToBegin();
    }else{
        mask_it = new ConstIteratorType(this->GetInputImage(m_Number_images), hilbert_path);
        mask_it->GoToBegin();
    }


    //Chains creation
    std::cout<<"Chain creation"<<std::endl;
    std::vector< std::vector<double> > Chain;
    std::vector< std::vector<double> > Chain_Atlas;
 
    Chains_Creation(Chain, Chain_Atlas, mask_it, inputs_its, atlas_its);
    //Chain_Atlas.clear();
    //Chain_Atlas.resize(Chain.size(), std::vector<double>(m_Number_images, 1.));


    //Chain initialisation
    std::cout<<"Chain initialisation"<<std::endl;    
    std::vector< std::vector<double> > Moyenne(m_Number_images, std::vector<double>(m_Number_classes, 0.));
    std::vector< std::vector< std::vector<double> > > Variance(m_Number_classes, std::vector< std::vector<double> >(m_Number_images, std::vector<double>(m_Number_images, 0.)));
    std::vector<double> Proba_ini(m_Number_classes, 0.);
    std::vector< std::vector<double> > aij(m_Number_classes, std::vector<double>(m_Number_classes, 0.));

    Chain_Itialisation(Chain, Moyenne, Variance, Proba_ini, aij);
    
    std::vector<unsigned int> Chain_seg(Chain.size(), 0);
    std::vector<unsigned int> Chain_outliers(Chain.size(), 0);

    //Expectation-Maximization
    std::cout<<"Expectation-Maximization"<<std::endl;        
    Chain_EM(Chain, Chain_Atlas, Moyenne, Variance, Proba_ini, aij, Chain_seg, Chain_outliers);

    //MPM
    std::cout<<"MPM"<<std::endl;
    Chain_MPM(Chain, Chain_Atlas, Moyenne, Variance, Proba_ini, aij, Chain_seg, Chain_outliers);

/*    std::vector<double> Chain_seg(Chain.size(), 0);
    std::vector<double> Chain_outliers(Chain.size(), 0);
    for(unsigned int i=0; i<Chain.size(); i++){
        Chain_seg[i]=Chain[i][0];
        Chain_outliers[i]=Chain[i][1];
    }*/
    

    //Output iterators construction nad output images generation
    std::cout<<"Images reconstruction"<<std::endl;
    typename OutputImageType::Pointer output_seg_image = this->GetOuputSegImage();
    output_seg_image->SetRegions(this->GetInputImage(0)->GetLargestPossibleRegion());
    output_seg_image->Allocate();
    IteratorType output_seg_it(output_seg_image, hilbert_path);
    output_seg_it.GoToBegin();


    typename OutputImageType::Pointer output_outliers_image = this->GetOuputOutliersImage();
    output_outliers_image->SetRegions(this->GetInputImage(0)->GetLargestPossibleRegion());
    output_outliers_image->Allocate();
    IteratorType output_outliers_it(output_outliers_image, hilbert_path);
    output_outliers_it.GoToBegin();

    Images_Reconstruction(Chain_seg, Chain_outliers, mask_it, output_seg_it, output_outliers_it);

    /*mask_it->GoToBegin();
    while(!output_seg_it.IsAtEnd()){
        output_seg_it.Set(atlas_its[0]->Get());
        ++output_seg_it;
        output_outliers_it.Set(atlas_its[1]->Get());
        ++output_outliers_it;
        ++(*atlas_its[0]);
        ++(*atlas_its[1]);
    }*/
    /*while(!output_outliers_it.IsAtEnd()){
        output_outliers_it.Set(1);
        ++output_outliers_it;
    }*/
/*
    mask_it->GoToBegin();
    for(unsigned int j=0; j!=m_Number_images; j++)
        inputs_its[j]->GoToBegin();
    while(!(mask_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            double toto=0.;
            double big_toto=1e25;
            for(k=0; k!=m_Number_classes; k++){
                for(unsigned int j=0; j!=m_Number_images; j++){
                    toto+=(inputs_its[j]->Get()-MuProv[j][k]);
                }
                if(toto<big_toto){
                    output_seg_it.Set(k);
                    big_toto=toto;
                }
            }
        }
        output_outliers_it.Set(mask_it->Get());
        ++output_outliers_it;
        ++(*mask_it);
        for(unsigned int j=0; j!=m_Number_images; j++){
            ++(*inputs_its[j]);
        }
        ++output_seg_it;
    }*/

    std::cout<<"End !!!"<<std::endl;

}


template<typename TInputImage, typename TOutputImage>
DataObject::Pointer
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::MakeOutput(unsigned int idx)
{
    DataObject::Pointer output;
 
    switch ( idx )
        {
        case 0:
            output = ( TOutputImage::New() ).GetPointer();
            break;
        case 1:
            output = ( TOutputImage::New() ).GetPointer();
            break;
        default:
            std::cerr << "No output " << idx << std::endl;
            output = NULL;
        break;
        }
    return output.GetPointer();
}

 
template<typename TInputImage, typename TOutputImage>
TOutputImage*
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::GetOuputSegImage()
{
    return dynamic_cast< TOutputImage * >(this->ProcessObject::GetOutput(0) );
}

 
template<typename TInputImage, typename TOutputImage>
TOutputImage*
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::GetOuputOutliersImage()
{
    return dynamic_cast< TOutputImage * >(this->ProcessObject::GetOutput(1) );
}


template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::Chains_Creation(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Chain_Atlas, ConstIteratorType *mask_it, std::vector<ConstIteratorType *> &inputs_its, std::vector<ConstIteratorType *> &atlas_its)
{
    //std::vector<double> normalisation(2);
    //normalisation[0] = pow(2, 3);
    //normalisation[1] = pow(2, 7);
    mask_it->GoToBegin();
    for(unsigned int i=0; i<m_Number_images; i++)
        inputs_its[i]->GoToBegin();
    while(!mask_it->IsAtEnd()){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            std::vector<double> temp(m_Number_images, 0.);
            for(unsigned int i=0; i<m_Number_images; i++)
                temp[i] = (inputs_its[i]->Get());// / normalisation[i];
            Chain.push_back(temp);
            if(m_Atlas_bool){
                std::vector<double> temp_atlas(m_Number_classes, 0.);
                for(unsigned int i=0; i<m_Number_classes; i++)
                    temp_atlas[i] = atlas_its[i]->Get();
                Chain_Atlas.push_back(temp_atlas);
            }  
        }
        ++(*mask_it);
        for(unsigned int i=0; i<m_Number_images; i++)
            ++(*inputs_its[i]);
        if(m_Atlas_bool){
            for(unsigned int i=0; i<m_Number_classes; i++)
                ++(*atlas_its[i]);
        }
    }
    if(!(m_Atlas_bool)){
        Chain_Atlas.clear();
        Chain_Atlas.resize(Chain.size(), std::vector<double>(m_Number_classes, 1.));
    }
}


template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::Chain_Itialisation(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Moyenne, std::vector< std::vector< std::vector<double> > > &Variance, std::vector<double> &Proba_ini,  std::vector< std::vector<double> > &aij)
{
    unsigned int classe1, classe2, length_Chain, length_histo;


    //Initialisation of Proba_ini//
    for(classe1=0; classe1<m_Number_classes; classe1++){
        Proba_ini[classe1] = 1. / m_Number_classes;
    }

    //Initialisation of aij//
    for(classe1=0; classe1<m_Number_classes; classe1++){
        for(classe2=0; classe2<m_Number_classes; classe2++){
            if(classe1==classe2){
                aij[classe1][classe2] = 3. / 4;
            }else{
                aij[classe1][classe2] = 1. / (4 * (m_Number_classes - 1));
            }
        }
    }

    //Initialisation of noise parameters//
    std::vector<double> Chain_mini(Chain.size(), 0.);
    for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
        Chain_mini[length_Chain] = Chain[length_Chain][0];
    }

    std::vector<double> Moyenne_c(m_Number_classes, 0.);
    std::vector< std::vector<double> > Standard_deviation_c(m_Number_classes, std::vector<double>(m_Number_classes, 0.));
    std::vector<double> Pi(m_Number_classes, 0.);

    //Min and Max of Chain_mini
    double min_c, max_c;
    min_c = Chain_mini[0];
    max_c = min_c;

    for(length_Chain=0; length_Chain<Chain_mini.size(); length_Chain++){
        if(Chain_mini[length_Chain]>max_c){
            max_c = Chain_mini[length_Chain];
        }
        if(Chain_mini[length_Chain]<min_c){
            min_c = Chain_mini[length_Chain];
        }
    }

    //Histogram of Chaine_mini
    unsigned int quantification = 128;
    std::vector<double> histo(quantification, 0.);
    
    for(length_Chain=0; length_Chain<Chain_mini.size(); length_Chain++){
        histo[(int)rint((Chain_mini[length_Chain] - min_c) * ((float)quantification - 1) / (max_c - min_c))]++;
    }

    for(length_histo=0; length_histo<histo.size(); length_histo++){
        histo[length_histo] /= Chain_mini.size();
    }

    //Evaluation of Mean_c (mean of Chain_mini)
    double gap = 1. / (2 * m_Number_classes);
    double temp;

    for(classe1=0; classe1<m_Number_classes; classe1++){
        temp = 0.;
        length_histo = 0.;
        while(temp < (gap * ( 2 * classe1 + 1))){
            temp += histo[length_histo];
            length_histo++;
        }
        Moyenne_c[classe1] = (double)length_histo;       
    }

    //Evaluation of Standard_deviation_c (standard deviation of Chain_mini)
    unsigned int res, res2;
    double Tinter, ecart;
    res = 0.;
    for(classe1=0; classe1<m_Number_classes; classe1++){
        Tinter = 0.;
        ecart = 0.;

        if(classe1<m_Number_classes-1){
            res2 = (int)(rint((Moyenne_c[classe1+1] - Moyenne_c[classe1]) / 2) + Moyenne_c[classe1]);
        }else{
            res2 = histo.size();
        }
        length_histo = res;

        while(length_histo<res2){
            ecart += (length_histo - Moyenne_c[classe1]) * (length_histo - Moyenne_c[classe1]) * histo[length_histo];
            Tinter += histo[length_histo];
            length_histo++;
        }

        if((Tinter!=0)&&(ecart!=0)){
            Pi[classe1] = Tinter;
            Standard_deviation_c[classe1][classe1] = sqrt(ecart);
        }else{
            Pi[classe1] = 1. / (m_Number_classes + 1);
            Standard_deviation_c[classe1][classe1] = 1. / sqrt(2 * math_Pi);
        }
        res=res2;
    }

    //EM on histogram
    unsigned int iter=0;
    double err=1000.;
    unsigned int quanti_lvl;
    std::vector< std::vector<double> > apriori(quantification, std::vector<double>(m_Number_classes, 0.));
    std::vector<double> apriori_sum(quantification, 0.);
    std::vector<double> data(1, 0.);
    std::vector<double> data_mean(1, 0.);
    std::vector< std::vector<double> > data_std(1, std::vector<double>(1, 0.));
    std::vector<double> Moyenne_c_temp(m_Number_classes, 0.);
    std::vector< std::vector<double> > Standard_deviation_c_temp(m_Number_classes, std::vector<double>(m_Number_classes, 0.));
    std::vector<double> Pi_temp(m_Number_classes, 0.);

    while((err>1e-5)&&(iter<1000)){
        iter++;
        temp=0.;

        for(quanti_lvl=0; quanti_lvl<quantification; quanti_lvl++){
            apriori_sum[quanti_lvl]=0.;
            for(classe1=0; classe1<m_Number_classes; classe1++){
                data[0] = (double) quanti_lvl;
                data_mean[0] = Moyenne_c[classe1];
                data_std[0][0] = Standard_deviation_c[classe1][classe1];
                apriori[quanti_lvl][classe1] = Pi[classe1] * Gaussian(data, data_mean, data_std);
                apriori_sum[quanti_lvl] += apriori[quanti_lvl][classe1];
            }
            temp += histo[quanti_lvl];
            for(classe1=0; classe1<m_Number_classes; classe1++){
                if((apriori_sum[quanti_lvl]>1.e-45) && (apriori[quanti_lvl][classe1]>1.e-45)){
                    apriori[quanti_lvl][classe1] /= apriori_sum[quanti_lvl];
                }else{
                    apriori[quanti_lvl][classe1] = 0.;
                }
            }
        }

        //New probability a priori (Pi_temp)
        for(classe1=0; classe1<m_Number_classes; classe1++){
            Pi_temp[classe1] = 0.;
            for(quanti_lvl=0; quanti_lvl<quantification; quanti_lvl++){
                Pi_temp[classe1] += histo[quanti_lvl] * apriori[quanti_lvl][classe1];
            }
            if(Pi_temp[classe1]!=0){
                Pi_temp[classe1] /= temp;
            }else{
                Pi_temp[classe1] = 1. / (m_Number_classes + 1);
            }
        }

        //New mean (Moyenne_c_temp)
        for(classe1=0; classe1<m_Number_classes; classe1++){
            Moyenne_c_temp[classe1] = 0.;
            temp = 0.;
            for(quanti_lvl=0; quanti_lvl<quantification; quanti_lvl++){
                Moyenne_c_temp[classe1] += histo[quanti_lvl] * quanti_lvl * apriori[quanti_lvl][classe1];
                temp += histo[quanti_lvl] * apriori[quanti_lvl][classe1];
            }
            Moyenne_c_temp[classe1] /= temp;
        }

        //New standard deviation (Standard_deviation_c_temp)
        for(classe1=0; classe1<m_Number_classes; classe1++){
            Standard_deviation_c_temp[classe1][classe1] = 0.;
            temp = 0.;
            for(quanti_lvl=0; quanti_lvl<quantification; quanti_lvl++){
                Standard_deviation_c_temp[classe1][classe1] += histo[quanti_lvl] * (quanti_lvl - Moyenne_c_temp[classe1]) * (quanti_lvl - Moyenne_c_temp[classe1]) * apriori[quanti_lvl][classe1];
                temp += histo[quanti_lvl] * apriori[quanti_lvl][classe1];
            }
            if((Standard_deviation_c_temp[classe1][classe1]/temp)>(1. / (2 * math_Pi))){
                Standard_deviation_c_temp[classe1][classe1] = sqrt(Standard_deviation_c_temp[classe1][classe1]/temp);
            }else{
                Standard_deviation_c_temp[classe1][classe1] = 1. / sqrt(2 * math_Pi);
            }
        } 

        //Error evaluation
        temp = 0.;
        err = 0.;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            err += (Pi[classe1] - Pi_temp[classe1]) * (Pi[classe1] - Pi_temp[classe1]) + (Moyenne_c[classe1] - Moyenne_c_temp[classe1]) * (Moyenne_c[classe1] - Moyenne_c_temp[classe1]) + (Standard_deviation_c[classe1][classe1] - Standard_deviation_c_temp[classe1][classe1]) * (Standard_deviation_c[classe1][classe1] - Standard_deviation_c_temp[classe1][classe1]);
            temp += Pi[classe1] + Moyenne_c[classe1] + Standard_deviation_c[classe1][classe1];
        }
        Pi = Pi_temp;
        Moyenne_c = Moyenne_c_temp;
        Standard_deviation_c = Standard_deviation_c_temp;
    }

    temp = (double)(max_c - min_c) / quantification;
    for(classe1=0; classe1<m_Number_classes; classe1++){
        Moyenne_c[classe1] = min_c + Moyenne_c[classe1] * temp;
        Standard_deviation_c[classe1][classe1] = Standard_deviation_c[classe1][classe1] * temp;
    }


    apriori.clear();
    apriori_sum.clear();

    //Segmentation chain from EM with maximum likelihood
    std::vector<unsigned int> Seg_Chain(Chain.size(), std::numeric_limits<unsigned int>::max());
    for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
        temp = -1.;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            data[0] = Chain_mini[length_Chain];
            data_mean[0] = Moyenne_c[classe1];
            data_std[0][0] = Standard_deviation_c[classe1][classe1];
            if(Gaussian(data, data_mean, data_std)>temp){
                temp = Gaussian(data, data_mean, data_std);
                Seg_Chain[length_Chain] = classe1;
            }
        }
    }

    

    //Evaluation of number of voxels for each classe thanks to segmentaion chain
    std::vector<unsigned int> NbVox(m_Number_classes, 0);
    for(classe1=0; classe1<m_Number_classes; classe1++){
        temp = 0;
        for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
            if(Seg_Chain[length_Chain]==classe1){
                temp++;
            }
        }
        NbVox[classe1] = temp;
    }

    //Evaluation of the mean of each classe
    unsigned int mod1;   
    for(classe1=0; classe1<m_Number_classes; classe1++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Moyenne[mod1][classe1] = 0.;
            for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
                if(Seg_Chain[length_Chain]==classe1){
                    Moyenne[mod1][classe1] += Chain[length_Chain][mod1];
                }
            }
            Moyenne[mod1][classe1] /= NbVox[classe1];
        }
    }

    //Evaluation of the variance of each classe
    for(classe1=0; classe1<m_Number_classes; classe1++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Variance[classe1][mod1][mod1] = 0.;
            for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
                if(Seg_Chain[length_Chain]==classe1){
                    Variance[classe1][mod1][mod1] += (Chain[length_Chain][mod1] - Moyenne[mod1][classe1]) * (Chain[length_Chain][mod1] - Moyenne[mod1][classe1]);
                }
            }
            Variance[classe1][mod1][mod1] = sqrt(Variance[classe1][mod1][mod1] / (NbVox[classe1] - 1.));
        }
    }

/*    std::cout<<"Moyenne = "<<std::endl;
    for(mod1=0; mod1<m_Number_images; mod1++){
        std::cout<<"\t"; 
        for(classe1=0; classe1<m_Number_classes; classe1++){
            std::cout<<Moyenne[mod1][classe1]<<"\t";
        }
        std::cout<<std::endl;
    }

    std::cout<<"Variance = "<<std::endl; 
    for(classe1=0; classe1<m_Number_classes; classe1++){
        std::cout<<"\tClasse "<<classe1+1<<std::endl;
        for(mod1=0; mod1<m_Number_images; mod1++){
            std::cout<<"\t\t"; 
            for(mod2=0; mod2<m_Number_images; mod2++){
                std::cout<<pow(Variance[classe1][mod1][mod2], 2)<<"\t";
            }
            std::cout<<std::endl;
        }
    }*/


    //K-means//
    std::cout<<"K-means"<<std::endl;

    //Centroid creation
    std::vector< std::vector<double> > Centroids(m_Number_classes, std::vector<double>(2 * m_Number_images, 0.));
    for(classe1=0; classe1<m_Number_classes; classe1++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Centroids[classe1][mod1] = Moyenne[mod1][classe1];
            Centroids[classe1][m_Number_images + mod1] = Variance[classe1][mod1][mod1];
        }
    }    

    //Data creation (mean and variance for a segment)
    unsigned int length_segment = 8;
    std::vector< std::vector<double> > Chain_segment(((Chain.size() / length_segment) + 1), std::vector<double>(2 * m_Number_images, 0.));

    for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Chain_segment[length_Chain/length_segment][mod1] += (double)(Chain[length_Chain][mod1] / (double)(length_segment));
        }
    }

    for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Chain_segment[length_Chain/length_segment][m_Number_images+mod1] += (double)((Chain[length_Chain][mod1] - Chain_segment[length_Chain/length_segment][mod1]) * (Chain[length_Chain][mod1] - Chain_segment[length_Chain/length_segment][mod1]) / (double)(length_segment));
        }
    }

    for(length_Chain=0; length_Chain<Chain_segment.size(); length_Chain++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Chain_segment[length_Chain][m_Number_images+mod1] = sqrt(Chain_segment[length_Chain][m_Number_images+mod1]);
        }
    }


    std::vector<double> Norme(2, 0.);
    Norme[0] = 1000.;
    iter = 0;
    unsigned int par;


    std::vector< std::vector<double> > Centroids_temp(m_Number_classes, std::vector<double>(2 * m_Number_images, 0.));


    while((fabs(Norme[0] - Norme[1])>1e-30) && (iter<1000)){

        par = iter % 2;

        NbVox.clear();
        NbVox.resize(m_Number_classes, 0);
   

        for(length_Chain=0; length_Chain<Chain_segment.size(); length_Chain++){
            min_c = std::numeric_limits<double>::max();
            classe2 = 0.;
            for(classe1=0; classe1<m_Number_classes; classe1++){
                temp = 0.;
                for(mod1=0; mod1<m_Number_images; mod1++){
                    temp += (Centroids[classe1][mod1] - Chain_segment[length_Chain][mod1]) * (Centroids[classe1][mod1] - Chain_segment[length_Chain][mod1]);
                }
                if(temp<min_c){
                    min_c = temp;
                    classe2 = classe1;
                }
            }
            NbVox[classe2]++;
            for(mod1=0; mod1<(2 * m_Number_images); mod1++){
                Centroids_temp[classe2][mod1] += Chain_segment[length_Chain][mod1];
            }
        }


        Norme[par] = 0.;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            for(mod1=0; mod1<(2 * m_Number_images); mod1++){
                Centroids_temp[classe1][mod1] = (1. / (double)NbVox[classe1]) * Centroids_temp[classe1][mod1];
                Norme[par] += (Centroids_temp[classe1][mod1] - Centroids[classe1][mod1]) * (Centroids_temp[classe1][mod1] - Centroids[classe1][mod1]);
                Centroids[classe1][mod1] = Centroids_temp[classe1][mod1];
                Centroids_temp[classe1][mod1] = 0.;
            }
        }
        Norme[par] = sqrt(Norme[par]);        
        iter++;

    }

    std::cout<<"iter = "<<iter<<std::endl;

    for(classe1=0; classe1<m_Number_classes; classe1++){
        for(mod1=0; mod1<m_Number_images; mod1++){
            Moyenne[mod1][classe1] = Centroids[classe1][mod1];
            Variance[classe1][mod1][mod1] = pow(Centroids[classe1][m_Number_images + mod1], 2);
        }
    } 

       
    /*std::cout<<"Pi = "<<std::endl; 
    for(classe1=0; classe1<m_Number_classes; classe1++){
        std::cout<<Pi[classe1]<<std::endl;
    }*/

    /*std::cout<<"Moyenne = "<<std::endl;
    for(mod1=0; mod1<m_Number_images; mod1++){
        std::cout<<"\t"; 
        for(classe1=0; classe1<m_Number_classes; classe1++){
            std::cout<<Moyenne[mod1][classe1]<<"\t";
        }
        std::cout<<std::endl;
    }

    std::cout<<"Variance = "<<std::endl; 
    for(classe1=0; classe1<m_Number_classes; classe1++){
        std::cout<<"\tClasse "<<classe1+1<<std::endl;
        for(mod1=0; mod1<m_Number_images; mod1++){
            std::cout<<"\t\t"; 
            for(unsigned mod2=0; mod2<m_Number_images; mod2++){
                std::cout<<Variance[classe1][mod1][mod2]<<"\t";
            }
            std::cout<<std::endl;
        }
    }*/

    
}


template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::Chain_EM(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Chain_Atlas, std::vector< std::vector<double> > &Moyenne, std::vector< std::vector< std::vector<double> > > &Variance, std::vector<double> &Proba_ini,  std::vector< std::vector<double> > &aij, std::vector<unsigned int> &Chain_seg , std::vector<unsigned int> &Chain_outliers)
{
    unsigned int iter, length_Chain, classe1, classe2, mod1, mod2;
    double change, normalisation, temp;

    std::vector<double> data(m_Number_images, 0.);
    std::vector<double> moyenne_c(m_Number_images, 0.);
    std::vector< std::vector<double> > variance_c(m_Number_images, std::vector<double>(m_Number_images, 0.));

    std::vector<double> data_s(1, 0.);
    std::vector<double> moyenne_c_s(1, 0.);
    std::vector< std::vector<double> > variance_c_s(1, std::vector<double>(1, 0.));

    std::vector< std::vector<double> > Forward(m_Number_classes, std::vector<double> (Chain.size(), 0.));
    std::vector< std::vector<double> > Backward(m_Number_classes, std::vector<double> (Chain.size(), 0.));
    std::vector< std::vector<double> > ForwardBackward(m_Number_classes, std::vector<double> (Chain.size(), 0.));

    std::vector< std::vector<double> > Moyenne2(m_Number_images, std::vector<double>(m_Number_classes, 0.));
    std::vector< std::vector< std::vector<double> > > Variance2(m_Number_classes, std::vector< std::vector<double> >(m_Number_images, std::vector<double>(m_Number_images, 0.)));
    std::vector<double> Proba_ini2(m_Number_classes, 0.);
    std::vector< std::vector<double> > aij2(m_Number_classes, std::vector<double>(m_Number_classes, 0.));

    iter = 0;
    change = 1000.;
    while((iter<m_Number_iter) && (change>0.01)){

        std::cout<<"iteration = "<<iter+1<<std::endl;

        //Forward probabilities estimation//
        normalisation = 0.;

        //length_Chain=0
        for(classe1=0; classe1<m_Number_classes; classe1++){

            for(mod1=0; mod1<m_Number_images; mod1++){
                data[mod1] = Chain[0][mod1];
                moyenne_c[mod1] = Moyenne[mod1][classe1];
                for(mod2=0; mod2<m_Number_images; mod2++){
                    variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
                }
            }

            if(Chain_outliers[0]==0){
                Forward[classe1][0] = Proba_ini[classe1] * Chain_Atlas[0][classe1] * Gaussian(data, moyenne_c, variance_c);
                normalisation += Forward[classe1][0];
            }else{
                Forward[classe1][0] = Proba_ini[classe1] * Chain_Atlas[0][classe1];
                normalisation += Forward[classe1][0];
            }
        }

        //Normalisation
        for(classe1=0; classe1<m_Number_classes; classe1++){
            if(normalisation!=0.){
                Forward[classe1][0] /= normalisation;
            }else{
                std::cout<<"Division by zero in Forward estimation"<<std::endl;
                Forward[classe1][0] = 1. / m_Number_classes;
            }
        }

        //length_Chain!=0
        for(length_Chain=1; length_Chain<Chain.size(); length_Chain++){

            normalisation = 0.;
            
            for(classe1=0; classe1<m_Number_classes; classe1++){
                temp = 0.;

                for(mod1=0; mod1<m_Number_images; mod1++){
                    data[mod1] = Chain[length_Chain][mod1];
                    moyenne_c[mod1] = Moyenne[mod1][classe1];
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
                    }
                }

                for(classe2=0; classe2<m_Number_classes; classe2++){
                    temp += Forward[classe2][length_Chain-1] * aij[classe2][classe1];
                }

                if(Chain_outliers[length_Chain]==0){
                    Forward[classe1][length_Chain] = temp * Chain_Atlas[length_Chain][classe1] * Gaussian(data, moyenne_c, variance_c);
                }else{
                    Forward[classe1][length_Chain] = temp * Chain_Atlas[length_Chain][classe1];
                }

                normalisation += Forward[classe1][length_Chain];
            }

            //Normalisation
            for(classe1=0; classe1<m_Number_classes; classe1++){
                if(normalisation!=0.){
                    Forward[classe1][length_Chain] /= normalisation;
                }else{
                    std::cout<<"Division by zero in Forward estimation"<<std::endl;
                    Forward[classe1][length_Chain] = 1. / m_Number_classes;
                }
            }
        }


        //Backward probabilities estimation//
        normalisation = 0.;

        //length_Chain = end
        for(classe1=0; classe1<m_Number_classes; classe1++){
            Backward[classe1][Chain.size()-1] = 1.;
        }

        //length_Chain != end
        for(int length_Chain=Chain.size()-2; length_Chain>=0; length_Chain--){

            normalisation = 0.;

            for(mod1=0; mod1<m_Number_images; mod1++){
                data[mod1] = Chain[length_Chain+1][mod1];
            }
            
            for(classe1=0; classe1<m_Number_classes; classe1++){
                temp = 0.;

                for(classe2=0; classe2<m_Number_classes; classe2++){

                    for(mod1=0; mod1<m_Number_images; mod1++){
                        moyenne_c[mod1] = Moyenne[mod1][classe2];
                        for(mod2=0; mod2<m_Number_images; mod2++){
                            variance_c[mod1][mod2] = Variance[classe2][mod1][mod2];
                        }
                    }

                    if(Chain_outliers[length_Chain+1]==0){
                        temp += Backward[classe2][length_Chain+1] * Chain_Atlas[length_Chain+1][classe2] * Gaussian(data, moyenne_c, variance_c) * aij[classe1][classe2];
                    }else{
                        temp += Backward[classe2][length_Chain+1] * Chain_Atlas[length_Chain+1][classe2] * aij[classe1][classe2];
                    }

                }
                
                Backward[classe1][length_Chain] = temp;
                normalisation += temp;               
                
            }

            //Normalisation
            for(classe1=0; classe1<m_Number_classes; classe1++){
                if(normalisation!=0.){
                    Backward[classe1][length_Chain] /= normalisation;
                }else{
                    std::cout<<"Division by zero in Forward estimation"<<std::endl;
                    Backward[classe1][length_Chain] = 1. / sqrt(m_Number_classes);
                }
            }
        }

        //Forward-Backward estimation, proba marginales a posterioiri//
        unsigned int length_Chain;

        for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
            normalisation = 0.;
            for(classe1=0; classe1<m_Number_classes; classe1++){
                ForwardBackward[classe1][length_Chain] = Forward[classe1][length_Chain] * Backward[classe1][length_Chain];
                normalisation += ForwardBackward[classe1][length_Chain];
            }
            for(classe1=0; classe1<m_Number_classes; classe1++){
                ForwardBackward[classe1][length_Chain] /= normalisation;
            }
        }


        //Sum joint probabilities a posteriori//
        std::vector< std::vector<double> > Psi(m_Number_classes, std::vector<double>(m_Number_classes, 0.));
        std::vector< std::vector<double> > temp_matrix(m_Number_classes, std::vector<double>(m_Number_classes, 0.));

        for(length_Chain=1; length_Chain<Chain.size()-1; length_Chain++){

            normalisation = 0.;

            for(mod1=0; mod1<m_Number_images; mod1++){
                data[mod1] = Chain[length_Chain+1][mod1];
            }
            
            for(classe1=0; classe1<m_Number_classes; classe1++){

                for(mod1=0; mod1<m_Number_images; mod1++){
                    moyenne_c[mod1] = Moyenne[mod1][classe1];
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
                    }
                }

                for(classe2=0; classe2<m_Number_classes; classe2++){
                    if(Chain_outliers[length_Chain+1]==0){
                        temp_matrix[classe1][classe2] = Forward[classe2][length_Chain] * aij[classe2][classe1] * Chain_Atlas[length_Chain+1][classe1] * Backward[classe1][length_Chain+1] * Gaussian(data, moyenne_c, variance_c); 
                    }else{
                        temp_matrix[classe1][classe2] = Forward[classe2][length_Chain] * aij[classe2][classe1] * Chain_Atlas[length_Chain+1][classe1] * Backward[classe1][length_Chain+1];
                    }
                    normalisation += temp_matrix[classe1][classe2];
                }

            }

            //Normalisation
            if(normalisation!=0.){
                for(classe1=0; classe1<m_Number_classes; classe1++){
                    for(classe2=0; classe2<m_Number_classes; classe2++){    
                        Psi[classe1][classe2] += temp_matrix[classe1][classe2]/normalisation;
                    }
                }
            }else{
                std::cout<<"Division by zero in Psi estimation"<<std::endl;
            } 
                
            

        }
        

        //Re-estimation of prior//
/*        std::cout<<"Proba_ini = "<<std::endl;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            Proba_ini[classe1] = ForwardBackward[classe1][0];
            std::cout<<"\t"<<Proba_ini[classe1]<<std::endl;
        }*/

        //Re-estimation of transition matrix//
        for(classe1=0; classe1<m_Number_classes; classe1++){
            temp = 0.;
            normalisation = 0.;
            for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
                temp += ForwardBackward[classe1][length_Chain];
            }
            for(classe2=0; classe2<m_Number_classes; classe2++){
                if(temp == 0.){
                    aij[classe1][classe2] = 0.;
                }else{
                    aij[classe1][classe2] = Psi[classe1][classe2] / temp;
                    normalisation += aij[classe1][classe2];
                }
            }
            for(classe2=0; classe2<m_Number_classes; classe2++){
                if(normalisation == 0.){
                    aij[classe1][classe2] = 1. / m_Number_classes;
                }else{
                    aij[classe1][classe2] /= normalisation;
                }
            }
        }
 /*       std::cout<<"aij = "<<std::endl;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            for(classe2=0; classe2<m_Number_classes; classe2++){
                std::cout<<"\t"<<aij[classe1][classe2];
            }
            std::cout<<std::endl;
        }*/

        //Re-estimation of means//
        for(mod1=0; mod1<m_Number_images; mod1++){
            for(classe1=0; classe1<m_Number_classes; classe1++){
                temp = 0.;
                normalisation = 0.;
                for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
                    if(Chain_outliers[length_Chain]==0){
                        temp += Chain[length_Chain][mod1] * ForwardBackward[classe1][length_Chain];
                        normalisation += ForwardBackward[classe1][length_Chain];
                    }
                }
                if(normalisation == 0.){
                    Moyenne[mod1][classe1] = 0.;
                }else{
                    Moyenne[mod1][classe1] = temp / normalisation;
                }
            }
        }
        
/*        std::cout<<"Moyenne = "<<std::endl;
        for(mod1=0; mod1<m_Number_images; mod1++){ 
            for(classe1=0; classe1<m_Number_classes; classe1++){
                std::cout<<"\t"<<Moyenne[mod1][classe1];
            }
            std::cout<<std::endl;
        }*/

        
        //Re-estimation of variance//
        for(classe1=0; classe1<m_Number_classes; classe1++){
            normalisation = 0.;
            for(mod1=0; mod1<m_Number_images; mod1++){
                for(mod2=0; mod2<m_Number_images; mod2++){
                    Variance[classe1][mod1][mod2] = 0.;
                }
            }
            for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
                normalisation += ForwardBackward[classe1][length_Chain];
                for(mod1=0; mod1<m_Number_images; mod1++){
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        Variance[classe1][mod1][mod2] += ForwardBackward[classe1][length_Chain] * (Chain[length_Chain][mod1] - Moyenne[mod1][classe1]) * (Chain[length_Chain][mod2] - Moyenne[mod2][classe1]);
                    }
                }
            }
            if(normalisation != 0.){
                for(mod1=0; mod1<m_Number_images; mod1++){
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        Variance[classe1][mod1][mod2] /= normalisation;
                    }
                }
            }
        }

 /*       std::cout<<"Variance = "<<std::endl; 
        for(classe1=0; classe1<m_Number_classes; classe1++){
            std::cout<<"\tClasse "<<classe1+1<<std::endl;
            for(mod1=0; mod1<m_Number_images; mod1++){
                std::cout<<"\t\t"; 
                for(mod2=0; mod2<m_Number_images; mod2++){
                    std::cout<<Variance[classe1][mod1][mod2]<<"\t";
                }
                std::cout<<std::endl;
            }
        }*/


        //Estimation of residus
        std::vector<double> Residu(Chain.size(), 0.);
        double ProbaMin, ProbaMax;
        std::vector<double> PX(m_Number_classes, 0.);
        std::vector<double> PXtemp(m_Number_classes, 0.);

        temp = 0.;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            if(m_Flair_bool){
                data_s[0] = Chain[0][m_Position_Flair];
                moyenne_c_s[0] = Moyenne[m_Position_Flair][classe1];
                variance_c_s[0][0] = Variance[classe1][m_Position_Flair][m_Position_Flair];
                temp += Proba_ini[classe1] * Gaussian(data_s, moyenne_c_s, variance_c_s) * Chain_Atlas[0][classe1];
            }else{
                for(mod1=0; mod1<m_Number_images; mod1++){
                    data[mod1] = Chain[0][mod1];
                    moyenne_c[mod1] = Moyenne[mod1][classe1];
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
                    }
                }
                temp += Proba_ini[classe1] * Gaussian(data, moyenne_c, variance_c) * Chain_Atlas[0][classe1];
            }
            PX[classe1] = Proba_ini[classe1];
        }

        Residu[0] = temp;
        ProbaMin = temp;
        ProbaMax = temp;

        for(length_Chain=1; length_Chain<Chain.size(); length_Chain++){

            for(classe1=0; classe1<m_Number_classes; classe1++){
                temp = 0.;
                for(classe2=0; classe2<m_Number_classes; classe2++){
                    temp += PX[classe2] * aij[classe2][classe1];
                }
                PXtemp[classe1] = temp;
            }

            temp = 0.;
            for(classe1=0; classe1<m_Number_classes; classe1++){
                if(m_Flair_bool){
                    data_s[0] = Chain[length_Chain][m_Position_Flair];
                    moyenne_c_s[0] = Moyenne[m_Position_Flair][classe1];
                    variance_c_s[0][0] = Variance[classe1][m_Position_Flair][m_Position_Flair];
                    PX[classe1] = PXtemp[classe1];
                    temp += PX[classe1] * Gaussian(data_s, moyenne_c_s, variance_c_s) * Chain_Atlas[length_Chain][classe1];
                }
                else{
                    for(mod1=0; mod1<m_Number_images; mod1++){
                        data[mod1] = Chain[length_Chain][mod1];
                        moyenne_c[mod1] = Moyenne[mod1][classe1];
                        for(mod2=0; mod2<m_Number_images; mod2++){
                            variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
                        }
                    }
                    PX[classe1] = PXtemp[classe1];
                    temp += PX[classe1] * Gaussian(data, moyenne_c, variance_c) * Chain_Atlas[length_Chain][classe1];
                }
            }

            Residu[length_Chain] = temp;
            
            if(temp<ProbaMin){
                ProbaMin = temp;
            }
            if(temp>ProbaMax){
                ProbaMax = temp;
            }            

        }


        for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
            Residu[length_Chain] = (Residu[length_Chain] - ProbaMin) / (ProbaMax - ProbaMin);
        }


        //Threshold estimation//
        double Threshold;
        if(m_Criterion_outliers==0){
            unsigned int NbOutliers = (unsigned int) Chain.size() * m_Criterion_outliers_value / 100;
            std::vector<double> temp_vector(Residu);
            std::stable_sort(temp_vector.begin(),temp_vector.end());
            Threshold = temp_vector[NbOutliers];
        }else{
            Threshold = m_Criterion_outliers_value;
        }
        //std::cout<<"Threshold = "<<Threshold<<std::endl;

        for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
            if((Residu[length_Chain] < Threshold) && (Chain[length_Chain][m_Position_Flair] > (Moyenne[m_Position_Flair][0] + 3 * sqrt(Variance[0][m_Position_Flair][m_Position_Flair])))){
                Chain_outliers[length_Chain]=1;
            }else{
                Chain_outliers[length_Chain]=0;
            }
        }

        normalisation = 0.;
        temp = 0.;

        if(iter>0){
            for(classe1=0; classe1<m_Number_classes; classe1++){
                temp += std::abs(Proba_ini[classe1] - Proba_ini2[classe1]);
                normalisation += Proba_ini2[classe1];
                for(classe2=0; classe2<m_Number_classes; classe2++){
                    temp += std::abs(aij[classe1][classe2] - aij2[classe1][classe2]);
                    normalisation += aij2[classe1][classe2];
                }
                for(mod1=0; mod1<m_Number_images; mod1++){
                    temp += std::abs(Moyenne[mod1][classe1] - Moyenne2[mod1][classe1]);
                    normalisation += Moyenne2[mod1][classe1];
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        temp += std::abs(Variance[classe1][mod1][mod2] - Variance2[classe1][mod1][mod2]);
                        normalisation += Variance2[classe1][mod1][mod2];
                    }
                }
            }
            change = temp / normalisation;
        }

        Moyenne2 = Moyenne;
        Variance2 = Variance;
        aij2 = aij;
        Proba_ini2 = Proba_ini;
        
        std::cout<<"\tchange = "<<change<<std::endl;

        iter++;

    }
}


template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::Chain_MPM(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Chain_Atlas, std::vector< std::vector<double> > &Moyenne, std::vector< std::vector< std::vector<double> > > &Variance, std::vector<double> &Proba_ini,  std::vector< std::vector<double> > &aij, std::vector<unsigned int> &Chain_seg, std::vector<unsigned int> &Chain_outliers)
{
    unsigned int length_Chain, classe1, classe2, mod1, mod2;
    double normalisation, temp;

    std::vector<double> data(m_Number_images, 0.);
    std::vector<double> moyenne_c(m_Number_images, 0.);
    std::vector< std::vector<double> > variance_c(m_Number_images, std::vector<double>(m_Number_images, 0.));

    std::vector< std::vector<double> > Forward(m_Number_classes, std::vector<double> (Chain.size(), 0.));
    std::vector< std::vector<double> > Backward(m_Number_classes, std::vector<double> (Chain.size(), 0.));
    std::vector< std::vector<double> > ForwardBackward(m_Number_classes, std::vector<double> (Chain.size(), 0.));

    //Forward probabilities estimation//
    normalisation = 0.;

    //length_Chain=0
    for(classe1=0; classe1<m_Number_classes; classe1++){

        for(mod1=0; mod1<m_Number_images; mod1++){
            data[mod1] = Chain[0][mod1];
            moyenne_c[mod1] = Moyenne[mod1][classe1];
            for(mod2=0; mod2<m_Number_images; mod2++){
                variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
            }
        }

        if(Chain_outliers[0]==0){
            Forward[classe1][0] = Proba_ini[classe1] * Chain_Atlas[0][classe1] * Gaussian(data, moyenne_c, variance_c);
            normalisation += Forward[classe1][0];
        }else{
            Forward[classe1][0] = Proba_ini[classe1] * Chain_Atlas[0][classe1];
            normalisation += Forward[classe1][0];
        }
    }

    //Normalisation
    for(classe1=0; classe1<m_Number_classes; classe1++){
        if(normalisation!=0.){
            Forward[classe1][0] /= normalisation;
        }else{
            std::cout<<"Division by zero in Forward estimation"<<std::endl;
            Forward[classe1][0] = 1. / m_Number_classes;
        }
    }

    //length_Chain!=0
    for(length_Chain=1; length_Chain<Chain.size(); length_Chain++){

        normalisation = 0.;
        
        for(classe1=0; classe1<m_Number_classes; classe1++){
            temp = 0.;

            for(mod1=0; mod1<m_Number_images; mod1++){
                data[mod1] = Chain[length_Chain][mod1];
                moyenne_c[mod1] = Moyenne[mod1][classe1];
                for(mod2=0; mod2<m_Number_images; mod2++){
                    variance_c[mod1][mod2] = Variance[classe1][mod1][mod2];
                }
            }

            for(classe2=0; classe2<m_Number_classes; classe2++){
                temp += Forward[classe2][length_Chain-1] * aij[classe2][classe1];
            }

            if(Chain_outliers[length_Chain]==0){
                Forward[classe1][length_Chain] = temp * Chain_Atlas[length_Chain][classe1] * Gaussian(data, moyenne_c, variance_c);
            }else{
                Forward[classe1][length_Chain] = temp * Chain_Atlas[length_Chain][classe1];
            }

            normalisation += Forward[classe1][length_Chain];
        }

        //Normalisation
        for(classe1=0; classe1<m_Number_classes; classe1++){
            if(normalisation!=0.){
                Forward[classe1][length_Chain] /= normalisation;
            }else{
                std::cout<<"Division by zero in Forward estimation"<<std::endl;
                Forward[classe1][length_Chain] = 1. / m_Number_classes;
            }
        }
    }


    //Backward probabilities estimation//
    normalisation = 0.;

    //length_Chain = end
    for(classe1=0; classe1<m_Number_classes; classe1++){
        Backward[classe1][Chain.size()-1] = 1.;
    }


    //length_Chain != end
    for(int length_Chain=Chain.size()-2; length_Chain>=0; length_Chain--){

        normalisation = 0.;

        for(mod1=0; mod1<m_Number_images; mod1++){
            data[mod1] = Chain[length_Chain+1][mod1];
        }
        
        for(classe1=0; classe1<m_Number_classes; classe1++){
            temp = 0.;

            for(classe2=0; classe2<m_Number_classes; classe2++){

                for(mod1=0; mod1<m_Number_images; mod1++){
                    moyenne_c[mod1] = Moyenne[mod1][classe2];
                    for(mod2=0; mod2<m_Number_images; mod2++){
                        variance_c[mod1][mod2] = Variance[classe2][mod1][mod2];
                    }
                }

                if(Chain_outliers[length_Chain+1]==0){
                    temp += Backward[classe2][length_Chain+1] * Chain_Atlas[length_Chain+1][classe2] * Gaussian(data, moyenne_c, variance_c) * aij[classe1][classe2];
                }else{
                    temp += Backward[classe2][length_Chain+1] * Chain_Atlas[length_Chain+1][classe2] * aij[classe1][classe2];
                }

            }
            
            Backward[classe1][length_Chain] = temp;
            normalisation += temp;               
            
        }

        //Normalisation
        for(classe1=0; classe1<m_Number_classes; classe1++){
            if(normalisation!=0.){
                Backward[classe1][length_Chain] /= normalisation;
            }else{
                std::cout<<"Division by zero in Forward estimation"<<std::endl;
                Backward[classe1][length_Chain] = 1. / sqrt(m_Number_classes);
            }
        }
    }


    //Forward-Backward estimation, proba marginales a posterioiri//
    for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
        normalisation = 0.;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            ForwardBackward[classe1][length_Chain] = Forward[classe1][length_Chain] * Backward[classe1][length_Chain];
            normalisation += ForwardBackward[classe1][length_Chain];
        }
        for(classe1=0; classe1<m_Number_classes; classe1++){
            ForwardBackward[classe1][length_Chain] /= normalisation;
        }
    }

    //MPM estimation//
    double max;
    unsigned int imax;
    for(length_Chain=0; length_Chain<Chain.size(); length_Chain++){
        max = 0.;
        imax = 0;
        for(classe1=0; classe1<m_Number_classes; classe1++){
            if(ForwardBackward[classe1][length_Chain]>max){
                max = ForwardBackward[classe1][length_Chain];
                imax = classe1+1;
            }
        }
        Chain_seg[length_Chain] = imax;
    }
    
}


template<typename TInputImage, typename TOutputImage>
void
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::Images_Reconstruction(std::vector<unsigned int> &Chain_seg, std::vector<unsigned int> &Chain_outliers, ConstIteratorType *mask_it, IteratorType &output_seg_it, IteratorType &output_outliers_it)
{
    mask_it->GoToBegin();
    output_seg_it.GoToBegin();
    output_outliers_it.GoToBegin();
    unsigned int i=0;
    while(!mask_it->IsAtEnd()){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            output_seg_it.Set(Chain_seg[i]);
            output_outliers_it.Set(Chain_outliers[i]);
            i++;
        }else{
            output_seg_it.Set(0);
            output_outliers_it.Set(0);
        }
        ++(*mask_it);
        ++output_seg_it;
        ++output_outliers_it;
    }
}


template<typename TInputImage, typename TOutputImage>
double
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::Gaussian(std::vector<double> Data, std::vector<double> Mu, std::vector< std::vector<double> > Sigma)
{
    double res=0.;
    double DetSigma;
    double arg;

    std::vector< std::vector<double> > InvSigma(Sigma);

    switch(Mu.size()){
        case 1:
            res=(1./(sqrt(2*math_Pi)*sqrt(Sigma[0][0])))*exp((Data[0]-Mu[0])*(Data[0]-Mu[0])/(-2*Sigma[0][0]));
            break;
  
        case 2:
            DetSigma=Sigma[0][0]*Sigma[1][1]-Sigma[1][0]*Sigma[0][1];
  
            if(DetSigma!=0.)
            {
                InvSigma[0][0]=Sigma[1][1]/DetSigma;
                InvSigma[0][1]=-Sigma[0][1]/DetSigma;
                InvSigma[1][0]=-Sigma[1][0]/DetSigma;
                InvSigma[1][1]=Sigma[0][0]/DetSigma;

                arg=(Data[0]-Mu[0])*((Data[0]-Mu[0])*InvSigma[0][0]+(Data[1]-Mu[1])*InvSigma[1][0])+(Data[1]-Mu[1])*((Data[0]-Mu[0])*InvSigma[0][1]+(Data[1]-Mu[1])*InvSigma[1][1]);

                res=(1./(pow(2*math_Pi,(double)Data.size()/2)*sqrt(fabs(DetSigma))))*exp(-0.5*arg);
            }else
                res=-1.;
            break;

        case 3:
            DetSigma=Sigma[0][0]*(Sigma[1][1]*Sigma[2][2]-Sigma[2][1]*Sigma[1][2]);
            DetSigma+=-Sigma[1][0]*(Sigma[0][1]*Sigma[2][2]-Sigma[2][1]*Sigma[0][2]);
            DetSigma+=Sigma[2][0]*(Sigma[0][1]*Sigma[1][2]-Sigma[1][1]*Sigma[0][2]);

            if(DetSigma!=0.)
            {
                InvSigma[0][0]=(Sigma[1][1]*Sigma[2][2]-Sigma[1][2]*Sigma[2][1])/DetSigma;
                InvSigma[0][1]=(Sigma[0][2]*Sigma[2][1]-Sigma[0][1]*Sigma[2][2])/DetSigma;
                InvSigma[0][2]=(Sigma[0][1]*Sigma[1][2]-Sigma[0][2]*Sigma[1][1])/DetSigma;

                InvSigma[1][0]=(Sigma[1][2]*Sigma[2][0]-Sigma[1][0]*Sigma[2][2])/DetSigma;
                InvSigma[1][1]=(Sigma[0][0]*Sigma[2][2]-Sigma[0][2]*Sigma[2][0])/DetSigma;
                InvSigma[1][2]=(Sigma[0][2]*Sigma[1][0]-Sigma[0][0]*Sigma[1][2])/DetSigma;

                InvSigma[2][0]=(Sigma[1][0]*Sigma[2][1]-Sigma[1][1]*Sigma[2][0])/DetSigma;
                InvSigma[2][1]=(Sigma[0][1]*Sigma[2][0]-Sigma[0][0]*Sigma[2][1])/DetSigma;
                InvSigma[2][2]=(Sigma[0][0]*Sigma[1][1]-Sigma[0][1]*Sigma[1][0])/DetSigma;

                arg=(Data[0]-Mu[0])*((Data[0]-Mu[0])*InvSigma[0][0]+(Data[1]-Mu[1])*InvSigma[1][0]+(Data[2]-Mu[2])*InvSigma[2][0]);
                arg+=(Data[1]-Mu[1])*((Data[0]-Mu[0])*InvSigma[0][1]+(Data[1]-Mu[1])*InvSigma[1][1]+(Data[2]-Mu[2])*InvSigma[2][1]);
                arg+=(Data[2]-Mu[2])*((Data[0]-Mu[0])*InvSigma[0][2]+(Data[1]-Mu[1])*InvSigma[1][2]+(Data[2]-Mu[2])*InvSigma[2][2]);

                res=(1./(pow(2*math_Pi,(double)Data.size()/2)*sqrt(fabs(DetSigma))))*exp(-0.5*arg);
            }else
                res=-1.;
            break;

        default:
            std::cout<<"Not implemented"<<std::endl;
            break;
    }
    return res;
}

}

#endif //_itkHiddenMarkovChainFilter_txx
