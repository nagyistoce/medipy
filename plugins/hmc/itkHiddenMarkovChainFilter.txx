#ifndef _itkHiddenMarkovChainFilter_txx
#define _itkHiddenMarkovChainFilter_txx

#include "itkImageToImageFilter.h"


#include<vector>


#include <itkObjectFactory.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkPathIterator.h>
#include <itkPathConstIterator.h>




namespace itk{


template<typename TInputImage, typename TOutputImage>
HiddenMarkovChainFilter<TInputImage, TOutputImage>
::HiddenMarkovChainFilter()
:m_Nb_images(0), m_Atlas_bool(true), m_Number_iter(5), m_Number_classes(3), m_Criterion_outliers(1), m_Criterion_outliers_value(0.05)
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
    os << indent << "Number of input images : " << this->m_Nb_images << "\n";
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
    typename PathType::Pointer hilbert_path = PathType::New();
    hilbert_path->SetHilbertOrder(m_Taille_cube);
    hilbert_path->Initialize();
    
    std::vector<ConstIteratorType *> inputs_its;
    std::vector<ConstIteratorType *> atlas_its;

    ConstIteratorType * mask_it;

    for(unsigned int i=0; i < m_Nb_images; i++){
        typename InputImageType::Pointer input = this->GetInputImage(i);
        ConstIteratorType *it = new ConstIteratorType(input, hilbert_path);
        it->GoToBegin();
        inputs_its.push_back(it);
    }

    if(m_Atlas_bool){
        for(unsigned int i=m_Nb_images; i < m_Nb_images+3; i++){
            typename InputImageType::Pointer input = this->GetInputImage(i);
            ConstIteratorType *it = new ConstIteratorType(input, hilbert_path);
            it->GoToBegin();
            inputs_its.push_back(it);
        }
        mask_it = new ConstIteratorType(this->GetInputImage(m_Nb_images+3), hilbert_path);
        mask_it->GoToBegin();
    }else{
        mask_it = new ConstIteratorType(this->GetInputImage(m_Nb_images), hilbert_path);
        mask_it->GoToBegin();
    }

    /*

    //Initialisation
    std::vector<float> PIi(m_Number_classes, 1./m_Number_classes);
    std::vector< std::vector<float> > aij(m_Number_classes, std::vector<float>(m_Number_classes, 1./(4*(m_Number_classes-1))));
    for(unsigned int i=0; i<m_Number_classes; i++){
        aij[i][i]=3./4;    
    }

    typename InputImageType::PixelType min, max;
    unsigned int quantification=128;
    float hist, intervalle, S, ecart, Tinter;
    double math_Pi=3.141592654;
    unsigned int taille, res, res2, i, k;
    
    min=NumericTraits<typename InputImageType::PixelType>::Zero;
    max=NumericTraits<typename InputImageType::PixelType>::Zero;
    taille=0;

    ConstIteratorType * input_it = inputs_its[0];
    input_it->GoToBegin();
    mask_it->GoToBegin();

    //recherche de min et max
    while(!(input_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            if(input_it->Get()>max){max=input_it->Get();}
            if(input_it->Get()<min){min=input_it->Get();}
            taille++;
        }
        ++(*input_it);
        ++(*mask_it);   
    }

    std::vector<float> histo(quantification, 0.);
    input_it->GoToBegin();
    mask_it->GoToBegin();
    while(!(input_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            hist=(input_it->Get()-min)*((float)quantification-1)/(max-min);
            res=(int)rint(hist);
            histo[res]++;
        }
        ++(*input_it);
        ++(*mask_it);         
    }

    for(i=0; i!=quantification; i++){
        histo[i]/=taille;
    }

    intervalle=1./(2*m_Number_classes);

    std::vector<float> Pi(m_Number_classes, 0.);
    std::vector<float> Mu1(m_Number_classes, 0.);
    std::vector<float> Sigma1(m_Number_classes, 0.);

    //moyennes d'initialisation
    for(k=0; k!=m_Number_classes; k++){
        S=0.;
        i=0;
        while(S<(intervalle*(2*k+1))){
            S+=histo[i];
            i++;
        }
        Mu1[k]=(float)i;
    }

    //variances d'initialisation
    res=0;
    for(k=0; k!=m_Number_classes-1; k++){
        Tinter=0.;
        ecart=0.;
        res2=(int)(rint((Mu1[k+1]-Mu1[k])/2)+Mu1[k]);

        for(i=res; i!=res2; i++){
            ecart+=(i-Mu1[k])*(i-Mu1[k])*histo[i];
            Tinter+=histo[i];
        }
        if((Tinter!=0)&&(ecart!=0)){
            Pi[k]=Tinter;
            Sigma1[k]=sqrt(ecart);
        }else{
            Pi[k]=1/(m_Number_classes+1);
            Sigma1[k]=1/sqrt(2*math_Pi);
        }
        res=res2;
    }

    //variance dernier intervalle
    Tinter=0.;
    ecart=0.;
    for(i=res; i!=quantification; i++){
        ecart+=(i-Mu1[k])*(i-Mu1[k])*histo[i];
        Tinter+=histo[i];
    }
    if((Tinter!=0)&&(ecart!=0)){
        Pi[m_Number_classes-1]=Tinter;
        Sigma1[m_Number_classes-1]=sqrt(ecart);
    }else{
        Pi[m_Number_classes-1]=1/(m_Number_classes+1);
        Sigma1[m_Number_classes-1]=1/sqrt(2*math_Pi);
    }

    //histo EM
    float err=1000., eps=1.e-45;
    int iter=0;
    std::vector< std::vector<float> > apriori(quantification, std::vector<float>(m_Number_classes, 0.));
    std::vector<float> aprioriDen(quantification, 0.);
    std::vector<float> Pi_new(m_Number_classes, 0.);
    std::vector<float> Mu_new(m_Number_classes, 0.);
    std::vector<float> Sigma_new(m_Number_classes, 0.);


    for(k=0; k!=m_Number_classes; k++){
        for(i=0; i!=quantification; i++){
            apriori[i][k]=Pi[k]*(1./(sqrt(2*math_Pi)*Sigma1[k]))*exp(((float)i-Mu1[k])*((float)i-Mu1[k])/(-2*Sigma1[k]*Sigma1[k]));
        }
    }

    while((err>1e-5)&&(iter<1000)){
        iter++;
        for(k=0; k!=m_Number_classes; k++){
            for(i=0; i!=quantification; i++){
                apriori[i][k]=Pi[k]*(1./(sqrt(2*math_Pi)*Sigma1[k]))*exp(((float)i-Mu1[k])*((float)i-Mu1[k])/(-2*Sigma1[k]*Sigma1[k]));
            }
        }
        S=0.;

        for(i=0; i!=quantification; i++){
            for(k=0; k!=m_Number_classes; k++)
                aprioriDen[i]+=apriori[i][k];
            S+=histo[i];
        }

        for(k=0; k!=m_Number_classes; k++)
            for(i=0; i!=quantification; i++){
                if((aprioriDen[i]>eps)&&(apriori[i][k]>eps))
                    apriori[i][k]/=aprioriDen[i];
                else
                    apriori[i][k]=0;
        }

        //nouvelle proba a priori
        for(k=0; k!=m_Number_classes; k++){
            for(i=0; i!=quantification; i++)
                Pi_new[k]+=histo[i]*apriori[i][k];

            if(Pi_new[k]!=0)
                Pi_new[k]/=S;
            else
                Pi_new[k]=1/(m_Number_classes+1);
        }

        //nouvelle moyenne
        for(k=0; k!=m_Number_classes; k++){
            S=0.;
            for(i=0; i!=quantification; i++){
                Mu_new[k]+=histo[i]*i*apriori[i][k];
                S+=histo[i]*apriori[i][k];
            }
            Mu_new[k]/=S;
        }

        //nouvel ecart-type
        for(k=0; k!=m_Number_classes; k++){
            S=0.;
            for(i=0;i!=quantification;i++){
                Sigma_new[k]+=histo[i]*(i-Mu_new[k])*(i-Mu_new[k])*apriori[i][k];
                S+=histo[i]*apriori[i][k];
            }
            if((Sigma_new[k]/S)>(1/(2*math_Pi)))
                Sigma_new[k]=sqrt(Sigma_new[k]/S);
            else
                Sigma_new[k]=1/(sqrt(2*math_Pi));
	    }

        S=0.;
        err=0.;

        //calcul de l'erreur
        for(k=0; k!=m_Number_classes; k++){
            err+=(Pi[k]-Pi_new[k])*(Pi[k]-Pi_new[k])+(Mu1[k]-Mu_new[k])*(Mu1[k]-Mu_new[k])+(Sigma1[k]-Sigma_new[k])*(Sigma1[k]-Sigma_new[k]);
            S+=Pi[k]+Mu1[k]+Sigma1[k];
        }

        err=sqrt(err/S);

        //passage des nouveaux parametres
        Pi=Pi_new;
        Mu1=Mu_new;
        Sigma1=Sigma_new;

        aprioriDen.clear();
        Pi_new.clear();
        Mu_new.clear();
        Sigma_new.clear();

    }

    apriori.clear();

    //transformation des valeurs pour qu'elles correspondent à la réalité
    for(k=0; k!=m_Number_classes; k++){
        Mu1[k]=min+Mu1[k]*(float)(max-min)/quantification;
        Sigma1[k]=Sigma1[k]*(float)(max-min)/quantification;
    }
    histo.clear();

    unsigned int TailleEch=8;
    float Max_V=0.;

    std::vector<unsigned int> MP(taille);
    std::vector<unsigned int> NbVox(m_Number_classes);
    std::vector< std::vector<double> > MuProv(m_Nb_images, std::vector<double>(m_Number_classes,0.));
    std::vector< std::vector<double> > SigmaProv(m_Nb_images, std::vector<double>(m_Number_classes,0.));


    //carte de segmentation issue de l'EM
    input_it->GoToBegin();
    mask_it->GoToBegin();
    i=0;
    while(!(input_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            Max_V=-1;
            MP[i]=-1;
            for(k=0; k!=m_Number_classes;k++){
                if((1./(sqrt(2*math_Pi)*Sigma1[k]))*exp((input_it->Get()-Mu1[k])*(input_it->Get()-Mu1[k])/(-2*Sigma1[k]*Sigma1[k]))>Max_V){
                    Max_V=(1./(sqrt(2*math_Pi)*Sigma1[k]))*exp((input_it->Get()-Mu1[k])*(input_it->Get()-Mu1[k])/(-2*Sigma1[k]*Sigma1[k]));
                    MP[i]=k;
                }
            }
            i++;
        }
        ++(*input_it);
        ++(*mask_it);           
    }

    //calcul du nb de voxels appartenant à chaque classe
    for(k=0; k!=m_Number_classes; k++){
        S=0;
        for(i=0; i!=taille; i++)
            if(MP[i]==k)
                S++;
        NbVox[k]=S;
    }

    //calcul de la moyenne pour chaque classe
    i=0;
    mask_it->GoToBegin();
    for(unsigned int j=0; j!=m_Nb_images; j++)
        inputs_its[j]->GoToBegin();
    while(!(mask_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            for(unsigned int j=0; j!=m_Nb_images; j++){
                MuProv[j][MP[i]]+=inputs_its[j]->Get();
            }
            i++;
        }
        ++(*mask_it);
        for(unsigned int j=0; j!=m_Nb_images; j++){
            ++(*inputs_its[j]);
        }
    }
    for(k=0; k!=m_Number_classes; k++){
        for(unsigned int j=0; j!=m_Nb_images; j++){
            MuProv[j][k]/=NbVox[k];
        }
    }


    //calcul de l'ecart-type 
    i=0;
    mask_it->GoToBegin();
    for(unsigned int j=0; j!=m_Nb_images; j++)
        inputs_its[j]->GoToBegin();
    while(!(mask_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            for(unsigned int j=0; j!=m_Nb_images; j++){
                SigmaProv[j][MP[i]]+=(inputs_its[j]->Get()-MuProv[j][MP[i]])*(inputs_its[j]->Get()-MuProv[j][MP[i]]);
            }
            i++;
        }
        ++(*mask_it);
        for(unsigned int j=0; j!=m_Nb_images; j++){
            ++(*inputs_its[j]);
        }
    }
    for(k=0; k!=m_Number_classes; k++){
        for(unsigned int j=0; j!=m_Nb_images; j++){
            SigmaProv[j][k]=sqrt(SigmaProv[j][k]/(NbVox[k]-1.));
        }
    }


/*
    std::cout<<"MuProv"<<std::endl;
    for(unsigned int j=0; j!=m_Nb_images; j++){
        for(k=0; k!=m_Number_classes; k++){
            std::cout<<MuProv[j][k]<<"\t";
        }
        std::cout<<std::endl;
    }

    std::cout<<"SigmaProv"<<std::endl;
    for(unsigned int j=0; j!=m_Nb_images; j++){
        for(k=0; k!=m_Number_classes; k++){
            std::cout<<SigmaProv[j][k]<<"\t";
        }
        std::cout<<std::endl;
    }*/


/*
//Kmoyennes
    int Ini=0;
    int m, par;
    float mini;
    //int Taille=Chain.size2();
    int n=taille/TailleEch;
    //int r=Taille%TailleEch;
    float *Norme = new float[2];
    std::vector< std::vector<double> > C(m_Number_classes, std::vector<double>(2*m_Nb_images,0.));
    std::vector< std::vector<double> > y(m_Nb_images, std::vector<double>(TailleEch,0.));
    std::vector<double> mu_prov(m_Nb_images,0.);
    std::vector<double> Sigma_prov(m_Nb_images,0.);
    std::vector< std::vector<double> > c;//(n, std::vector<double>(2*m_Nb_images,0.));
    std::vector<double> tab(m_Number_classes);
    std::vector<int> ind(n);
    std::vector<int> card(m_Number_classes);
 
    //on insere les resultats de l'em dans les kmeans comme centroides initiaux
    for(unsigned int k=0;k!=m_Number_classes;k++){
        for(unsigned int j=0;j!=m_Nb_images;j++){
	        C[k][j]=MuProv[j][k];
	        C[k][j+m_Nb_images]=SigmaProv[j][k];
        }
    }
    i=0;
    mask_it->GoToBegin();
    for(unsigned int j=0; j!=m_Nb_images; j++)
        inputs_its[j]->GoToBegin();
    while(!(mask_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            for(unsigned int j=0; j!=m_Nb_images; j++){
                y[j][i]=inputs_its[j]->Get();
                mu_prov[j]+=y[j][i]/((double)TailleEch);
            }
            i++;
        }
        ++(*mask_it);
        for(unsigned int j=0; j!=m_Nb_images; j++){
            ++(*inputs_its[j]);
        }
        if(i==TailleEch){
            i=0;
            for(unsigned int k=0;k!=TailleEch;k++){
                for(unsigned int j=0; j!=m_Nb_images; j++){
                    Sigma_prov[j]+=(y[j][k]-mu_prov[j])*(y[j][k]-mu_prov[j])/((double)TailleEch);
                }
            }
            std::vector<double> temp(mu_prov);
            temp.insert(temp.end(), Sigma_prov.begin(), Sigma_prov.end());
            c.push_back(temp);
            std::fill(mu_prov.begin(), mu_prov.end(), 0.);
            std::fill(Sigma_prov.begin(), Sigma_prov.end(), 0.);
        }
    }
    for(unsigned int k=0;k!=i;k++){
        for(unsigned int j=0; j!=m_Nb_images; j++){
            Sigma_prov[j]+=(y[j][k]-mu_prov[j])*(y[j][k]-mu_prov[j])/((double)TailleEch);
        }
    }
    for(unsigned int j=0; j!=m_Nb_images; j++)
        Sigma_prov[j]=sqrt(Sigma_prov[j]);
    std::vector<double> temp(mu_prov);
    temp.insert(temp.end(), Sigma_prov.begin(), Sigma_prov.end());
    c.push_back(temp);


    Norme[0]=1000.;
    Norme[1]=0.;


    while((fabs(Norme[0]-Norme[1])>1e-30)&&(Ini!=1000)){
        par=Ini%2;

        //ini card
        for(unsigned int k=0;k!=m_Number_classes;k++)
            card[k]=0;
        for(int p=0;p!=n;p++){

            //ini tab
            std::fill(tab.begin(), tab.end(), 0.);

            //calcul tab
	  
            for(unsigned int j=0;j!=m_Number_classes;j++)
                for(unsigned int k=0;k!=m_Nb_images;k++)
                    tab[j]+=(C[j][k]-c[p][k])*(C[j][k]-c[p][k]);
	  
	        for(unsigned int j=0;j!=m_Number_classes;j++)
	            tab[j]=sqrt(tab[j]);
	  
	        mini=tab[0];
	        ind[p]=0;
	  
	        for(unsigned int j=0;j!=m_Number_classes;j++){
	            if(tab[j]<=mini){
		            mini=tab[j];
		            ind[p]=j;
	            }
            }
	 
	        card[ind[p]]++;
	    }

        //ini Mu et Sigma
        for(unsigned int p=0;p!=m_Nb_images;p++){
            for(unsigned int j=0;j!=m_Number_classes;j++){
                MuProv[p][j]=0.;
                SigmaProv[p][j]=0.;
            }
        }

        for(unsigned int j=0;j!=m_Number_classes;j++){
	        std::vector<int> masque(card[j]);
	        m=0;
	        for(int p=0;p!=n;p++){
	            if(ind[p]==j){
                    masque[m]=p;
		            m++;
	            }
            }

	        for(unsigned int p=0;p!=m_Nb_images;p++){
	            for(int k=0;k!=card[j];k++){
		            MuProv[p][j]+=(1./(float)card[j])*c[masque[k]][p];
		            SigmaProv[p][j]+=(1./(float)card[j])*c[masque[k]][p+m_Nb_images];
	            }
            }
            masque.clear();
        }

        //calcul changement  
   
        Norme[par]=0;
        for(unsigned int p=0;p!=m_Nb_images;p++)
	        for(unsigned int j=0;j!=m_Number_classes;j++)
	            Norme[par]+=(MuProv[p][j]-C[j][p])*(MuProv[p][j]-C[j][p])+(SigmaProv[p][j]-C[j][p+m_Nb_images])*(SigmaProv[p][j]-C[j][p+m_Nb_images]);

        Norme[par]=sqrt(Norme[par]);

        for(unsigned int p=0;p!=m_Nb_images;p++)
	        for(unsigned int j=0;j!=m_Number_classes;j++){
	            C[j][p]=MuProv[p][j];
	            C[j][p+m_Nb_images]=SigmaProv[p][j];
	        }
        Ini++;
    }

    for(unsigned int p=0;p!=m_Nb_images;p++)
        for(unsigned int j=0;j!=m_Number_classes;j++){
	        MuProv[p][j]=C[j][p];
	        SigmaProv[p][j]=C[j][p+m_Nb_images];
    }


    std::vector< std::vector< std::vector<double> > > Sigma;
    
    Sigma.resize(m_Number_classes);
    for(unsigned int j=0;j<m_Number_classes;j++){
        Sigma[j].resize(m_Nb_images);
        for(unsigned int p=0;p!=m_Nb_images;p++){
            Sigma[j][p].resize(m_Nb_images);
        }
    }

    for(unsigned int j=0;j!=m_Number_classes;j++)
        for(unsigned int p=0;p!=m_Nb_images;p++)
            Sigma[j][p][p]=SigmaProv[p][j]*SigmaProv[p][j];


    std::cout<<"MuProv"<<std::endl;
    for(unsigned int j=0; j!=m_Nb_images; j++){
        for(k=0; k!=m_Number_classes; k++){
            std::cout<<MuProv[j][k]<<"\t";
        }
        std::cout<<std::endl;
    }

    std::cout<<"SigmaProv"<<std::endl;
    for(unsigned int j=0; j!=m_Nb_images; j++){
        for(k=0; k!=m_Number_classes; k++){
            std::cout<<SigmaProv[j][k]<<"\t";
        }
        std::cout<<std::endl;
    }


*/


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

    //mask_it->GoToBegin();
    while(!output_seg_it.IsAtEnd()){
        output_seg_it.Set(1);
        ++output_seg_it;
        //output_outliers_it.Set(mask_it->Get());
        //++output_outliers_it;
        //++(*mask_it);
    }
    while(!output_outliers_it.IsAtEnd()){
        output_outliers_it.Set(1);
        ++output_outliers_it;
    }
/*
    mask_it->GoToBegin();
    for(unsigned int j=0; j!=m_Nb_images; j++)
        inputs_its[j]->GoToBegin();
    while(!(mask_it->IsAtEnd())){
        if(mask_it->Get()!=NumericTraits<typename InputImageType::PixelType>::Zero){
            double toto=0.;
            double big_toto=1e25;
            for(k=0; k!=m_Number_classes; k++){
                for(unsigned int j=0; j!=m_Nb_images; j++){
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
        for(unsigned int j=0; j!=m_Nb_images; j++){
            ++(*inputs_its[j]);
        }
        ++output_seg_it;
    }*/

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


}

#endif //_itkHiddenMarkovChainFilter_txx
