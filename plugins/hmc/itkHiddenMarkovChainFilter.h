#ifndef _itkHiddenMarkovChainFilter_h
#define _itkHiddenMarkovChainFilter_h

#include <itkImageToImageFilter.h>
#include "itkHilbertPath.h"
#include <itkPathIterator.h>
#include <itkPathConstIterator.h>

#include<vector>


namespace itk{

template<typename TInputImage, typename TOutputImage>
class HiddenMarkovChainFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public :
    /** Standard class typedefs. */
    typedef HiddenMarkovChainFilter Self;
    typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
    typedef SmartPointer<Self> Pointer;
    typedef SmartPointer<Self const> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(HiddenMarkovChainFilter, ImageToImageFilter);

    /** Type of the input image. */
    typedef TInputImage InputImageType;

    /** Type of the output image. */
    typedef TOutputImage OutputImageType;
    typedef typename OutputImageType::IndexType OutputImageIndexType;


    typedef HilbertPath<unsigned int, 3> PathType;
    typedef typename PathType::HilbertOrderType PathOrderType;
    typedef PathConstIterator<InputImageType, PathType> ConstIteratorType;
    typedef PathIterator<OutputImageType, PathType> IteratorType;

    /** Number of input images*/
    itkGetMacro(Number_images, unsigned int);
    itkSetMacro(Number_images, unsigned int);

    /** Boolean for atlas*/
    itkGetMacro(Atlas_bool, bool);
    itkSetMacro(Atlas_bool, bool);

    /** Number of iterations*/
    itkGetMacro(Number_iter, unsigned int);
    itkSetMacro(Number_iter, unsigned int);

    /** Number of classes*/
    itkGetMacro(Number_classes, unsigned int);
    itkSetMacro(Number_classes, unsigned int);

    /** Type of criterion for outliers (0 for percentage and 1 for threshold)*/
    itkGetMacro(Criterion_outliers, unsigned int);
    itkSetMacro(Criterion_outliers, unsigned int);

    /** Value for the criterion for outliers*/
    itkGetMacro(Criterion_outliers_value, double);
    itkSetMacro(Criterion_outliers_value, double);

    /** Taille cube*/
    itkGetMacro(Taille_cube, PathOrderType);
    itkSetMacro(Taille_cube, PathOrderType);

    /** Outliers only on Flair (0 for all images and 1 for Flair only)*/
    itkGetMacro(Flair_bool, bool);
    itkSetMacro(Flair_bool, bool);

    /** Position of the Flair*/
    itkGetMacro(Position_Flair, unsigned int);
    itkSetMacro(Position_Flair, unsigned int);

    TOutputImage* GetOuputSegImage();
    TOutputImage* GetOuputOutliersImage();

    typename TInputImage::Pointer GetInputImage(unsigned int idx);
    void SetInputImage(unsigned int idx, TInputImage* image);

protected :
    HiddenMarkovChainFilter();
    ~HiddenMarkovChainFilter();
    void PrintSelf(std::ostream& os, itk::Indent indent) const;
    void GenerateData();

    /**  Create the Output */
    DataObject::Pointer MakeOutput(unsigned int idx);


private :

    unsigned int m_Number_images;
    bool m_Atlas_bool;
    unsigned int m_Number_iter;
    unsigned int m_Number_classes;
    unsigned int m_Criterion_outliers;
    double m_Criterion_outliers_value;
    PathOrderType m_Taille_cube;
    bool m_Flair_bool;
    unsigned int m_Position_Flair;



    HiddenMarkovChainFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    void Chains_Creation(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Chain_Atlas, ConstIteratorType *mask_it, std::vector<ConstIteratorType *> &inputs_its, std::vector<ConstIteratorType *> &atlas_its);
    void Chain_Itialisation(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Moyenne, std::vector< std::vector< std::vector<double> > > &Variance, std::vector<double> &Proba_ini,  std::vector< std::vector<double> > &aij);
    void Chain_EM(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Chain_Atlas, std::vector< std::vector<double> > &Moyenne, std::vector< std::vector< std::vector<double> > > &Variance, std::vector<double> &Proba_ini,  std::vector< std::vector<double> > &aij, std::vector<unsigned int> &Chain_seg , std::vector<unsigned int> &Chain_outliers);
    void Chain_MPM(std::vector< std::vector<double> > &Chain, std::vector< std::vector<double> > &Chain_Atlas, std::vector< std::vector<double> > &Moyenne, std::vector< std::vector< std::vector<double> > > &Variance, std::vector<double> &Proba_ini,  std::vector< std::vector<double> > &aij, std::vector<unsigned int> &Chain_seg, std::vector<unsigned int> &Chain_outliers);
    void Images_Reconstruction(std::vector<unsigned int> &Chain_seg, std::vector<unsigned int> &Chain_outliers, ConstIteratorType *mask_it, IteratorType &output_seg_it, IteratorType &output_outliers_it);
    double Gaussian(std::vector<double> Data, std::vector<double> Mu, std::vector< std::vector<double> > Sigma);
};

}

#include "itkHiddenMarkovChainFilter.txx"

#endif //_itkHiddenMarkovChainFilter_h
