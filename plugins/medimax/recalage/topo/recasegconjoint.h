/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

#ifndef RECASEGCONJOINT_H
#define RECASEGCONJOINT_H

#include<vector>



#ifdef __cplusplus
extern "C" {
#endif
int show_picture_3d (int pos);
void lancementrecaseg();
void matching_Bspline_topo_primitive_non_seg_3d();

#ifdef __cplusplus
}
#endif
	  
typedef struct vector_of_3D_points_int 
 {
  int x;
  int y;
  int z;
 }Points_3dInt; 
 



  typedef struct pts{
    unsigned int x;
    unsigned int y;
    unsigned int z;
  } pt;


void Segmentation   (float ***ImageSeuilHaut     , //Image representant le seuil haut
                     float ***ImageSeuilBas      , //Image representant le seuil bas
                     grphic3d *ImageASegmenter    , //Image a segmenter
                     grphic3d *ImageSegmentee      //Image resultat
		     ); 


double minimisation22(
		    grphic3d *ImageTemplateOsDil  , 
		    grphic3d *ImageASegmenter  ,
		    grphic3d *ImageResultat    ,
                    int nb_param	       ,
		    double *param	       ,
		    double seuilBas            ,
		    double seuilHaut           ,
		    int,grphic3d *ImageTemplateOs   );		       



			 
void TranslateBase3dToWholeImage(float ***ImageRes     ,
				 double *param          ,
			   	 int numImage           ,	 
			    	 int nb_param      	,
				 double seuil          );

void TranslateBase3dToImage(float ***ImageRes     ,
			    double *param          ,
			    int numImage           ,	 
			    int nb_param           ,			    
			    int topi               ,
			    int topj               ,
			    int topk               ,
			    double seuil          );



int imx_segmentation_p2(grphic3d *ImageTemplate , 
			grphic3d *ImageTemplate2 , 
		       grphic3d *ImageASegmenter,
		       grphic3d *ImageResultat  ,
		       double   *param	        ,
		       int      nb_param        ,
		       double   seuilHaut       ,
		       double   seuilBas	,                       
		       field3d *chres,
		       int );

int imx_matching_Bspline_topo_primitive_non_seg_3d(int im_ref,int im_ref2, int im_reca, int im_res, int func_type, int dist_type, int reg_type, int inter_type, int min_type, int save_type, char *nomfichres, int
resolf, double Jmin, double Jmax);

int imx_matching_Bspline_topo_primitive_non_seg_3d_p(grphic3d *imref,grphic3d *imref2,  grphic3d *imreca, grphic3d *imres,
                                int func_type, int dist_type, int reg_type, int inter_type,
                                int min_type, int save_type, char *nomfichres,
                                int resolf, double Jmin, double Jmax);


int imx_convert_fmatrix3d_grphic3d(float ***mat,grphic3d* im);
double gaussienneValue(double x,double mu,double sigma);
double **alloc_matrice_2d_doubleT(int width, int height);
void free_matrice_2d_doubleT(double **D);

int ChoixSeuil(std::vector<unsigned int> cumul, int min,double p1);
void filtrage(double *y, long int nbpts);
void EstimeGauss(double *y, long int nbpts);


//PROTO NICOLAGET_PS
void MarqueCCObjet2(int x,int y,int z,char tampon[3][3][3]);
int CCObjet2(char I[3][3][3]);
void MarqueCCObjet(int x,int y,int z,char tampon[3][3][3]);
int calcT1(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres);
int calcT2(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres);
bool simple(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres);
bool simple2(unsigned int x,unsigned int y,unsigned int z,grphic3d *imres);
void reduction_homotopique(grphic3d *imin,short***imth,short***imtb,grphic3d *imres,int x,int y,int z);
void generation_tunnels(grphic3d *imin,short***imth,short***imtb,grphic3d *imres);
void separation_composantes(grphic3d *imin,short***imth,short***imtb,grphic3d *imres);
void segmentationNicolas(grphic3d*imin,short***imth,short***imtb,grphic3d*imres);
void segmentationNicolas2(short***imth,short***imtb,grphic3d*imin,grphic3d*impr,grphic3d*imres,
            Points_3dInt*pointBG,Points_3dInt*pointHD,t_ptrlist*add,t_ptrlist*remove);
int free_matrix_3dSHORT(short int ***m);
short int ***alloc_matrix_3dSHORT(int nb_lig, int nb_col, int nb_dep);


//
void DetermineParameter(double *param	             ,
			int nb_param,
		  	grphic3d * ImageTemplateOsDilate,
			float *** ImageSeuilBas      ,
			double seuilBas               ,
		  	int ***masque_param	      );
int OptimisationLigne(int topi	     	         ,
			 int topj	     	         ,
			 int topk	     	         ,
			 double *param	            	 ,
		    	 double *moyenne		 ,
			 int nb_param		         , 
		    	 grphic3d * ImageTemplateOsDilate,
		    	 float *** ImageSeuilBas	 ,
		    	 double seuilBas                 ,		    
		    	 int ***masque_param	         ,
			 double pas                      ) ;
double UpdateFunctionCout(int topi,int topj,int topk,
		    double *param	            ,
		    double *moyenne,
		    int nb_param		    ,
		    grphic3d * ImageTemplateOsDilate,
		    float *** ImageSeuilBas	    ,
		    double seuilBas                 ,		    
		    int ***masque_param	            );
double FunctionCoutZZ(double *moyenne	              ,
		      double *param                   ,
		      int nb_param		      ,
		      grphic3d * ImageTemplateOsDilate,
		      float *** ImageSeuilBas	      ,
		      double seuilBas                 ,		    
		      int ***masque_param	      );
#endif
