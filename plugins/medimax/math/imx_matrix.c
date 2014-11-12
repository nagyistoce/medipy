/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/*!-----------------------------------------------------------------------
***	
***	\file:       imx_matrix.c
***
***	project:    Imagix 2.01 
***			
***
***	\brief description:	operation matrice .
***	
***	
***	Copyright (c) 1993, ULP-IPB Strasbourg.
***	All rights are reserved.
***
***     Last user action by: Mr. ARMSPACH on Aug 17th 1994
***     Last user action by: Mr. NIKOU    on Dec 13th 1995
***
***---------------------------------------------------------------------*/

#include <config.h>
#include "noyau/imagix.h"
#include "math/imx_matrix.h"
#include <time.h>


#define ITMAX 100
#define TINY 1.0e-20 
#define GOLD 1.618034
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GLIMIT 100.0
#define TOL 2.0e-4
#define EPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define DSQR(x)   ((x)*(x))
#define SIGNX(a,b) ((b)>=0.0 ? fabs(a) : -fabs(a))
#define DMAX(a,b) ((a)>=(b) ? a : b)
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	                  a[k][l]=h+s*(g-h*tau);

int ncom;
double *pcom,*xicom,(*nrfunc)();


#if defined (WIN32) || defined (CYGWIN) // CORRECTION RC. Ajout CYGWIN.
/********************************************************
**
** fonction type xrand48 non disponible sous NT
** source trouvees sur le site :
** http://www.freebsd.org/cgi/cvsweb.cgi/src/lib/libc/gen
**
*********************************************************/
#define	RAND48_SEED_0	(0x330e)
#define	RAND48_SEED_1	(0xabcd)
#define	RAND48_SEED_2	(0x1234)
#define	RAND48_MULT_0	(0xe66d)
#define	RAND48_MULT_1	(0xdeec)
#define	RAND48_MULT_2	(0x0005)
#define	RAND48_ADD		(0x000b)

unsigned short _rand48_seed[3];
unsigned short _rand48_mult[3];
unsigned short _rand48_add = RAND48_ADD;
#endif

/************************************************************************/

/******************************************************************
***    Sous programme de manipulation et calcul
**     matriciel
***/
/*********************************************************************/

/*!********************************************************************
**
**	\brief allocation d'un vecteur d'int
**	\param n : nbre d'element
**	\retval int *  : vecteur malloc'e
**
**********************************************************************/
int *alloc_ivector(int n) /* n=nb des elements du vecteur */

 {
 int *v;

 v=(int *)malloc((size_t)(n*sizeof(int)));

 return(v);
}

/*!********************************************************************
**
**	\brief allocation d'un vecteur de float
**	\param n : nbre d'element
**	\retval float *  : vecteur malloc'e
**
**********************************************************************/
 float *alloc_vector(int n) /* n=nb des elements du vecteur */

 {
  float *v;

  v=(float *)malloc((size_t)(n*sizeof(float)));

  return (v);
 }


/*!********************************************************************
**
**	\brief allocation d'un vecteur de double
**	\param n : nbre d'element
**	\retval double *  : vecteur malloc'e
**
**********************************************************************/
double *alloc_dvector(int n) /* n=nb des elements du vecteur */

 {
  double *v;

  v=(double *)malloc((size_t)(n*sizeof(double)));

  return (v);

 }

/*!********************************************************************
**
**	\brief rallocation d'un vecteur de double
**	\param vector : le vecteur
**	\param n : nbre d'element
**	\retval double *  : vecteur realloc'e
**
**********************************************************************/
double *realloc_dvector(double *vector, int n) /* n=nb des elements du vecteur */
 {

  vector=(double *)realloc((double *)vector,(size_t)((n+1)*sizeof(double)));

  return (vector);

 }

/*********************************************************************/
int *realloc_ivector(int *vector, int n) /* n=nb des elements du vecteur */


 {

  vector=(int *)realloc((int *)vector,(size_t)((n+1)*sizeof(int)));

  return (vector);

 }

/*********************************************************************/
float **alloc_matrix(int nb_lig, int nb_col)
{
  float **m;
  int i;

  m=(float **) malloc((size_t)((nb_lig+1)*sizeof(float *)));

  if (!m) { printf("\n allocation failed 1"); exit(1); }

  m[0]=(float *)malloc((size_t)((nb_lig*nb_col)*sizeof(float)));

  if (!m[0])  { printf("\n allocation failed 2 m[0]"); exit(1); }

  for(i=1;i<nb_lig;i++)
     m[i]=m[i-1]+nb_col;

  return m;
 }

/*********************************************************************/
int **alloc_imatrix(int nb_lig, int nb_col)
{
  int **m;
  int i;

  m=(int **) malloc((size_t)((nb_lig)*sizeof(int *)));

  if (!m){ printf("\n allocation failed 1"); exit(1);}

  m[0]=(int *)malloc((size_t)((nb_lig*nb_col)*sizeof(int)));

  if (!m[0]) {printf("\n allocation failed 2 m[0]="); exit(1);}

  for(i=1;i<nb_lig;i++)
     m[i]=m[i-1]+nb_col;


  return (m);
 }

/*!********************************************************************
**
**	\brief allocation d'une matrice de double
**	\param nb_lig, nb_col : dimension de la matrice
**	\retval la matrice
**  \remark En cas d'echec exit(1)
**
**********************************************************************/
double **alloc_dmatrix(int nb_lig, int nb_col)
{
  double **m;
  int i;

  m=(double **) malloc((size_t)((nb_lig+1)*sizeof(double *)));

  if (!m){ printf("\n allocation failed 1"); exit(1);}

  m[0]=(double *)malloc((size_t)((nb_lig*nb_col)*sizeof(double)));

  if (!m[0]) {printf("\n allocation failed 2 m[0]=0"); exit(1);}

  for(i=1;i<nb_lig;i++)
     m[i]=m[i-1]+nb_col;

  return (m);
 }

/*********************************************************************/

void free_matrix(float **m, int co, int li)
{
  free((char *)(m[0]));
  free((char*)(m));
 }

/*!********************************************************************
**
**	\brief desallocation d'un vecteur de float
**	\param v : le vecteur a desallouer (free)
**	\param n : nbre d'element
**
**********************************************************************/
void free_vector(float *v, int n)
{

  free((char *)(v));

 }

/*!********************************************************************
**
**	\brief desallocation d'un vecteur d'int
**	\param v : le vecteur a desallouer (free)
**	\param n : nbre d'element
**
**********************************************************************/
void free_ivector(int *v, int n)
{

  free(v);

 }

/*********************************************************************/

void free_imatrix(int **m, int co, int li)
{
  free((char *) (m[0]));
  free((char*)(m));
 }


/*!********************************************************************
**
**	\brief desallocation d'une matrice de double
**	\param m : la matrice a desallouer (free)
**
**********************************************************************/
void free_dmatrix(double **m, int co, int li)
{
  free((char *) (m[0]));
  free((char*)(m));
 }

/*!********************************************************************
**
**	\brief desallocation d'un vecteur de double
**	\param v : le vecteur a desallouer (free)
**	\param n : nbre d'element
**
**********************************************************************/

void free_dvector(double *v, int n)
{
  free((char *)(v));
 }

/*!********************************************************************
**
**	\brief multiplication de 2 matrices de doubles
**	\param a : matrice
**	\param l1 : nbre de ligne de a
**	\param C1 : nbre de colonnes de a
**	\param b : matrice
**	\param l1 : nbre de ligne de b
**	\param C1 : nbre de colonnes de b
**	\retval Une matrice nouvellement allouee (alloc_dmatrix)
**********************************************************************/
 double **multiply_matrix(double **a, int l1, int c1, double **b, int l2, int c2)
{
  int i,j,k;
  double **c,sum;

  if(c1!=l2) {printf("\n multiplication impossible\n"); exit(1);}

  c=alloc_dmatrix(l1,c2);

 if(c==NULL)
    printf("allocation impossible !\n");

  for(i=0;i<l1;i++)
     for(j=0;j<c2;j++)
        {
         sum=0.0;
           for(k=0;k<c1;k++)
              sum+=a[i][k]*b[k][j];

         c[i][j]=sum;
        }

  return(c);

 }

/*********************************************************************/

 double **add_matrix(double **a, int l, int c, double **b)
{
  double **d;
  int i,j;

  d=alloc_dmatrix(l,c);

  for(i=0;i<l;i++)
     for(j=0;j<c;j++)
        d[i][j]=a[i][j]+b[i][j];

  return(d);
 }

/*********************************************************************/

 double **sub_matrix(double **a, int l, int c, double **b)
{
  double **d;
  int i,j;

  d=alloc_dmatrix(l,c);

  for(i=0;i<l;i++)
     for(j=0;j<c;j++)
        d[i][j]=a[i][j]-b[i][j];

  return(d);
 }


/*********************************************************************/
/*! 
**	muliplie une matrice par un vecteur
**
**	\param aa : matrice  a[l][c]
**	\param l : nbre de lignes
**  \param c : nbre de colonnes
**  \param b : vecteur b[c]
**  \retval le vecteur resultat
*****************************************************************/
 double *matrix_multiply_vector(double **a, int l, int c, double *b)
{
  int i,j;
  double *d;

  d=alloc_dvector(l);

     for(j=0;j<l;j++)
        d[j]=0.0;

  for(i=0;i<l;i++)
     for(j=0;j<c;j++)
        d[i]+=a[i][j]*b[j];

  return(d);

 }

/*********************************************************************/
/*! 
**	\brief calcul d'une transposee de matrice
**	\param a : matrice
**	\param li : nbre de lignes
**	\param co : nbre de colonnes
**	\retval Une matrice nouvellement allouee ( alloc_dmatrix )
**********************************************************************/
 double **matrix_transp(double **a, int li, int co)
{
  double **c;
  int i,j;

  c=alloc_dmatrix(co,li);

  for(i=0;i<li;i++)
     for(j=0;j<co;j++)
        c[j][i]=a[i][j];

  return(c);

 }

/*********************************************************************/

 double *substract_vectors(double *a, double *b, int n)
{
  double *c;
  int i;

  c=alloc_dvector(n);

  for(i=0;i<n;i++)
     c[i]=a[i]-b[i];

  return(c);

 }


/*********************************************************************/

 double *add_vectors(double *a, double *b, int n)
{
  double *c;
  int i;

  c=alloc_dvector(n);

  for(i=0;i<n;i++)
     c[i]=a[i]+b[i];

  return(c);

 }

/*********************************************************************/

 double **matrix_multiply_matrix_transp(double **a, int li1, int co1, double **b, int li2, int co2)
{
  double **c,**t;


  t=matrix_transp(b,li2,co2);

  c=multiply_matrix(a,li1,co1,t,co2,li2);

  return(c);

}

/*********************************************************************/

 double **cons_multiply_matrix(double c, double **a, int li, int co)
{
  int i,j;
  double **d;

  d=alloc_dmatrix(li,co);

  for(i=0;i<li;i++)
     for(j=0;j<co;j++)
        d[i][j]=c*a[i][j];

  return(d);

 }

/*********************************************************************/

 double *cons_multiply_vector(double c, double *v, int n)
{
  int i;
  double *d;

  d=alloc_dvector(n);

  for(i=0;i<n;i++)
     d[i]=c*v[i];

  return(d);

 }


/*********************************************************************/

 double **vector_multiply_vector_transp(double *v, int n)
{
  int i,j;
  double **x;

  x=alloc_dmatrix(n,n);

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        x[i][j]=v[i]*v[j];

  return (x);

 }

/*********************************************************************/

void lu_decomposition(double **a, int n, int *index, double *d)
{
  int i=0,imax=0,j=0,k=0;
  double big,dum,sum,temp;
  double *v;
  double tiny=1.0e-20;
  double **b;

  b=alloc_dmatrix(n+1,n+1);
  v=alloc_dvector(n+1);

  for(i=1;i<n+1;i++)
     for(j=1;j<n+1;j++)
        b[i][j]=a[i-1][j-1];

  *d=1.0;

  for(i=1;i<=n;i++)
     {
      big=0.0;

      for(j=1;j<=n;j++)
         if((temp=fabs(b[i][j]))>big)
           big=temp;

       if(big==0.0)   // Modif JPA 2/2003 big=0.0
         printf("la matrice a une ligne de zeros \n");

       v[i]=1.0/big;
      }

  for(j=1;j<=n;j++)
     {

       for(i=1;i<j;i++)
          {
           sum=b[i][j];

	   for(k=1;k<i;k++)
              sum-=b[i][k]*b[k][j];

	   b[i][j]=sum;
          }

        big=0.0;

	for(i=j;i<=n;i++)
           {
            sum=b[i][j];

	    for(k=1;k<j;k++)
               sum-=b[i][k]*b[k][j];

	    b[i][j]=sum;

	    if((dum=v[i]*fabs(sum))>big)
              {
               big=dum;
               imax=i;
              }
            }

          if(j!=imax)
            {
             for(k=1;k<=n;k++)
                {
                 dum=b[imax][k];
                 b[imax][k]=b[j][k];
                 b[j][k]=dum;
                }
             *d=-(*d);
             v[imax]=v[j];
            }

           index[j]=imax;

	   if(b[j][j]==0.0)
             b[j][j]=tiny;

	    if(j!=n)
              {
               dum=1.0/(b[j][j]);

	       for(i=j+1;i<=n;i++)
                  b[i][j]*=dum;
              }

         }

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        a[i][j]=b[i+1][j+1];


  free_dmatrix(b,n+1,n+1);
  free_dvector(v,n+1);


 }

/*****************************************************************/

void back_subst(double **a, int n, int *index, double *b) /*************** Attenton b[1...n] *******************/
 			/*************** Attenton a[0...n-1][0...n-1] *******************/


 {
  int i,ii=0,ip,j;
  double sum;
  double **c,*x;

  c=alloc_dmatrix(n+1,n+1);
  x=alloc_dvector(n+1);

  for(i=1;i<n+1;i++)
     for(j=1;j<n+1;j++)
        c[i][j]=a[i-1][j-1];

  for(i=1;i<n+1;i++)
     x[i]=b[i];

  for(i=1;i<=n;i++)
     {
      ip=index[i];
      sum=x[ip];
      x[ip]=x[i];

      if(ii)
        for(j=ii;j<=i-1;j++)
           sum-=c[i][j]*x[j];

      else if(sum)
             ii=i;

      x[i]=sum;

     }

  for(i=n;i>=1;i--)
     {
      sum=x[i];

      for(j=i+1;j<=n;j++)
         sum-=c[i][j]*x[j];

      x[i]=sum/c[i][i];
     }

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        a[i][j]=c[i+1][j+1];


  for(i=1;i<=n;i++)
     b[i]=x[i];

 free_dmatrix(c,n+1,n+1);
 free_dvector(x,n+1);

 }

/*****************************************************************/

void linear_system(double **a, int n, double *b)
{
  int *index;
  double d;

  index=alloc_ivector(n+1);

  lu_decomposition(a,n,index,&d);
  back_subst(a,n,index,b);

  free_ivector(index,n+1);

 }

/*****************************************************************/

 double **matrix_inversion(double **a, int n) /* ATTENTION the input matrix a is destroyed */
 {
  double **y,d,*col,**c;
  int i,j,*index;
  int flag;

  y=alloc_dmatrix(n+1,n+1);
  c=alloc_dmatrix(n,n);
  col=alloc_dvector(n+1);
  index=alloc_ivector(n+1);


  for(i=0;i<n;i++)
     {
      flag=0;

      for(j=0;j<n;j++)
         if(a[i][j]==0.0)

      flag++;

      if(flag==n)
        a[i][0]=1.0e-10;
     }



  lu_decomposition(a,n,index,&d);

  for(j=1;j<=n;j++)
     {
      for(i=1;i<=n;i++)
         col[i]=0.0;

      col[j]=1.0;

      back_subst(a,n,index,col);

      for(i=1;i<=n;i++)
         y[i][j]=col[i];
     }

   for(i=0;i<n;i++)
      for(j=0;j<n;j++)
         c[i][j]=y[i+1][j+1];

  free_dmatrix(y,n+1,n+1);
  free_dvector(col,n+1);
  free_ivector(index,n+1);

  return(c);

 }

/****************************************************************
**	 matrix_inversion_sym(a,n)
** inversion d'une matrice symetrique definie positive par
** la methode de Cholesky
** Si la matrice a n'est pas symetrique, rine n'est signale,
** la procedure ne travaillant que sur la partie superieure de la
** matrice
** Si la matrice n'est pas definie positive, un message d'erreur
** est affiche
** ATTENTION: la partie inferieur de la matrice a est detruite
*****************************************************************/
 double **matrix_inversion_sym(double **a, int n)
{
 double *p,sum,**c;
 int i,j,k;

 c=alloc_dmatrix(n,n);
 for (i=0;i<n;i++) for (j=0;j<n;j++) if (i==j) c[i][j]=0.0; else c[i][j]=1.0;

 p=CALLOC(n,double); /*vecteur contenant les valeurs propre*/

 /*decomposition de Choleski: L.LT=A
 L mis dans la partie inferieure de A*/
 for (i=0;i<n;i++)
  for (j=i;j<n;j++)
  {
   sum=a[i][j];
   for (k=i-1;k>=0;k--) sum-=a[i][k]*a[j][k];
   if (i==j)
    {
     if (sum<=0.0)
      {printf("La matrice n'est pas definie positive RESULTAT= MATRICE IDENTITE\n");return(c);}
     p[i]=sqrt(sum);
    }
   else
    a[j][i]=sum/p[i];
  }

 /*calcul de l'inverse de la matrice triangulaire L*/
 for (i=0;i<n;i++)
  {
   a[i][i]=1.0/p[i];
   for (j=i+1;j<n;j++)
    {
     sum=0.0;
     for (k=i;k<j;k++) sum-=a[j][k]*a[k][i];
     a[j][i]=sum/p[j];
    }
  }

  /*resultat de l'inversion dans c: la matrice inverse de a est le produit
  de L-1 par LT-1*/

  for (i=0;i<n;i++)
   for (j=0;j<=i;j++)
   {
    sum=0.0;
    for (k=i;k<n;k++) sum+=a[k][i]*a[k][j];
    c[i][j]=sum;
   }
  /*la matrice c est symetrique*/
  for (i=0;i<n;i++)
   for (j=0;j<i;j++)
    c[j][i]=c[i][j];

 FREE(p);

 return(c);
}

/*****************************************************************/

void matrix_tridiagonal_reduction(double **aa, int n, double *d, double *e, double **x)
                         /* d contains the diagonal elements and e the off-diagonal */

 {
  int i,j,k,l;
  double scale,hh,h,f,g,**a;

  a=alloc_dmatrix(n+1,n+1);

  for(i=1;i<n+1;i++)
     for(j=1;j<n+1;j++)
        a[i][j]=aa[i-1][j-1];

  for(i=n;i>=2;i--)
     {
      l=i-1;
      h=scale=0.0;

	 if(l>1)
           {
            for(k=1;k<=l;k++)
               scale+=fabs(a[i][k]);

	    if(scale==0.0)
              e[i]=a[i][l];

	    else
	        {
		 for(k=1;k<=l;k++)
		    {
		     a[i][k]/=scale;;
		     h+=a[i][k]*a[i][k];
		     }

		   f=a[i][l];
		   g=(f>0 ? -sqrt(h) : sqrt(h));
		   e[i]=scale*g;
		   h-=f*g;
		   a[i][l]=f-g;
		   f=0.0;

		   for(j=1;j<=l;j++)
		      {
		       a[j][i]=a[i][j]/h;
		       g=0.0;

		       for(k=1;k<=j;k++)
			  g+=a[j][k]*a[i][k];

		       for(k=j+1;k<=l;k++)
			  g+=a[k][j]*a[i][k];

		       e[j]=g/h;
		       f+=e[j]*a[i][j];
		      }

		    hh=f/(h+h);

		    for(j=1;j<=l;j++)
		      {
			f=a[i][j];
			e[j]=g=e[j]-hh*f;

			for(k=1;k<=j;k++)
			  a[j][k]-=(f*e[k]+g*a[i][k]);
		      }
		  }/* end else */
	     }

      else
          e[i]=a[i][l];

      d[i]=h;
   }

  d[1]=0.0;
  e[1]=0.0;

  for(i=1;i<=n;i++)
     {
      l=i-1;
      if(d[i])
        {
	 for(j=1;j<=l;j++)
	    {
	     g=0.0;

	     for(k=1;k<=l;k++)
	        g+=a[i][k]*a[k][j];

   	     for(k=1;k<=l;k++)
	        a[k][j]-=g*a[k][i];
	    }
         }

      d[i]=a[i][i];
      a[i][i]=1.0;

      for(j=1;j<=l;j++)
         a[j][i]=a[i][j]=0.0;
     }

  for(i=1;i<=n;i++)
     for(j=1;j<=n;j++)
        x[i][j]=a[i][j];

 free_dmatrix(a,n+1,n+1);
 }

/***************************************************************/

 double pythag(double a, double b)
{
  double absa,absb,a2,b2;

  absa=fabs(a);
  absb=fabs(b);
  a2=absa*absa;
  b2=absb*absb;

  if(absa>absb)
    return(absa*sqrt(1.0+(b2/a2)));

  else
     return(absb==0 ? 0.0 : (absb*sqrt(1.0+(a2/b2))));

 }

/*************************************************************************/

 double **eigen_ql(double *d, double *e, int n, double **x, int sort)
{
  int m,l,iter,i,k,j;
  double s,r,p,g,f,dd,c,b,pp;
  double **z,**v;

  z=alloc_dmatrix(n+1,n+1);
  v=alloc_dmatrix(n+1,n);



  for(i=1;i<=n;i++)
     for(j=1;j<=n;j++)
        z[i][j]=x[i][j];

  for(i=2;i<=n;i++)
     e[i-1]=e[i];

  e[n]=0.0;

  for(l=1;l<=n;l++)
     {
      iter=0;
		do
		  {
		   for(m=1;m<=n-1;m++)
			{
			  dd=fabs(d[m])+fabs(d[m+1]);
			  if((double)(fabs(e[m])+dd)==dd) break;
			}
		   if(m!=1)
			 {
			  if(iter++==30)
			  printf("Too many iterations (%d) in eigen_ql algorithm\n",iter);
                          g=(d[l+1]-d[l])/(2.0*e[l]);
			  r=pythag(g,1.0);
			  g=d[m]-d[l]+e[l]/(g+SIGNX(r,g));
			  s=c=1.0;
			  p=0.0;

				for(i=m-1;i>=l;i--)
				   {
				     f=s*e[i];
				     b=c*e[i];
				     e[i+1]=(r=pythag(f,g));
                                           if(r==0.0)
					     {
					       d[i+1]-=p;
					       e[m]=0.0;
					       break;
					     }


			s=f/r;
			c=g/r;
			g=d[i+1]-p;
			r=(d[i]-g)*s+2.0*c*b;
			d[i+1]=g+(p=s*r);
			g=c*r-b;

				for(k=1;k<=n;k++)
				   {
				    f=z[k][i+1];
				    z[k][i+1]=s*z[k][i]+c*f;
				    z[k][i]=c*z[k][i]-s*f;
				   }
				}


		if(r==0.0 && i) continue;
		d[l]-=p;
	        e[l]=g;
		e[m]=0.0;
	      }

	}

	while(m!=1);

    }



/* sorting eigenvalues and eigenvectors */


   if(sort!=0)
     {
      // printf("sorting eigenvalues\n");

      for(i=1;i<n;i++)
         {
          pp=d[k=i];

          for(j=i+1;j<=n;j++)
             if(d[j]>=pp)
		pp=d[k=j];

          if(k!=i)
            {
	     d[k]=d[i];
             d[i]=pp;

	       for(j=1;j<=n;j++)
                  {
 	           pp=z[j][i];
	           z[j][i]=z[j][k];
	           z[j][k]=pp;
                  }
             }
          }

      } /* end sorting */

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        v[i][j]=z[i+1][j+1];


   for(j=0;j<n;j++)
      v[n][j]=d[j+1];


  free_dmatrix(z,n+1,n+1);

  return(v);


 }




/***************************************************************************/

 double **eigen_values(double **a, int n, int sort)
{
  double **v;
  double **x;
  double *d;
  double *e;

  d=alloc_dvector(n+1);
  e=alloc_dvector(n+1);
  x=alloc_dmatrix(n+1,n+1);

  matrix_tridiagonal_reduction(a,n,d,e,x);

  v=eigen_ql(d,e,n,x,sort);

  free_dvector(d,n+1);
  free_dvector(e,n+1);
  free_dmatrix(x,n+1,n+1);

  return (v);


 }

/***************************************************************************/

 double **eigen_values_jacobi(double **a, int n, int sort)
{
  int  i,j,p,q,k;
  double  seuil,theta,tau,t,sm,s,h,g,c,*b,*z;
  double **v, *d;
  double **eigen;
  double pp;

  b=alloc_dvector(n);
  z=alloc_dvector(n);
  d=alloc_dvector(n);
  v=alloc_dmatrix(n,n);
  eigen=alloc_dmatrix(n+1,n);


  for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) v[i][j] = 0.0;
      v[i][i] = 1.0;
     }
  for (j = 0; j < n ; j++) {
      b[j] = d[j] = a[j][j];
      z[j] = 0.0;
     }
  for (i = 0; i < 50; i++) {
      sm = 0.0;
      for (p = 0; p < n-1; p++) {
          for (q = p+1; q < n; q++) sm += fabs(a[p][q]);
         }


      if ( sm == 0.0) {

                         if(sort!=0)
                           {
                            // printf("sorting eigenvalues\n");

                            for(i=0;i<n;i++)
                               {
                                pp=d[k=i];

                                for(j=i+1;j<n;j++)
                                   if(d[j]>=pp)
		                     pp=d[k=j];

                                if(k!=i)
                                  {
	                           d[k]=d[i];
                                   d[i]=pp;

	                           for(j=0;j<n;j++)
                                      {
 	                               pp=v[j][i];
	                               v[j][i]=v[j][k];
	                               v[j][k]=pp;
                                      }
                                  }
                               }

                            } /* end sorting */

		       for(i=0;i<n;i++)
                          for(j=0;j<n;j++)
	                     eigen[i][j]=v[i][j];

                       for(i=0;i<n;i++)
                          eigen[n][i]=d[i];


                       free_dvector(b,n);
                       free_dvector(z,n);
                       free_dvector(d,n);
                       free_dmatrix(v,n,n);

                       return(eigen);
                      }

     if (i < 3) seuil = 0.2*sm/(n*n);
     else seuil = 0.0;
     for (p = 0; p < n-1; p++) {
         for (q = p+1; q < n; q++) {
            g = 100.0*fabs(a[p][q]);
            if (i > 3 && (float)(fabs(d[p])+g) == (float)fabs(d[p])
                      && (float)(fabs(d[q])+g) == (float)fabs(d[q]))
               a[p][q] = 0.0;
            else if (fabs(a[p][q]) > seuil) {
               h = d[q] - d[p];
               if ((float)(fabs(h)+g) == (float)fabs(h))
                   t = (a[p][q])/h;
               else {
                   theta = 0.5*h/(a[p][q]);
                   t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                   if (theta < 0.0) t = -t;
                  }
               c = 1.0/sqrt(1+t*t);
               s = t*c;
               tau = s/(1.0+c);
               h = t*a[p][q];
               z[p] -= h; z[q] += h;
               d[p] -= h; d[q] += h;
               a[p][q] = 0.0;
               for (j = 0; j < p; j++) {
                   ROTATE(a,j,p,j,q)
                   }
               for (j = p+1; j < q; j++) {
                   ROTATE(a,p,j,j,q)
                   }
               for (j = q+1; j < n; j++) {
                   ROTATE(a,p,j,q,j)
                   }
               for (j = 0; j < n; j++) {
                   ROTATE(v,j,p,j,q)
                   }
              }
          }
      }
      for (p = 0; p < n; p++) {
          b[p] += z[p];
          d[p] = b[p];
          z[p] = 0.0;
         }
    }


/* sorting eigenvalues and eigenvectors */


   if(sort!=0)
     {
      // printf("sorting eigenvalues\n");

      for(i=0;i<n;i++)
         {
          pp=d[k=i];

          for(j=i+1;j<n;j++)
             if(d[j]>=pp)
		pp=d[k=j];

          if(k!=i)
            {
	     d[k]=d[i];
             d[i]=pp;

	       for(j=0;j<n;j++)
                  {
 	           pp=v[j][i];
	           v[j][i]=v[j][k];
	           v[j][k]=pp;
                  }
             }
          }

      } /* end sorting */


   for(i=0;i<n;i++)
      for(j=0;j<n;j++)
	 eigen[i][j]=v[i][j];

   for(i=0;i<n;i++)
      eigen[n][i]=d[i];


  free_dvector(b,n);
  free_dvector(z,n);
  free_dvector(d,n);
  free_dmatrix(v,n,n);

  return(eigen);

  printf("Too many iterations in routine eigen_values_jacobi()");
 }

/****************************************************************************/

 double *vecteur_moyen(double **vectors, int t, int n)

          /* t=nombre des vecteurs, n=dimension des vecteurs */
 {
  int i,k;
  double *vm;

  vm=alloc_dvector(n);

  for(i=0;i<n;i++)
     vm[i]=0.0;

  for(i=0;i<n;i++)
     for(k=0;k<t;k++)
	vm[i]+=(vectors[k][i]/((double)t));

  return(vm);

 }

/****************************************************************************/

 double **covariance_matrix(double **vectors, double *vm, int t, int n)

          /* t=nombre des vecteurs, n=dimension des vecteurs */
 {
  int i,j,k;
  double **cov;
  double vki,vkj;

  cov=alloc_dmatrix(n,n);

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
	cov[i][j]=0.0;

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        for(k=0;k<t;k++)
           {
            vki=(vectors[k][i]-vm[i]);
            vkj=(vectors[k][j]-vm[j]);
            cov[i][j]+=(vki*vkj);
	   }

  for(i=0;i<n;i++)
     for(j=0;j<n;j++)
        cov[i][j]/=((double)t);

  return(cov);

 }

/****************************************************************************/

 double **covariance_of_small_training_set(double **vectors, double *vm, int t, int n)

          /* t=nombre des vecteurs, n=dimension des vecteurs */
 {
  int i,j,k;
  double **cov;
  double vki,vkj;


  cov=alloc_dmatrix(t,t);

  for(i=0;i<t;i++)
     for(j=0;j<t;j++)
	cov[i][j]=0.0;

  for(i=0;i<t;i++)
     for(j=0;j<t;j++)
        for(k=0;k<n;k++)
           {
            vki=(vectors[i][k]-vm[k]);
            vkj=(vectors[j][k]-vm[k]);
            cov[i][j]+=(vki*vkj);
	   }

  for(i=0;i<t;i++)
     for(j=0;j<t;j++)
        cov[i][j]/=((double)t);

  return(cov);

 }

/*********************************************************************/

 double brent(double ax, double bx, double cx, double (*f) (/* ??? */), double tol, double *xmin)
{
	int iter;
	double a=0,b=0,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGNX(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGNX(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	printf("Too many iterations in brent\n");
	*xmin=x;
	return fx;
 }

/***************************************************************************************************/

 void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func) (/* ??? */))
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGNX(DMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

/***************************************************************************************************/

 double f1dim(double x)
{
  int j;
  double f,*xt;

  xt=alloc_dvector(ncom);

  for(j=0;j<ncom;j++)
     xt[j]=pcom[j]+x*xicom[j];

  f=(*nrfunc)(xt);

  free_dvector(xt,ncom);

  return(f);

 }

/***************************************************************************************************/

 void linmin(double *p, double *xi, int n, double *fret, double (*func) (/* ??? */))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=alloc_dvector(n);
	xicom=alloc_dvector(n);
	nrfunc=func;
	for (j=0;j<n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=0;j<n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
 }

/***************************************************************************************************/

 void powell(double *p, double **xi, int n, double ftol, int *iter, double *fret, double (*func) (/* ??? */))
{
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;

	pt=alloc_dvector(n);
	ptt=alloc_dvector(n);
	xit=alloc_dvector(n);
	*fret=(*func)(p);
	for (j=0;j<n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			return;
		}
		if (*iter == ITMAX) printf("powell exceeding maximum iterations.\n");
		for (j=0;j<n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*DSQR(fp-(*fret)-del)-del*DSQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=0;j<n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
/***************************************************************************************************/

void powell_registration(double *p, int n, double (*func) (/* ??? */), double *fret, int *iter)
{
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;
	double ftol;
        double **xi;

	xi=alloc_dmatrix(n,n);

	for(i=0;i<n;i++)
	   for(j=0;j<n;j++)
              xi[i][j]=0.0;

	      for(i=0;i<n;i++)
		 xi[i][i]=1.0;

        ftol=TOL;

	pt=alloc_dvector(n);
	ptt=alloc_dvector(n);
	xit=alloc_dvector(n);
	*fret=(*func)(p);
	for (j=0;j<n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			return;
		}
		if (*iter == ITMAX) printf("powell exceeding maximum iterations.\n");
		for (j=0;j<n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*DSQR(fp-(*fret)-del)-del*DSQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=0;j<n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}

/***************************************************************************************************/

 void frprmn(double *p, int n, double ftol, int *iter, double *fret, double (*func) (/* ??? */), double (*dfunc) (/* ??? */))
{
	int j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

	g=alloc_dvector(n);
	h=alloc_dvector(n);
	xi=alloc_dvector(n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=0;j<n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=0;j<n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			return;
		}
		gam=dgg/gg;
		for (j=0;j<n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	printf("Too many iterations in frprmn\n");
 }

/*******************************************************************************************/

/******************************************************************************/
/*     affiche_dmatrix(m,nbl,nbc,format)                                      */
/*                                                                            */
/*     	affiche la matrix de double m de taille nbl*nbc avec un format        */
/*	specifie par la chaine format (idem printf)                           */
/******************************************************************************/
int	affiche_dmatrix(double **m, int nbl, int nbc, char *format)
{
 int i,j;

 for (i=0;i<nbl;i++)
 {
  for (j=0;j<nbc;j++) printf(format,m[i][j]);
  printf("\n");
 }


 return(0);
}

/******************************************************************************/
/*     affiche_imatrix(m,nbl,nbc,format)                                      */
/*                                                                            */
/*     	affiche la matrix de integer m de taille nbl*nbc avec un format        */
/*	specifie par la chaine format (idem printf)                           */
/******************************************************************************/
int	affiche_imatrix(int **m, int nbl, int nbc, char *format)
{
 int i,j;

 for (i=0;i<nbl;i++)
 {
  for (j=0;j<nbc;j++) printf(format,m[i][j]);
  printf("\n");
 }


 return(0);
}



//-----------------------------------------------------------------------------//
//                                                                             //
//                                 MATRICES 3D                                 //
//                                                                             //
//-----------------------------------------------------------------------------//

/*********************************************************************/
/*     alloc_dmatrix_3d(nb_lig,nb_col,nb_dep)                        */
/*                                                                   */
/*!     allocation memoire d'une matrice 3D
       le resultat est un triple pointeur sur un double
	   \param nb_lig, nb_col, nb_dep : dimension de la matrice
	   \retval la matrice alloue             						*/
/*********************************************************************/

double ***alloc_dmatrix_3d(int nb_lig, int nb_col, int nb_dep)
{
  double ***m;
  double *mt;
  int i,j;

  m=(double ***) malloc((size_t)((nb_lig+1)*sizeof(double **)));

  if (!m){ printf("\n allocation failed 1 in alloc_dmatrix_3d"); exit(1);}

  m[0]=(double **)malloc((size_t)((nb_lig*nb_col)*sizeof(double *)));

  if (!m[0]) {printf("\n allocation failed 2 in alloc_dmatrix_3d"); exit(1);}

  for(i=1;i<nb_lig;i++)
     m[i]=m[i-1]+nb_col;

  m[0][0]=(double *)malloc((size_t)((nb_lig*nb_col*nb_dep)*sizeof(double)));
  mt=m[0][0];

  if (!mt) {printf("\n allocation failed 3 in alloc_dmatrix_3d"); exit(1);}

  for (i=0;i<nb_lig;i++)
   for (j=0;j<nb_col;j++)
   {m[i][j]=mt;mt=mt+nb_dep;}


  return (m);
 }

/*********************************************************************/
/*     free_dmatrix_3d(m)                        		     		 */
/*                                                                   */
/*!     liberation memoire d'une matrice double 3D
		\param m : la matrice a liberer
		\retval 1												     */
/*********************************************************************/
int free_dmatrix_3d(double ***m)
{
  free((char *) (m[0][0]));
  free((char*)(m[0]));
  free((char*)m);

  return(1);
 }

/*********************************************************************/
/*     alloc_imatrix_3d(nb_lig,nb_col,nb_dep)                        */
/*                                                                   */
/*!     allocation memoire d'une matrice 3D
**      le resultat est un triple pointeur sur un integer
**		\param nb_lig : nbre de lignes
**		\param nb_col : nbre de col.
**		\param nb_dep : profondeur
**		\retval Une matrice allouee (a desallouer avec free_imatrix_3d)
*********************************************************************/
int ***alloc_imatrix_3d(int nb_lig, int nb_col, int nb_dep)
{
  int ***m;
  int *mt;
  int i,j;

  m=(int ***) malloc((size_t)((nb_lig+1)*sizeof(int **)));

  if (!m){ printf("\n allocation failed 1 in alloc_imatrix_3d"); exit(1);}

  m[0]=(int **)malloc((size_t)((nb_lig*nb_col)*sizeof(int *)));

  if (!m[0]) {printf("\n allocation failed 2 in alloc_imatrix_3d"); exit(1);}

  for(i=1;i<nb_lig;i++)
     m[i]=m[i-1]+nb_col;

  m[0][0]=(int *)malloc((size_t)((nb_lig*nb_col*nb_dep)*sizeof(int)));
  mt=m[0][0];

  if (!mt) {printf("\n allocation failed 3 in alloc_imatrix_3d"); exit(1);}

  for (i=0;i<nb_lig;i++)
   for (j=0;j<nb_col;j++)
   {m[i][j]=mt;mt=mt+nb_dep;}


  return (m);
 }


/*********************************************************************/
/*     free_imatrix_3d(m)                        		     */
/*                                                                   */
/*!     liberation memoire d'une matrice int 3D
**		\param m : matrice a liberer
**		\retval 1
*********************************************************************/
int free_imatrix_3d(int ***m)
{
  free((char *) (m[0][0]));
  free((char*)(m[0]));
  free((char*)m);

  return(1);
}

/*
  alloc_matrix_3d()

  allocation memoire d'une matrice 3D de floats

  nb_lig, nb_col, nb_dep : nombre de lignes, colones et plans de la matrice

  renvoie un pointeur sur un tableau de nb_lig pointeurs sur des tableaux
  de nb_col pointeurs sur des tableaux de nb_dep floats chaque (ouf !).

  Cet OVNI se manipule avec une syntaxe comparable a un
  "float toto[nb_lig][nb_col][nb_dep]". Les contraintes de contiguite memoire
  sont comparables (i.e. il vaut mieux iterer la derniere composante ("z") pour
  chaque valeur des premieres ("x" et "y") que l'inverse), a la difference pres
  qu'on doit effectuer en plus de multiples indirections par pointeurs.
*/
float ***alloc_matrix_3d(int nb_lig, int nb_col, int nb_dep)
{
  float ***m;
  float *mt;
  int i, j;

  m = malloc(sizeof(m[0])*nb_lig);

  if (! m)
  {
    return NULL;
  }

  m[0] = malloc(sizeof(m[0][0])*nb_lig*nb_col);

  if (! m[0])
  {
    free(m);
    return NULL;
  }

  for(i=1; i<nb_lig; i++)
     m[i] = m[i-1]+nb_col;

  mt = malloc(sizeof(m[0][0][0])*nb_lig*nb_col*nb_dep);

  if (! mt)
  {
    free(m[0]);
    free(m);
    return NULL;
  }

  for (i=0; i<nb_lig; i++)
    for (j=0; j<nb_col; j++)
    {
      m[i][j] = mt;
      mt += nb_dep;
    }

  return m;
}

/*
  free_matrix_3d()

  liberation memoire d'une matrice float 3D allouee par alloc_matrix_3d()
*/
int free_matrix_3d(float ***m)
{
  free(m[0][0]);
  free(m[0]);
  free(m);

  return(0);
}

/****************************************************************************/
void  sort_by_straight_insertion(double *vector, int n)
{
  int i,j;
  double a;

 for(j=1;j<n;j++)
    {
     a=vector[j];
     i=j-1;

     while(i>=0 && vector[i]>a)
          {
           vector[i+1]=vector[i];
           i--;
          }

      vector[i+1]=a;
    }

 }

/*************************************************************************************/
double  vector_median(double *vector, int n)
{
  double median;

  sort_by_straight_insertion(vector,n);

  if(n%2!=0)
    median=vector[n/2];

  else
     median=0.5*(vector[n/2]+vector[(n/2)-1]);

  return(median);

 }

/*************************************************************************************/
void arrange_vector(double *arr, int n, int *index)
{
  int i;
  double *wkspc;

  wkspc=alloc_dvector(n);

  for(i=0;i<n;i++)
      wkspc[i]=arr[i];

  for(i=0;i<n;i++)
      arr[i]=wkspc[index[i]];

   free_dvector(wkspc,n);

}

/********************************************************************
**     update_matrix_BFGS(Hk1,Hk2,paramk1,paramk2,gradk1,gradk2,nb_param)                 
*/                                                                   
/*!    mise a jour de la matrice Hbfgs 
**     \param Hk1: matrice a l'iteration k1
**     \param Hk2: matrice a l'iteration k2 (E/S)
**     \param paramk1: valeur des parametre a l'iteration k
**     \param gradk1: valeur du gradient a l'iteration k
**     \param paramk2: valeur des parametre a l'iteration k+1
**     \param gradk1: valeur du gradient a l'iteration k+1
**
**	   \retval 1
**********************************************************************/
int  update_matrix_BFGS(double **Hk1, double **Hk2, double *paramk1, double *paramk2, double *gradk1, double *gradk2, int nb_param)
{
 double *S,*Y,*AY,*U;
 double a,b;
 int i,j;
 
 /*allocation et remplissage des vecteur S,Y*/
 S=alloc_dvector(nb_param);
 /*S est le vecteur Xk+1-Xk*/
 for (i=0;i<nb_param;i++) S[i]=paramk2[i]-paramk1[i]; 
 Y=alloc_dvector(nb_param);
 /*Y est le vecteur Gk+1-Gk*/
 for (i=0;i<nb_param;i++) Y[i]=gradk2[i]-gradk1[i]; 
 
 /*calcul du scalaire a=S.Y*/
 a=0.0;
 for (i=0;i<nb_param;i++) a=a+S[i]*Y[i];

 /*calcul du vecteur Hk1Y*/ 
 AY=matrix_multiply_vector(Hk1,nb_param,nb_param,Y); 
 /*calcul du scalaire b=Y.Hk1.Y*/
 b=0.0;
 for (i=0;i<nb_param;i++) b=b+Y[i]*AY[i];
 /*allocation et calcul du vecteur U*/
 U=alloc_dvector(nb_param);
 for (i=0;i<nb_param;i++) U[i]=S[i]/a-AY[i]/b;
 
 /*calcul de la matrice Hk2*/
 if (a!=0) 
  {
   for (i=0;i<nb_param;i++)
    for (j=0;j<nb_param;j++)
     {
      Hk2[i][j]=Hk1[i][j]+S[i]*S[j]/a-AY[i]*AY[j]/b+b*U[i]*U[j]; 
     }
  }
 else
  {
//  aff_log("a=0\n");
   for (i=0;i<nb_param;i++)
    for (j=0;j<nb_param;j++)
     {
      Hk2[i][j]=Hk1[i][j]; 
     }
  } 

 
 /*liberation memoire des variables*/
 free_dvector(S,nb_param);free_dvector(Y,nb_param);
 free_dvector(AY,nb_param);free_dvector(U,nb_param);
 
 return(1);
}

int choldc(double **a, double *diag, int n)
{
 double sum;
 int i,j,k;

 for (i=0;i<n;i++)
  for (j=i;j<n;j++)
  {
   sum=a[i][j];
   for (k=i-1;k>=0;k--) sum-=a[i][k]*a[j][k];
   if (i==j)
   {
    if (sum<=0.0) { fprintf(stderr, "La matrice n'est pas definie positive RESULTAT= MATRICE IDENTITE\n"); return 1; }
    diag[i]=sqrt(sum);
   }
   else
    a[j][i]=sum/diag[i];
  }  

 return 0;
}

/*********************************************************************/
/*     Determinant(double **a,int n)                        		     */
/*                                                                   */
/*!     calcul recursif du determinant d'une matrice 2D
**		\param a : matrice dont le determinant est a calculer
**    \param n : taille de la matrice
**		\retval le determinant calcule
*********************************************************************/
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0.0;
   double **m = NULL;

   if (n < 1) { /* Error */
      fprintf (stderr, "erreur dans Determinant\n");
   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0.0;
      m = malloc((n-1)*sizeof(double *));
      for (i=0;i<n-1;i++) m[i] = malloc((n-1)*sizeof(double));
      for (j1=0;j1<n;j1++) {
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
      }
      for (i=0;i<n-1;i++) free(m[i]);
      free(m);
   }
   return(det);
}

/*******************************************
** 		  imx_convert_dmatrix3d_grphic3d() 
**
**  Conversion d'une matrice de double dans
**  une image grphic3d
**
** renvoie 0 si probl�e et 1 sinon
**  
********************************************/

int imx_convert_dmatrix3d_grphic3d(double ***mat,grphic3d* im) {
  int i,j,k;
  int wdth=im->width,hght=im->height,dpth=im->depth;
  double max,min,rcoeff;
  
  
  /* recherche du min et du max */
  min = HUGE_VAL; max	= -HUGE_VAL;
  
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++) {
	if (mat[i][j][k]<min)
          min=mat[i][j][k];
	
	if (mat[i][j][k]>max)
          max=mat[i][j][k];
      }
  
  
  imx_brukermax_3d(max, min, im);
  rcoeff=im->rcoeff;
  
  /* remplissage de l'image */
  for (i=0;i<wdth;i++)
    for (j=0;j<hght;j++)
      for (k=0;k<dpth;k++) {
	im->mri[i][j][k]=(int)(mat[i][j][k]/rcoeff);
      }
  
  return(1);
}


