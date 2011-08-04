/*************************************************************************
 * MediPy - Copyright (C) Universite de Strasbourg, 2011
 * Distributed under the terms of the CeCILL-B license, as published by
 * the CEA-CNRS-INRIA. Refer to the LICENSE file or to
 * http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
 * for details.
 ************************************************************************/

/********************************************************************************
**/	
/*!	   \file:        imx_math.c
***
***    project:    Imagix 2.01
***
***
***    \brief description:    fcts math. FFT, convolution, etc...
***
***    Fonctions decrites dans "Numerical recipes in C"
***
***---------------------------------------------------------------------*/
#include <config.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "noyau/imagix.h"
#include "math/imx_math.h"

#define NR_END 1

static int FFTRADIX (double *Re, double *Im, unsigned int nTotal, unsigned int nPass, unsigned int nSpan, int iSign, int max_factors, int max_perm);

/*
  Cette fonction initialise le generateur de nombres aleatoires a une valeur
  a peu pres aleatoire.

  Elle ne doit etre appelee qu'une seule fois, vers le debut du programme.
*/
void init_rand_seed(void)
{
//  struct tm when;
  time_t now, deadline;
  unsigned short the_seed[3];

  /* first set the random seed to the current date & time */
  now = time(NULL);
  srand((unsigned int) now);
#ifndef WIN32
  srandom((unsigned int) now);
#endif
  the_seed[0] = (unsigned short)now;
  the_seed[1] = (unsigned short)now >> (sizeof(unsigned short)); // correction berst 16;
  the_seed[2] = (unsigned short)now >> (sizeof(unsigned short)*2); // correction berst 32;
#if !defined( WIN32 ) && !defined( CYGWIN ) // CORRECTION RC
  seed48(the_seed);
#endif
  /* Usual UNIX implementations of time() only return a second count.
  In order to add further randommess, we need to take the fraction of second
  into account.  To do so, we run rand() repeatedly until the clock has
  incremented. */
//  localtime_r(& now,& when);
  /* wait from 1 through 2 seconds */
//  when.tm_sec += 2; 	/* result is in the range [2,61] */
  deadline = now+2;//mktime(&when);
  while (difftime(deadline, time(0)) > 0)
  {
    rand();
#if !defined( WIN32 ) && !defined( CYGWIN ) // CORRECTION RC
    random();
    drand48();
#endif
  }
}


void nrerror(char *msg)
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

float *vectorF(long int nl, long int nh)
{
	float *v;

	v=(float*)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
	if(!v)
		return NULL;

	return (float*)(v-nl+NR_END);
}

void free_vectorF(float *v, long int nl, long int nh)
{
	free((char*)(v+nl-NR_END));
}

void four1(float *data, long unsigned int nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta; /*Double precision for the trigonometric recurrences. */
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) /* This is the bit-reversal section of the routine. */
	{
		if (j > i)
		{
			SWAPF(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAPF(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m)
		{
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	/* Here begins the Danielson-Lanczos section of the routine. */
	mmax=2;
	while (n > mmax) /* Outer loop executed log 2 nn times. */
	{
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /*Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) /* Here are the two nested inner loops. */
		{
			for (i=m;i<=n;i+=istep)
			{
				j=i+mmax; /* This is the Danielson-Lanczos formula: */
				tempr=(float)(wr*data[j]-wi*data[j+1]);
				tempi=(float)(wr*data[j+1]+wi*data[j]);
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}



/*
Donnes deux tableaux de float en entree (data1[1..n] et data2[1..n]), cette fonction appelle
four1() et retourne deux tableaux (complexes, fft1[1..2n] et fft2[1..2n]), chacun de longueur
n (dimension reelle 2*n), et qui contiennent la TF discrete des tableaux de donnees respectifs.
n DOIT etre un entier puissance de 2
*/
void twofft(float *data1, float *data2, float *fft1, float *fft2, long unsigned int n)
{
	unsigned long nn3,nn2,jj,j;
	float rep,rem,aip,aim;
	nn3=1+(nn2=2+n+n);

	for (j=1,jj=2;j<=n;j++,jj+=2) /* Pack the two real arrays into one complex array */
	{
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1); /* Transform the complex array. */
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2)
	{
		rep=(float)(0.5*(fft1[j]+fft1[nn2-j])); /* Use symmetries to separate the two transforms */
		rem=(float)(0.5*(fft1[j]-fft1[nn2-j]));
		aip=(float)(0.5*(fft1[j+1]+fft1[nn3-j]));
		aim=(float)(0.5*(fft1[j+1]-fft1[nn3-j]));
		fft1[j]=rep; /* Ship them out in two complex arrays. */
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}

void convlv(float *data, long unsigned int n, float *respns, long unsigned int m, int isign, float *ans)
{
	unsigned long i,no2;
	float dum,mag2,*fft;

	fft=vectorF(1,n<<1);

	for (i=1;i<=(m-1)/2;i++) /* Put respns in array of length n. */
		respns[n+1-i]=respns[m+1-i];

	for (i=(m+3)/2;i<=n-(m-1)/2;i++) /* Pad with zeros. */
		respns[i]=0.0;

	twofft(data,respns,fft,ans,n); /* FFT both at once. */
	no2=n>>1;

	for (i=2;i<=n+2;i+=2)
	{
		if(isign == 1)
		{
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2; /* Multiply FFTs to convolve. */
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		}
		else if(isign == -1)
		{
			if ((mag2=(float)(SQR(ans[i-1])+SQR(ans[i]))) == 0.0)
				nrerror("Deconvolving at response zero in convlv");

			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2; /* Divide FFTs to deconvolve. */
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		}
		else
			nrerror("No meaning for isign in convlv");
	}

	ans[2]=ans[n+1]; /* Pack last element with rst for realft. */
	realft(ans,n,-1); /* Inverse transform back to time domain. */
	free_vectorF(fft,1,n<<1);
}


/*
Calculates the Fourier transform of a set of n real-valued data points. Replaces this data (which
is stored in array data[1..n]) by the positive frequency half of its complex Fourier transform.
The real-valued first and last components of the complex transform are returned as elements
data[1] and data[2], respectively. n must be a power of 2. This routine also calculates the
inverse transform of a complex data array if it is the transform of real data. (Result in this case
must be multiplied by 2/n.)

Calcule la TF d'un tableaux de n elements de type float. Remplace les donnees stockees dans data
par
*/

void realft(float *data, long unsigned int n, int isign)
{
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta; /* Double precision for the trigonometric recurrences. */

	theta=3.141592653589793/(double) (n>>1); /* Initialize the recurrence. */
	if (isign == 1)
	{
		c2 = -0.5;
		four1(data,n>>1,1); /* The forward transform is here. */
	}
	else
	{
		c2=0.5; /* Otherwise set up for an inverse transform. */
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;

	for(i=2;i<=(n>>2);i++)
	{
		/* Case i=1 done separately below. */
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]); /* The two separate transforms are separated out of data. */
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=(float)(h1r+wr*h2r-wi*h2i); /* Here they are recombined to form the true transform of the original real data. */
		data[i2]=(float)(h1i+wr*h2i+wi*h2r);
		data[i3]=(float)(h1r-wr*h2r+wi*h2i);
		data[i4]=(float)(-h1i+wr*h2i+wi*h2r);
		wr=(wtemp=wr)*wpr-wi*wpi+wr; /* The recurrence. */
		wi=wi*wpr+wtemp*wpi+wi;
	}

	if (isign == 1)
	{
		data[1] = (h1r=data[1])+data[2]; /* Squeeze the first and last data together to get them all within the original array. */
		data[2] = h1r-data[2];
	}
	else
	{
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1); /* This is the inverse transform for the case isign=-1. */
	}
}

/*--------------------------------*-C-*---------------------------------*
 * File:
 *	fftn.c
 *
 * Public:
 *	fft_free ();
 *	fftn / fftnf ();
 *
 * Private:
 *	fftradix / fftradixf ();
 *
 * Descript:
 *	multivariate complex Fourier transform, computed in place
 *	using mixed-radix Fast Fourier Transform algorithm.
 *
 *	Fortran code by:
 *	RC Singleton, Stanford Research Institute, Sept. 1968
 *
 *	translated by f2c (version 19950721).
 *
 * Revisions:
 * 26 July 95	John Beale
 *	- added maxf and maxp as parameters to fftradix()
 *
 *
 *	- added fft_free() to provide some measure of control over
 *	  allocation/deallocation.
 *
 *	- added fftn() wrapper for multidimensional FFTs
 *
 *	- use -DFFT_NOFLOAT or -DFFT_NODOUBLE to avoid compiling that
 *	  precision. Note suffix `f' on the function names indicates
 *	  float precision.
 *
 *	- revised documentation
 *
 * 31 July 95	Mark Olesen <olesen@me.queensu.ca>
 *	- added GNU Public License
 *	- more cleanup
 *	- define SUN_BROKEN_REALLOC to use malloc() instead of realloc()
 *	  on the first pass through, apparently needed for old libc
 *	- removed #error directive in favour of some code that simply
 *	  won't compile (generate an error that way)
 *
 * 1 Aug 95	Mark Olesen <olesen@me.queensu.ca>
 *	- define FFT_RADIX4 to only have radix 2, radix 4 transforms
 *	- made fftradix /fftradixf () static scope, just use fftn()
 *	  instead.  If you have good ideas about fixing the factors
 *	  in fftn() please do so.
 *
 * ======================================================================*
 * NIST Guide to Available Math Software.
 * Source for module FFT from package GO.
 * Retrieved from NETLIB on Wed Jul  5 11:50:07 1995.
 * ======================================================================*
 *
 *-----------------------------------------------------------------------*
 *
 * int fftn (int ndim, const int dims[], REAL Re[], REAL Im[],
 *	    int iSign, double scaling);
 *
 * NDIM = the total number dimensions
 * DIMS = a vector of array sizes
 *	if NDIM is zero then DIMS must be zero-terminated
 *
 * RE and IM hold the real and imaginary components of the data, and return
 * the resulting real and imaginary Fourier coefficients.  Multidimensional
 * data *must* be allocated contiguously.  There is no limit on the number
 * of dimensions.
 *
 * ISIGN = the sign of the complex exponential (ie, forward or inverse FFT)
 *	the magnitude of ISIGN (normally 1) is used to determine the
 *	correct indexing increment (see below).
 *
 * SCALING = normalizing constant by which the final result is *divided*
 *	if SCALING == -1, normalize by total dimension of the transform
 *	if SCALING <  -1, normalize by the square-root of the total dimension
 *
 * example:
 * tri-variate transform with Re[n1][n2][n3], Im[n1][n2][n3]
 *
 *	int dims[3] = {n1,n2,n3}
 *	fftn (3, dims, Re, Im, 1, scaling);
 *
 *-----------------------------------------------------------------------*
 * int fftradix (REAL Re[], REAL Im[], size_t nTotal, size_t nPass,
 *		 size_t nSpan, int iSign, size_t max_factors,
 *		 size_t max_perm);
 *
 * RE, IM - see above documentation
 *
 * Although there is no limit on the number of dimensions, fftradix() must
 * be called once for each dimension, but the calls may be in any order.
 *
 * NTOTAL = the total number of complex data values
 * NPASS  = the dimension of the current variable
 * NSPAN/NPASS = the spacing of consecutive data values while indexing the
 *	current variable
 * ISIGN - see above documentation
 *
 * example:
 * tri-variate transform with Re[n1][n2][n3], Im[n1][n2][n3]
 *
 *	fftradix (Re, Im, n1*n2*n3, n1,       n1, 1, maxf, maxp);
 *	fftradix (Re, Im, n1*n2*n3, n2,    n1*n2, 1, maxf, maxp);
 *	fftradix (Re, Im, n1*n2*n3, n3, n1*n2*n3, 1, maxf, maxp);
 *
 * single-variate transform,
 *    NTOTAL = N = NSPAN = (number of complex data values),
 *
 *	fftradix (Re, Im, n, n, n, 1, maxf, maxp);
 *
 * The data can also be stored in a single array with alternating real and
 * imaginary parts, the magnitude of ISIGN is changed to 2 to give correct
 * indexing increment, and data [0] and data [1] used to pass the initial
 * addresses for the sequences of real and imaginary values,
 *
 * example:
 *	REAL data [2*NTOTAL];
 *	fftradix ( &data[0], &data[1], NTOTAL, nPass, nSpan, 2, maxf, maxp);
 *
 * for temporary allocation:
 *
 * MAX_FACTORS	>= the maximum prime factor of NPASS
 * MAX_PERM	>= the number of prime factors of NPASS.  In addition,
 *	if the square-free portion K of NPASS has two or more prime
 *	factors, then MAX_PERM >= (K-1)
 *
 * storage in FACTOR for a maximum of 15 prime factors of NPASS. if NPASS
 * has more than one square-free factor, the product of the square-free
 * factors must be <= 210 array storage for maximum prime factor of 23 the
 * following two constants should agree with the array dimensions.
 *
 *-----------------------------------------------------------------------*
 *
 * void fft_free (void);
 *
 * free-up allocated temporary storage after finished all the Fourier
 * transforms.
 *
 *----------------------------------------------------------------------*/

#define size_t unsigned int

/* double precision routine
static int fftradix (double Re[], double Im[],
	  size_t nTotal, size_t nPass, size_t nSpan, int isign,
	  int max_factors, int max_perm);*/

/* float precision routine
static int fftradixf (float Re[], float Im[],
	   size_t nTotal, size_t nPass, size_t nSpan, int isign,
	   int max_factors, int max_perm);*/

/* parameters for memory management */

static size_t SpaceAlloced = 0;
static size_t MaxPermAlloced = 0;

/* temp space, (void *) since both float and double routines use it */
static void *Tmp0 = NULL;	/* temp space for real part */
static void *Tmp1 = NULL;	/* temp space for imaginary part */
static void *Tmp2 = NULL;	/* temp space for Cosine values */
static void *Tmp3 = NULL;	/* temp space for Sine values */
static int  *Perm = NULL;	/* Permutation vector */

#define NFACTOR	11
static int factor [NFACTOR];

void fft_free (void)
{
   SpaceAlloced = MaxPermAlloced = 0;
   if (Tmp0 != NULL)	{ free (Tmp0);	Tmp0 = NULL; }
   if (Tmp1 != NULL)	{ free (Tmp1);	Tmp1 = NULL; }
   if (Tmp2 != NULL)	{ free (Tmp2);	Tmp2 = NULL; }
   if (Tmp3 != NULL)	{ free (Tmp3);	Tmp3 = NULL; }
   if (Perm != NULL)	{ free (Perm);	Perm = NULL; }
}

/*
 *
 */
int fftn (int ndim, int *dims, double *Re, double *Im, int iSign)
{
   size_t nSpan, nPass;
   int ret, i , nTotal;
   unsigned int  max_factors,max_perm;
   REAL scaling;

   /*
    * tally the number of elements in the data array
    * and determine the number of dimensions
    */
   nTotal = 1;
   if (ndim && dims [0])
     {
	for (i = 0; i < ndim; i++)
	  {
	     if (dims [i] <= 0)
	       {
		  fputs ("Error: fftn () - dimension error\n", stderr);
		  fft_free ();	/* free-up memory */
		  return -1;
	       }
	     nTotal *= dims [i];
	  }
     }
   else
     {
	ndim = 0;
	for (i = 0; dims [i]; i++)
	  {
	     if (dims [i] <= 0)
	       {
		  fputs ("Error: fftn() - dimension error\n", stderr);
		  fft_free ();	/* free-up memory */
		  return -1;
	       }
	     nTotal *= dims [i];
	     ndim++;
	  }
     }

   /* determine maximum number of factors and permuations */
#if 1
   /*
    * follow John Beale's example, just use the largest dimension and don't
    * worry about excess allocation.  May be someone else will do it?
    */
   max_factors = max_perm = 1;
   for (i = 0; i < ndim; i++)
     {
	nSpan = dims [i];
	if (nSpan > max_factors) max_factors = nSpan;
	if (nSpan > max_perm) max_perm = nSpan;
     }
#endif
   /* loop over the dimensions: */
   nPass = 1;
   for (i = 0; i < ndim; i++)
     {
	printf("passage %d\n", i);
	nSpan = dims [i];
	nPass *= nSpan;
	ret = FFTRADIX (Re, Im, nTotal, nSpan, nPass, iSign,
			max_factors, max_perm);
	/* exit, clean-up already done */
	if (ret)
	  return ret;
     }

   /* Divide through by the normalizing constant: */
	if (iSign < 0) {
  	   scaling = 1.0 / nTotal;	/* multiply is often faster */
	   for (i = 0; i < nTotal; i += 1)
	     {
	       Re [i] *= scaling;
	       Im [i] *= scaling;
	     }
     }
   return 0;
}

/*
 * singleton's mixed radix routine
 *
 * could move allocation out to fftn(), but leave it here so that it's
 * possible to make this a standalone function
 */
static int FFTRADIX (double *Re, double *Im, unsigned int nTotal, unsigned int nPass, unsigned int nSpan, int iSign, int max_factors, int max_perm)
{
   int ii, mfactor, kspan, ispan, inc;
   int j, jc, jf, jj, k, k1, k2, k3, k4, kk, kt, nn, ns, nt;

   REAL radf;
   REAL c1, c2, c3, cd, aa, aj, ak, ajm, ajp, akm, akp;
   REAL s1, s2, s3, sd, bb, bj, bk, bjm, bjp, bkm, bkp;

   REAL *Rtmp = NULL;	/* temp space for real part*/
   REAL *Itmp = NULL;	/* temp space for imaginary part */
   REAL *Cos = NULL;	/* Cosine values */
   REAL *Sin = NULL;	/* Sine values */

   REAL s60 = SIN60;		/* sin(60 deg) */
   REAL c72 = COS72;		/* cos(72 deg) */
   REAL s72 = SIN72;		/* sin(72 deg) */
   REAL pi2 = M_PI;		/* use PI first, 2 PI later */

   /* gcc complains about k3 being uninitialized, but I can't find out where
    * or why ... it looks okay to me.
    *
    * initialize to make gcc happy
    */
   k3 = 0;

   /* gcc complains about c2, c3, s2,s3 being uninitialized, but they're
    * only used for the radix 4 case and only AFTER the (s1 == 0.0) pass
    * through the loop at which point they will have been calculated.
    *
    * initialize to make gcc happy
    */
   c2 = c3 = s2 = s3 = 0.0;

   /* Parameter adjustments, was fortran so fix zero-offset */
   Re--;
   Im--;

   if (nPass < 2)
     return 0;

   /*  allocate storage */
   if (SpaceAlloced < max_factors * sizeof (REAL))
     {
#ifdef SUN_BROKEN_REALLOC
	if (!SpaceAlloced)	/* first time */
	  {
	     SpaceAlloced = max_factors * sizeof (REAL);
	     Tmp0 = (void *) malloc (SpaceAlloced);
	     Tmp1 = (void *) malloc (SpaceAlloced);
	     Tmp2 = (void *) malloc (SpaceAlloced);
	     Tmp3 = (void *) malloc (SpaceAlloced);
	  }
	else
	  {
#endif
	     SpaceAlloced = max_factors * sizeof (REAL);
	     Tmp0 = (void *) realloc (Tmp0, SpaceAlloced);
	     Tmp1 = (void *) realloc (Tmp1, SpaceAlloced);
	     Tmp2 = (void *) realloc (Tmp2, SpaceAlloced);
	     Tmp3 = (void *) realloc (Tmp3, SpaceAlloced);
#ifdef SUN_BROKEN_REALLOC
	  }
#endif
     }
   else
     {
	/* allow full use of alloc'd space */
	max_factors = SpaceAlloced / sizeof (REAL);
     }
   if ((int)MaxPermAlloced < max_perm)
     {
#ifdef SUN_BROKEN_REALLOC
	if (!MaxPermAlloced)	/* first time */
	  Perm = (int *) malloc (max_perm * sizeof(int));
	else
#endif
	  Perm = (int *) realloc (Perm, max_perm * sizeof(int));
	MaxPermAlloced = max_perm;
     }
   else
     {
	/* allow full use of alloc'd space */
	max_perm = MaxPermAlloced;
     }
   if (Tmp0 == NULL || Tmp1 == NULL || Tmp2 == NULL || Tmp3 == NULL
       || Perm == NULL)
     goto Memory_Error_Label;

   /* assign pointers */
   Rtmp = (REAL *) Tmp0;
   Itmp = (REAL *) Tmp1;
   Cos  = (REAL *) Tmp2;
   Sin  = (REAL *) Tmp3;

   /*
    * Function Body
    */
   inc = iSign;
   if (iSign < 0) {
      s72 = -s72;
      s60 = -s60;
      pi2 = -pi2;
      inc = -inc;		/* absolute value */
   }

   /* adjust for strange increments */
   nt = inc * nTotal;
   ns = inc * nSpan;
   kspan = ns;

   nn = nt - inc;
   jc = ns / nPass;
   radf = pi2 * (double) jc;
   pi2 *= 2.0;			/* use 2 PI from here on */

   ii = 0;
   jf = 0;
   /*  determine the factors of n */
   mfactor = 0;
   k = nPass;
   while (k % 16 == 0) {
      mfactor++;
      factor [mfactor - 1] = 4;
      k /= 16;
   }
   j = 3;
   jj = 9;
   do {
      while (k % jj == 0) {
	 mfactor++;
	 factor [mfactor - 1] = j;
	 k /= jj;
      }
      j += 2;
      jj = j * j;
   }
   while (jj <= k);

   if (k <= 4) {
      kt = mfactor;
      factor [mfactor] = k;
      if (k != 1)
	      mfactor++;
    }
   else {
      if (k - (k / 4 << 2) == 0) {
	      mfactor++;
	      factor [mfactor - 1] = 2;
	      k /= 4;
       }
      kt = mfactor;
      j = 2;
      do {
	 if (k % j == 0) {
	    mfactor++;
	    factor [mfactor - 1] = j;
	    k /= j;
	 }
	 j = ((j + 1) / 2 << 1) + 1;
      } while (j <= k);
   }
   if (kt) {
      j = kt;
      do {
	 mfactor++;
	 factor [mfactor - 1] = factor [j - 1];
	 j--;
      } while (j);
   }

   /* test that mfactors is in range */
   if (mfactor > NFACTOR)
     {
	fputs ("Error: fftradix() - exceeded number of factors\n", stderr);
	goto Memory_Error_Label;
      }

   /* compute fourier transform */
   for (;;) {
      sd = radf / (double) kspan;
      cd = sin(sd);
      cd = 2.0 * cd * cd;
      sd = sin(sd + sd);
      kk = 1;
      ii++;

      switch (factor [ii - 1]) {
       case 2:
	 /* transform for factor of 2 (including rotation factor) */
	 kspan /= 2;
	 k1 = kspan + 2;
	 do {
	    do {
	       k2 = kk + kspan;
	       ak = Re [k2];
	       bk = Im [k2];
	       Re [k2] = Re [kk] - ak;
	       Im [k2] = Im [kk] - bk;
	       Re [kk] += ak;
	       Im [kk] += bk;
	       kk = k2 + kspan;
	    } while (kk <= nn);
	    kk -= nn;
	 } while (kk <= jc);
	 if (kk > kspan)
	   goto Permute_Results_Label;		/* exit infinite loop */
	 do {
	    c1 = 1.0 - cd;
	    s1 = sd;
	    do {
	       do {
		  do {
		     k2 = kk + kspan;
		     ak = Re [kk] - Re [k2];
		     bk = Im [kk] - Im [k2];
		     Re [kk] += Re [k2];
		     Im [kk] += Im [k2];
		     Re [k2] = c1 * ak - s1 * bk;
		     Im [k2] = s1 * ak + c1 * bk;
		     kk = k2 + kspan;
		  } while (kk < nt);
		  k2 = kk - nt;
		  c1 = -c1;
		  kk = k1 - k2;
	       } while (kk > k2);
	       ak = c1 - (cd * c1 + sd * s1);
	       s1 = sd * c1 - cd * s1 + s1;
	       c1 = 2.0 - (ak * ak + s1 * s1);
	       s1 *= c1;
	       c1 *= ak;
	       kk += jc;
	    } while (kk < k2);
	    k1 += inc + inc;
	    kk = (k1 - kspan) / 2 + jc;
	 } while (kk <= jc + jc);
	 break;

       case 4:			/* transform for factor of 4 */
	 ispan = kspan;
	 kspan /= 4;

	 do {
	    c1 = 1.0;
	    s1 = 0.0;
	    do {
	       do {
		  k1 = kk + kspan;
		  k2 = k1 + kspan;
		  k3 = k2 + kspan;
		  akp = Re [kk] + Re [k2];
		  akm = Re [kk] - Re [k2];
		  ajp = Re [k1] + Re [k3];
		  ajm = Re [k1] - Re [k3];
		  bkp = Im [kk] + Im [k2];
		  bkm = Im [kk] - Im [k2];
		  bjp = Im [k1] + Im [k3];
		  bjm = Im [k1] - Im [k3];
		  Re [kk] = akp + ajp;
		  Im [kk] = bkp + bjp;
		  ajp = akp - ajp;
		  bjp = bkp - bjp;
		  if (iSign < 0) {
		     akp = akm + bjm;
		     bkp = bkm - ajm;
		     akm -= bjm;
		     bkm += ajm;
		  } else {
		     akp = akm - bjm;
		     bkp = bkm + ajm;
		     akm += bjm;
		     bkm -= ajm;
		  }
		  /* avoid useless multiplies */
		  if (s1 == 0.0) {
		     Re [k1] = akp;
		     Re [k2] = ajp;
		     Re [k3] = akm;
		     Im [k1] = bkp;
		     Im [k2] = bjp;
		     Im [k3] = bkm;
		  } else {
		     Re [k1] = akp * c1 - bkp * s1;
		     Re [k2] = ajp * c2 - bjp * s2;
		     Re [k3] = akm * c3 - bkm * s3;
		     Im [k1] = akp * s1 + bkp * c1;
		     Im [k2] = ajp * s2 + bjp * c2;
		     Im [k3] = akm * s3 + bkm * c3;
		  }
		  kk = k3 + kspan;
	       } while (kk <= nt);

	       c2 = c1 - (cd * c1 + sd * s1);
	       s1 = sd * c1 - cd * s1 + s1;
	       c1 = 2.0 - (c2 * c2 + s1 * s1);
	       s1 *= c1;
	       c1 *= c2;
	       /* values of c2, c3, s2, s3 that will get used next time */
	       c2 = c1 * c1 - s1 * s1;
	       s2 = 2.0 * c1 * s1;
	       c3 = c2 * c1 - s2 * s1;
	       s3 = c2 * s1 + s2 * c1;
	       kk = kk - nt + jc;
	    } while (kk <= kspan);
	    kk = kk - kspan + inc;
	 } while (kk <= jc);
	 if (kspan == jc)
	   goto Permute_Results_Label;		/* exit infinite loop */
	 break;

       default:
	 /*  transform for odd factors */
#ifdef FFT_RADIX4
	 fputs ("Error: " FFTRADIXS "(): compiled for radix 2/4 only\n", stderr);
	 fft_free ();		/* free-up memory */
	 return -1;
	 break;
#else	/* FFT_RADIX4 */
	 k = factor [ii - 1];
	 ispan = kspan;
	 kspan /= k;

	 switch (k) {
	  case 3:	/* transform for factor of 3 (optional code) */
	    do {
	       do {
		  k1 = kk + kspan;
		  k2 = k1 + kspan;
		  ak = Re [kk];
		  bk = Im [kk];
		  aj = Re [k1] + Re [k2];
		  bj = Im [k1] + Im [k2];
		  Re [kk] = ak + aj;
		  Im [kk] = bk + bj;
		  ak -= 0.5 * aj;
		  bk -= 0.5 * bj;
		  aj = (Re [k1] - Re [k2]) * s60;
		  bj = (Im [k1] - Im [k2]) * s60;
		  Re [k1] = ak - bj;
		  Re [k2] = ak + bj;
		  Im [k1] = bk + aj;
		  Im [k2] = bk - aj;
		  kk = k2 + kspan;
	       } while (kk < nn);
	       kk -= nn;
	    } while (kk <= kspan);
	    break;

	  case 5:	/*  transform for factor of 5 (optional code) */
	    c2 = c72 * c72 - s72 * s72;
	    s2 = 2.0 * c72 * s72;
	    do {
	       do {
		  k1 = kk + kspan;
		  k2 = k1 + kspan;
		  k3 = k2 + kspan;
		  k4 = k3 + kspan;
		  akp = Re [k1] + Re [k4];
		  akm = Re [k1] - Re [k4];
		  bkp = Im [k1] + Im [k4];
		  bkm = Im [k1] - Im [k4];
		  ajp = Re [k2] + Re [k3];
		  ajm = Re [k2] - Re [k3];
		  bjp = Im [k2] + Im [k3];
		  bjm = Im [k2] - Im [k3];
		  aa = Re [kk];
		  bb = Im [kk];
		  Re [kk] = aa + akp + ajp;
		  Im [kk] = bb + bkp + bjp;
		  ak = akp * c72 + ajp * c2 + aa;
		  bk = bkp * c72 + bjp * c2 + bb;
		  aj = akm * s72 + ajm * s2;
		  bj = bkm * s72 + bjm * s2;
		  Re [k1] = ak - bj;
		  Re [k4] = ak + bj;
		  Im [k1] = bk + aj;
		  Im [k4] = bk - aj;
		  ak = akp * c2 + ajp * c72 + aa;
		  bk = bkp * c2 + bjp * c72 + bb;
		  aj = akm * s2 - ajm * s72;
		  bj = bkm * s2 - bjm * s72;
		  Re [k2] = ak - bj;
		  Re [k3] = ak + bj;
		  Im [k2] = bk + aj;
		  Im [k3] = bk - aj;
		  kk = k4 + kspan;
	       } while (kk < nn);
	       kk -= nn;
	    } while (kk <= kspan);
	    break;

	  default:
	    if (k != jf) {
	       jf = k;
	       s1 = pi2 / (double) k;
	       c1 = cos(s1);
	       s1 = sin(s1);
	       if (jf > max_factors)
		 goto Memory_Error_Label;
	       Cos [jf - 1] = 1.0;
	       Sin [jf - 1] = 0.0;
	       j = 1;
	       do {
		  Cos [j - 1] = Cos [k - 1] * c1 + Sin [k - 1] * s1;
		  Sin [j - 1] = Cos [k - 1] * s1 - Sin [k - 1] * c1;
		  k--;
		  Cos [k - 1] = Cos [j - 1];
		  Sin [k - 1] = -Sin [j - 1];
		  j++;
	       } while (j < k);
	    }
	    do {
	       do {
		  k1 = kk;
		  k2 = kk + ispan;
		  ak = aa = Re [kk];
		  bk = bb = Im [kk];
		  j = 1;
		  k1 += kspan;
		  do {
		     k2 -= kspan;
		     j++;
		     Rtmp [j - 1] = Re [k1] + Re [k2];
		     ak += Rtmp [j - 1];
		     Itmp [j - 1] = Im [k1] + Im [k2];
		     bk += Itmp [j - 1];
		     j++;
		     Rtmp [j - 1] = Re [k1] - Re [k2];
		     Itmp [j - 1] = Im [k1] - Im [k2];
		     k1 += kspan;
		  } while (k1 < k2);
		  Re [kk] = ak;
		  Im [kk] = bk;
		  k1 = kk;
		  k2 = kk + ispan;
		  j = 1;
		  do {
		     k1 += kspan;
		     k2 -= kspan;
		     jj = j;
		     ak = aa;
		     bk = bb;
		     aj = 0.0;
		     bj = 0.0;
		     k = 1;
		     do {
			k++;
			ak += Rtmp [k - 1] * Cos [jj - 1];
			bk += Itmp [k - 1] * Cos [jj - 1];
			k++;
			aj += Rtmp [k - 1] * Sin [jj - 1];
			bj += Itmp [k - 1] * Sin [jj - 1];
			jj += j;
			if (jj > jf) {
			   jj -= jf;
			}
		     } while (k < jf);
		     k = jf - j;
		     Re [k1] = ak - bj;
		     Im [k1] = bk + aj;
		     Re [k2] = ak + bj;
		     Im [k2] = bk - aj;
		     j++;
		  } while (j < k);
		  kk += ispan;
	       } while (kk <= nn);
	       kk -= nn;
	    } while (kk <= kspan);
	    break;
	 }
	 /*  multiply by rotation factor (except for factors of 2 and 4) */
	 if (ii == mfactor)
	   goto Permute_Results_Label;		/* exit infinite loop */
	 kk = jc + 1;
	 do {
	    c2 = 1.0 - cd;
	    s1 = sd;
	    do {
	       c1 = c2;
	       s2 = s1;
	       kk += kspan;
	       do {
		  do {
		     ak = Re [kk];
		     Re [kk] = c2 * ak - s2 * Im [kk];
		     Im [kk] = s2 * ak + c2 * Im [kk];
		     kk += ispan;
		  } while (kk <= nt);
		  ak = s1 * s2;
		  s2 = s1 * c2 + c1 * s2;
		  c2 = c1 * c2 - ak;
		  kk = kk - nt + kspan;
	       } while (kk <= ispan);
	       c2 = c1 - (cd * c1 + sd * s1);
	       s1 += sd * c1 - cd * s1;
	       c1 = 2.0 - (c2 * c2 + s1 * s1);
	       s1 *= c1;
	       c2 *= c1;
	       kk = kk - ispan + jc;
	    } while (kk <= kspan);
	    kk = kk - kspan + jc + inc;
	 } while (kk <= jc + jc);
	 break;
#endif	/* FFT_RADIX4 */
      }
   }

/*  permute the results to normal order---done in two stages */
/*  permutation for square factors of n */
Permute_Results_Label:
   Perm [0] = ns;
   if (kt) {
      k = kt + kt + 1;
      if (mfactor < k)
	k--;
      j = 1;
      Perm [k] = jc;
      do {
	 Perm [j] = Perm [j - 1] / factor [j - 1];
	 Perm [k - 1] = Perm [k] * factor [j - 1];
	 j++;
	 k--;
      } while (j < k);
      k3 = Perm [k];
      kspan = Perm [1];
      kk = jc + 1;
      k2 = kspan + 1;
      j = 1;
      if (nPass != nTotal) {
/*  permutation for multivariate transform */
Permute_Multi_Label:
	 do {
	    do {
	       k = kk + jc;
	       do {
		  /* swap Re [kk] <> Re [k2], Im [kk] <> Im [k2] */
		  ak = Re [kk]; Re [kk] = Re [k2]; Re [k2] = ak;
		  bk = Im [kk]; Im [kk] = Im [k2]; Im [k2] = bk;
		  kk += inc;
		  k2 += inc;
	       } while (kk < k);
	       kk += ns - jc;
	       k2 += ns - jc;
	    } while (kk < nt);
	    k2 = k2 - nt + kspan;
	    kk = kk - nt + jc;
	 } while (k2 < ns);
	 do {
	    do {
	       k2 -= Perm [j - 1];
	       j++;
	       k2 = Perm [j] + k2;
	    } while (k2 > Perm [j - 1]);
	    j = 1;
	    do {
	       if (kk < k2)
		 goto Permute_Multi_Label;
	       kk += jc;
	       k2 += kspan;
	    } while (k2 < ns);
	 } while (kk < ns);
      } else {
/*  permutation for single-variate transform (optional code) */
Permute_Single_Label:
	 do {
	    /* swap Re [kk] <> Re [k2], Im [kk] <> Im [k2] */
	    ak = Re [kk]; Re [kk] = Re [k2]; Re [k2] = ak;
	    bk = Im [kk]; Im [kk] = Im [k2]; Im [k2] = bk;
	    kk += inc;
	    k2 += kspan;
	 } while (k2 < ns);
	 do {
	    do {
	       k2 -= Perm [j - 1];
	       j++;
	       k2 = Perm [j] + k2;
	    } while (k2 > Perm [j - 1]);
	    j = 1;
	    do {
	       if (kk < k2)
		 goto Permute_Single_Label;
	       kk += inc;
	       k2 += kspan;
	    } while (k2 < ns);
	 } while (kk < ns);
      }
      jc = k3;
   }

   if ((kt << 1) + 1 >= mfactor)
     return 0;
   ispan = Perm [kt];
   /* permutation for square-free factors of n */
   j = mfactor - kt;
   factor [j] = 1;
   do {
      factor [j - 1] *= factor [j];
      j--;
   } while (j != kt);
   kt++;
   nn = factor [kt - 1] - 1;
   if (nn > max_perm)
     goto Memory_Error_Label;
   j = jj = 0;
   for (;;) {
      k = kt + 1;
      k2 = factor [kt - 1];
      kk = factor [k - 1];
      j++;
      if (j > nn)
	break;				/* exit infinite loop */
      jj += kk;
      while (jj >= k2) {
	 jj -= k2;
	 k2 = kk;
	 k++;
	 kk = factor [k - 1];
	 jj += kk;
      }
      Perm [j - 1] = jj;
   }
   /*  determine the permutation cycles of length greater than 1 */
   j = 0;
   for (;;) {
      do {
	 j++;
	 kk = Perm [j - 1];
      } while (kk < 0);
      if (kk != j) {
	 do {
	    k = kk;
	    kk = Perm [k - 1];
	    Perm [k - 1] = -kk;
	 } while (kk != j);
	 k3 = kk;
      } else {
	 Perm [j - 1] = -j;
	 if (j == nn)
	   break;		/* exit infinite loop */
      }
   }
   max_factors *= inc;
   /*  reorder a and b, following the permutation cycles */
   for (;;) {
      j = k3 + 1;
      nt -= ispan;
      ii = nt - inc + 1;
      if (nt < 0)
	break;			/* exit infinite loop */
      do {
	 do {
	    j--;
	 } while (Perm [j - 1] < 0);
	 jj = jc;
	 do {
	    kspan = jj;
	    if (jj > max_factors) {
	       kspan = max_factors;
	    }
	    jj -= kspan;
	    k = Perm [j - 1];
	    kk = jc * k + ii + jj;
	    k1 = kk + kspan;
	    k2 = 0;
	    do {
	       k2++;
	       Rtmp [k2 - 1] = Re [k1];
	       Itmp [k2 - 1] = Im [k1];
	       k1 -= inc;
	    } while (k1 != kk);
	    do {
	       k1 = kk + kspan;
	       k2 = k1 - jc * (k + Perm [k - 1]);
	       k = -Perm [k - 1];
	       do {
		  Re [k1] = Re [k2];
		  Im [k1] = Im [k2];
		  k1 -= inc;
		  k2 -= inc;
	       } while (k1 != kk);
	       kk = k2;
	    } while (k != j);
	    k1 = kk + kspan;
	    k2 = 0;
	    do {
	       k2++;
	       Re [k1] = Rtmp [k2 - 1];
	       Im [k1] = Itmp [k2 - 1];
	       k1 -= inc;
	    } while (k1 != kk);
	 } while (jj);
      } while (j != 1);
   }
   return 0;			/* exit point here */

   /* alloc or other problem, do some clean-up */
Memory_Error_Label:
   fputs ("Error: fftradix() - insufficient memory.\n", stderr);
   fft_free ();			/* free-up memory */
   return -1;
}
    /* _FFTN_C */
/* ---------------------- end-of-file (c source) ---------------------- */


void test_conv(void)
{
/*	int id,k,n;
	float *x;
	float *vals;
	float *ans,*resp;

	n = 90*2;

	x = (float*)malloc(sizeof(float)*n);
	vals = (float*)malloc(sizeof(float)*n);
	ans = (float*)malloc(sizeof(float)*n);
	resp = (float*)malloc(sizeof(float)*n);

	n /=2;

	id=OpenGraph();
	ClearGraph(id);

	for(k=0 ; k<n ; k++)
	{
		x[k] = k+1;
	}

	for(k=0 ; k<30 ; k++)
	{
		vals[k] = 0.0;
	}
	for(k=30 ; k<60 ; k++)
	{
		vals[k] = 10.0;
	}
	for(k=60 ; k<90 ; k++)
	{
		vals[k] = 0.0;
	}

	for(k=0 ; k<n ; k++)
	{
		resp[k] = vals[k];
	}

	convlv(vals, n, resp, 45, 1, ans);

	PLOT_LINE(id, x, ans, n, "X", "Y", "", "Test convolution");

	return;*/
}

/**************************************************/


/*---------------------------------------------------------------------*/


/***********************************************************************
** -- init_simplex ()
**
**     Preparation des tableaux necessaires au simplex
**
**   Utilisation : donner a la fonction un pointeur fun sur une structure
** deja existante de type Parametre (la structure doit deja etre allouee)
**
**                 le nombre de points de la fonction dentree
**
**                 un pointeur sur un vertex (en fait un tableau de
**                     vertices deja alloue)
**
**                 le nombre de vertices
**
**                 le nombre de parametres a estimer dans chaque vertex
**
*************************************************************************/
void init_simplex (Param_simplex *fun, int nb_pts, Vertex *vtx, int nb_vtx, int nb_param)
{

	int i;

	/* Allocation memoire de la structure Parametre suivant le nb de pts de la courbe voulu */
	fun->nb_point = nb_pts;
	if ((fun->tab_val_in = (float *) malloc(nb_pts * sizeof(float))) == NULL) exit(-1);
	if ((fun->tab_val_out = (float *) malloc(nb_pts * sizeof(float))) == NULL) exit(-1);
	if ((fun->x = (float *) malloc(nb_pts * sizeof(float))) == NULL) exit(-1);


	/* Allocation des tableaux des vertices */
	for (i = 0 ; i < nb_vtx ; i++) {
		vtx[i].nb_valeur_estimer = nb_param;
		if ((vtx[i].valeur_estimer = (float *) malloc(vtx[i].nb_valeur_estimer * sizeof(float))) == NULL) exit(-2);
	}


}

/***************************************************************************
** -- free_simplex()
**
**     Libere les tableaux dynamiques utilises par le simplex
**
**     Utilisation : A la fin de la routine simplex
**
**                   passer un pointeur sur les Parametre de la fonction a minimiser
**
**                   passer un tableau de vertices
**
**                   et le nombre de vertices de ce tableau
**
****************************************************************************/
void free_simplex(Param_simplex *fun, Vertex *vtx, int nb_vtx)
{

	int i;

	/* Liberation des courbes */
	if (fun->tab_val_in) free(fun->tab_val_in);
	if (fun->tab_val_out) free(fun->tab_val_out);
	if (fun->x) free(fun->x);

	/* Liberation des vertices */
	for (i = 0 ; i < nb_vtx ; i++)
		if (&vtx[i])
			if (vtx[i].valeur_estimer) free(vtx[i].valeur_estimer);

}
/*********************************************************************
** -- swap_vtx()
**
**      Echange les donnees de 2 vertices passes en arguments
**
**********************************************************************/
void swap_vtx(Vertex *v1, Vertex *v2)
{

	float tmp;
	int i;

	for (i = 0 ; i < v1->nb_valeur_estimer ; i++) {
		tmp = v1->valeur_estimer[i];
		v1->valeur_estimer[i] = v2->valeur_estimer[i];
		v2->valeur_estimer[i] = tmp;
	}

	tmp = v1->r;
	v1->r = v2->r;
	v2->r = tmp;

}

/***********************************************************************
** -- simplex_gen()
**
**     A partir d'un panel de solutions donnees par lutilisateur, le simplex
**   va rechercher une solution qui minimise la difference
**
**   Utilisation :
**
**        Il faut lui fournir des Parametres de fonction remplis, des vertices
**      initialises avec les premieres valeurs , la fonction de minimisation,
**      et le nombre de vtx
**
************************************************************************/
void simplex_gen(void (*function) (/* ??? */), Param_simplex f, Vertex *vtx, int nbv)
{

	int t, i, j;
	Vertex vtx_ref, vtx_ctr, vtx_exp, vtx_mid;
	int dimension = vtx[0].nb_valeur_estimer;
	int nb_vertices = nbv;

	vtx_ref.nb_valeur_estimer = dimension;
	vtx_ref.valeur_estimer = (float *) malloc(vtx_ref.nb_valeur_estimer * sizeof(float));

	vtx_ctr.nb_valeur_estimer = dimension;
	vtx_ctr.valeur_estimer = (float *) malloc(vtx_ctr.nb_valeur_estimer * sizeof(float));

	vtx_exp.nb_valeur_estimer = dimension;
	vtx_exp.valeur_estimer = (float *) malloc(vtx_exp.nb_valeur_estimer * sizeof(float));

	vtx_mid.nb_valeur_estimer = dimension;
	vtx_mid.valeur_estimer = (float *) malloc(vtx_mid.nb_valeur_estimer * sizeof(float));

	for (j = 0 ; j < nb_vertices ; j++) (*function)(&vtx[j], f);

	for (t = 0 ; t < 100 ; t++) {

		/*
		for (i = 0 ; i < nb_vertices ; i++)
			printf("vtx[%d].r = %10.5f\n", i, vtx[i].r);
		*/

		/* Ordonner les vertices par resultat croissant ! */
		for(i = 0 ; i < nb_vertices ; i++)
			for (j = 0 ; j < nb_vertices - i -1; j++)
				if (vtx[j].r > vtx[j+1].r) swap_vtx(&vtx[j], &vtx[j+1]);

		/*
		for (i = 0 ; i < nb_vertices ; i++)
			printf("vtx[%d].r = %10.5f\n", i, vtx[i].r);
		*/

		/* Calcul du point de reflection */
		for (j = 0 ; j < dimension ; j++) {
			vtx_mid.valeur_estimer[j] = 0;
			for (i = 0 ; i < nb_vertices ; i++) vtx_mid.valeur_estimer[j] += vtx[i].valeur_estimer[j];
			vtx_mid.valeur_estimer[j] /= nb_vertices;
		}

		/*
		printf("--------- mid vtx ----------\n");
		for (i = 0 ; i < vtx[0].nb_valeur_estimer ; i++)
			printf("vtx_mid.tab_valeur_estimer[%d] = %10.5f\n", i, vtx_mid.valeur_estimer[i]);
		*/

		/* Calcul du reflected vertex */
		for (j = 0 ; j < dimension ; j++)
			vtx_ref.valeur_estimer[j] = 2* vtx_mid.valeur_estimer[j] - vtx[nb_vertices-1].valeur_estimer[j];
		(*function)(&vtx_ref, f);

		/*
		printf("--------- ref vtx ----------\n");
		for (i = 0 ; i < vtx[0].nb_valeur_estimer ; i++)
			printf("vtx_ref.tab_valeur_estimer[%d] = %10.5f\n", i, vtx_ref.valeur_estimer[i]);
		printf("vtx_ref.r                     = %10.5f\n\n", vtx_ref.r);
		*/

		/* Si le reflected vertex a donne de meilleurs resultat
		a la minimisation que le precedent meilleur */
		if (vtx_ref.r < vtx[0].r) {
			/* Calcul de l'expanded vertex */
			for (j=  0 ; j < dimension ; j++)
				vtx_exp.valeur_estimer[j] = 3 * vtx_mid.valeur_estimer[j] - 2 * vtx[nb_vertices-1].valeur_estimer[j];
			(*function)(&vtx_exp, f);

			/*
			printf("--------- exp vtx ----------\n");
			for (i = 0 ; i < vtx[0].nb_valeur_estimer ; i++)
				printf("vtx_exp.tab_valeur_estimer[%d] = %10.5f\n", i, vtx_exp.valeur_estimer[i]);
			printf("vtx_exp.r                     = %10.5f\n\n", vtx_exp.r);
			*/

			/* Si le nouveau vertex est meilleur que le reflected */
			if (vtx_exp.r < vtx_ref.r) {
				/* on place l'expanded vtx a la fin de la liste des vtx */
				for (j = 0 ; j < dimension ; j++)
					vtx[nb_vertices-1].valeur_estimer[j] = vtx_exp.valeur_estimer[j];
				vtx[nb_vertices-1].r = vtx_exp.r;
			}
			else {
			  	/* On place le ref vtx a la place du dernier vtx */
				for (j = 0 ; j < dimension ; j++)
					vtx[nb_vertices-1].valeur_estimer[j] = vtx_ref.valeur_estimer[j];
				vtx[nb_vertices-1].r = vtx_ref.r;
			}
		}
		/* le ref vtx est moins bon que le premier vtx */
		else {
			if (vtx_ref.r < vtx[nb_vertices-1].r) {
				for (j = 0 ; j < dimension ; j++)
					vtx[nb_vertices-1].valeur_estimer[j] = vtx_ref.valeur_estimer[j];
				vtx[nb_vertices-1].r = vtx_ref.r;
			}
			else {
				/* Calcul du vertex contracted */
				for (j = 0 ; j < dimension ; j++)
					vtx_ctr.valeur_estimer[j] = vtx_mid.valeur_estimer[j]/2 + vtx[nb_vertices-1].valeur_estimer[j]/2;
				(*function)(&vtx_ctr, f);

				/*
				printf("--------- ctr vtx ----------\n");
				for (i = 0 ; i < vtx[0].nb_valeur_estimer ; i++)
					printf("vtx_ctr.tab_valeur_estimer[%d] = %10.5f\n", i, vtx_ctr.valeur_estimer[i]);
				printf("vtx_ctr.r                     = %10.5f\n\n", vtx_ctr.r);
				*/

				if (vtx_ctr.r < vtx[nb_vertices-1].r) {
					for (j = 0 ; j < dimension ; j++)
						vtx[nb_vertices-1].valeur_estimer[j] = vtx_ctr.valeur_estimer[j];
					vtx[nb_vertices-1].r = vtx_ctr.r;
				}
				else {
					/* contracter tous les vertices sauf le meilleur (le premier !) */
					for (i = 1; i < nb_vertices; i++) {
						for (j = 0 ; j < dimension ; j++)
							vtx[i].valeur_estimer[j] = (vtx[i].valeur_estimer[j] - vtx[0].valeur_estimer[j]) / 2;
					}
				}
			}
		}
	}

	free(vtx_ref.valeur_estimer);
	free(vtx_ctr.valeur_estimer);
	free(vtx_exp.valeur_estimer);
	free(vtx_mid.valeur_estimer);

}


/***** Pour OTSU entre autre    ****/
/********************************************************************/
/*!
**	calcul de omega   
**
**	\param prob : distribution de proba.
**	\param k1,k2 :  \f$omega =\sum_{k1\le i \le k2}(prob[i])   \f$
**	\retval omega
**
********************************************************************/
double calc_omega(double *prob, int k1, int k2)
{
  int i;
  double omega; 

  omega=0.0;

  for(i=k1;i<=k2;i++)
     omega+=prob[i];

  return(omega);
 
 }


/*******************************************
 ** --  ynist() ----------------------------
**
** Initialisation de l'histogramme (reso)
********************************************/
float   ynist(float aminh, float amaxh, long int Nbpas)
{
 float reso;

  if (Nbpas<=0) return(0.0);
  reso=(float) (amaxh-aminh)/(Nbpas);
  return(reso);

}
/*******************************************
** --  ystog() ----------------------------
**
** Calcul du tableau ibuf pour l'histogramme
**
**
**   Nbpas : Nombre max de points dans le tableau ibuf
**
********************************************/
int  ystog(float *ibuf, long int Nbpas, float aminh, float reso, float valeur, float duree)
{
  int i,NbpasM1;
 
  NbpasM1=Nbpas-1;
  i=(int)floor((valeur-aminh)/reso);
  i=MINI(NbpasM1,MAXI(0,i));
  ibuf[i]=ibuf[i]+duree;
  return(1);
}


/********************************************************************/
/*!
**	calcul de mu   
**
**	\param prob : distribution de proba.
**	\param k1,k2 :  \f$mu =(\sum_{k1\le i \le k2}(i*prob[i]))/omega   \f$
**	\retval mu
**
********************************************************************/
double calc_mu(double *prob, int k1, int k2)
{
  int i;
  double mu; 
  double omega;

  mu=0.0;
  omega=0.0;

  omega=calc_omega(prob,k1,k2);

  for(i=k1;i<=k2;i++)
     mu+=(double)i*prob[i];


   if(omega!=0.0)
     mu/=omega;

       else
         mu=0.0;

  return(mu);
 
 }

/*--------------------   get_gaussian              ---------------------*/
/*----------------------------------------------------------------------*/

/* Calculate a vectorial Gaussian mask with a standard deviation of sigma. */
/* Cut off when tails reach 1% of the maximum value. */
void get_gaussian(float s, float *y, int *len)
{
  int r,i;
  float gaussian_r;
   
  r=1;
  gaussian_r=1;
  *len=0;
  y[*len]=gaussian_r;
  while(gaussian_r >= 0.01){
	gaussian_r=(float)exp((double)(-0.5*(r*r)/(s*s)));
	if (gaussian_r>=0.01){
          for (i=*len;i>=0;i--) y[i+1]=y[i];
 	  i=*len;
	  y[i+2]=gaussian_r;
	  y[0]=gaussian_r;
	  r=r+1;
	  *len=i+2;
	
	}/*end if*/
  
  }/*end while*/

}/*end get_gaussian*/


/*--------------------   get_derigaussian          ---------------------*/
/*----------------------------------------------------------------------*/

/* Calculate a vectorial derivate of a Gaussian mask with a standar deviation */
/* of sigma. Cut off when tails reach 1% of the maximum value.*/

void get_derigaussian(float s, float *y, int *len)
{
  int r,i/*,j*/;
  float deriv_gaussian_r;
  float max_val;

  r=1;
  max_val=0;
  *len=0;
  deriv_gaussian_r=max_val;
  y[*len]=deriv_gaussian_r;
  while (deriv_gaussian_r>=0.01*max_val){
	deriv_gaussian_r=(float)(r/((s*s))*exp((double)(-0.5*(r*r)/(s*s))));
	  for(i=*len;i>=0;i--) y[i+1]=y[i];
	  i=*len;
	  y[i+2]= (- deriv_gaussian_r);
	  y[0]=deriv_gaussian_r;
	  r+=1;
	  *len=i+2;
	  if (deriv_gaussian_r > max_val) max_val=deriv_gaussian_r;

  }/*end while*/
  
 /* top = max_val;*/
}/*end get_derigaussian*/



