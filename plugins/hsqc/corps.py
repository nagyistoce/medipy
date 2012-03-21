#import segmentation as seg
from medipy.base import Image
import numpy as n
from medipy.io import load
import os
import fonc_util
from scipy.optimize import leastsq
import logging
import sys

def interquartile_range(values,scaletonormal=False):
    """
    Computes the interquartile range for the provided sequence of values, a
    more robust estimator than the variance.
    """
    from scipy.stats import scoreatpercentile
    from scipy.special import erfinv
     
    x = n.array(values,copy=False).ravel()
    res = scoreatpercentile(x,75) - scoreatpercentile(x,25)
    if scaletonormal:
        nrm = 8**0.5*erfinv(.5)
        return res/nrm
    else:
        return res

def test():
	root='/base_image/nmr/colon-HSQC/47hutann-C-HSQC/4/pdata/1'
	lolo= load(os.path.join(root,"2rr"))
	H=lolo.data
	M=n.size(H[:,1])
	N=n.size(H[1,:])
	im1=n.zeros([M,N],dtype=n.single)
	#im2=1*n.ones([2000,2000],dtype=n.single)
	#im2[9,9]=2;
	Ima1=Image(data=im1)
	Ima2=Image(data=H)
	#print im2
	#seg.threshold(Ima2,5,1,5,7,Ima1)
	im=Ima1.data
	ind=fonc_util.find(im,im==1)
	#M=n.size(lolo.data[:,1])
	#N=n.size(lolo.data[1,:])
	#M1=n.size(Ima1.data[:,1])
	#N1=n.size(Ima1.data[1,:])
	#print M,N,M1,N1
	#test=lolo.data*Ima1.data
	#ind=fonc_util.find(test,test<0)
	#print ind[0]
	#print im[ind[0][0],ind[1][0]]
	for m in range(20):
		#try:			
			zayneb=H[ind[0][m]-5:ind[0][m]+6,ind[1][m]-5:ind[1][m]+6]*1.
			#print ind
			#break
			#print H[ind[0][m],ind[1][m]]
			#print zayneb[5,5]
			zz=zayneb*1.
			k=fonc_util.subpix(zz,ind[0][m],ind[1][m])
			#z1=fonc_util.cauchy(11,[0.3,0.3])
			zo=zayneb*1.
			zo[zo<0]=z1[zo<0]
			e = lambda v,zo,: n.sum(n.abs(zo-zo[5,5]*fonc_util.cauchy(11,v)),1)
			vi=[0.3,0.3]
			#z[z<0]=0
			v, success = leastsq(e, vi, args=(zo), maxfev=10000)
			#print v
#	v=np.abs(v)
		#except:
			#print z
			#break
	
