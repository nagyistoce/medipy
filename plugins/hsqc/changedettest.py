"""
Validate
"""
import svdd as dd
import corps
import logging
import sys
import warnings
import numpy as np
from medipy.io import load
import os
#from sklearn import svm
import fonc_util as f
#import segmentation as seg
from medipy.base import Image
import gen_in
import lawij
import ser_xml
import MCMC as m
import xml_gen
import locale
#from sklearn import metrics
#from sklearn.cluster import KMeans
import func
#if __name__=="__main__":
    #root='/home/miv/belghith/Bureau/test/PH/5/pdata/1'
def detfark(root):    
    tabp,tabc=gen_in.ini(root)
    inm,inp=f.ser_mat(tabp)
    
    Xt = np.random.rand(10, np.size(tabp,0))
    Xo = np.random.rand(10, np.size(tabp,0))+5
    iter=-1
    res=[]
    for i in range(len(inm)):
        r1=inp[i]
        xti=Xt[:,r1[0]:r1[len(r1)-1]+1]
        xti[5]=0
        xoi=Xo[:,r1[0]:r1[len(r1)-1]+1]
        z=f.find(xti,xti[:,0]!=0)
        xti=xti[z,:]
        xti=xti[0]
        z=f.find(xoi,xoi[:,0]!=0)
        xoi=xoi[z,:]
        xoi=xoi[0]
        #print xti,xoi,np.size(xti,0)
        if np.size(xti,0)>4:
            s=np.std(xti,0)
            S=np.kron(np.ones((np.size(xti,0),1)),s)
            xti=xti/S
            S=np.kron(np.ones((np.size(xoi,0),1)),s)
            xoi=xoi/S
            resini=func.k_means(xti)['clusters']
            #print resini
            if np.minimum(len(resini[0]),len(resini[1]))>0:             
                if float(np.minimum(len(resini[0]),len(resini[1])))/ float(np.size(xti,0))<0.12:
                    #print float(np.minimum(len(resini[0]),len(resini[1])))/ float(np.size(xti,0))
                    if len(resini[0])<len(resini[1]):
                        nxti=xti[resini[1],:]
                    else:
                        xti=xti[resini[0],:]
            s=np.std(xti,0)
            S=np.kron(np.ones((np.size(xti,0),1)),s)
            xti=xti/S        #print nxti, xti     
            S=np.kron(np.ones((np.size(xoi,0),1)),s) 
            xoi=xoi/S
            if np.size(xti,1)>0:
                ctm=0
                me=0
                ram=0
                for ii in range(np.size(xti,1)):
                    X1 = np.array(xti[:,ii])
                    X2 = np.array(xoi[:,ii])
                    X1=f.tab2tab(X1)
                    X2=f.tab2tab(X2)
                    X1=X1.T
                    X2=X2.T
                    #raw_input()
                    C = 5
                    #clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
                    #clf.fit(X1.T)
                    K1 = np.dot(X1.T, X1) # the linear kernel
                    K2 = np.dot(X2.T, X2) # the linear kernel
                    m=0
                    model1 = dd.svdd(X1, K1, C)
                    mX2=[]
                    for m in range(np.size(X2,1)):
                        #y_pred = clf.predict(X2[:, m])
                        diff = dd.svdd_test(model1, X2[:, m])
                        #print diff[0]
                        mX2.append((diff[0]))
                        
                    moy=np.mean(mX2)
                    med=np.mean(mX2)
                    #print (0.5*(moy+med))
                    if (0.5*(moy+med)>2) :
                        ctm+=1
                        me=me+(moy+med)/2
                        ram=ram+np.mean(X1)/np.mean(X2)
                        #print ram,me,ctm
                        #raw_input()
                #print ctm/np.size(xti,1)
                if ctm/np.size(xti,1)>0.3:
                    #print 'ok',func.flou(2,1.05,2.3,4)
                    iter+=1
                    if np.mean(X1)>np.mean(X2):
                        eta='dec'
                    else:
                        eta='inc'
                    flo=0.4*func.flou(me/ctm,10,200,500)+0.6*func.flou(ram/ctm,0.8,1.2,4)    
                    res.append((inm[i],flo,eta))

            print res
            #raw_input()
