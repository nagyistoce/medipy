"""
Functions used in Metabolite annotation
"""
import numpy as np
from scipy.interpolate import interpolate
from scipy.optimize import leastsq
def find(a,cond):
    b=np.nonzero(cond)
    return b
    
def cauchy(p2,p3):
    siz=(p2-1)/2*(np.ones(2))
    xxrange = np.arange(-siz[1],siz[1]+1)
    yyrange = np.arange(-siz[1],siz[1]+1)
    X,Y = np.meshgrid(xxrange,yyrange)
    arg=((1/(p3[0]*3.14159))/((1/p3[0])**2+(X*X)))*((1/(p3[1]*3.14159))/((1/p3[1])**2+(Y*Y)))
    eps=2.2204*10**(-16)
    h=arg
    h[h<(eps*np.amax(h))]=0
    sumh=np.sum(h)
    if sumh!=0:
        h=h/sumh
    h=h/np.amax(h)
    return h
    
def expon(p2,p3):
    siz=(p2-1)/2*(np.ones(2))
    xxrange = np.arange(-siz[1],siz[1]+1)
    yyrange = np.arange(-siz[1],siz[1]+1)
    X,Y = np.meshgrid(xxrange,yyrange)
    #arg=((1/p3[0])/((1/p3[0])**2+(X*X)))*((1/p3[1])/((1/p3[1])**2+(Y*Y)))
    arg=(1/6.28)*(1/p3[0])*np.exp(-X*X/(2*p3[0]**2))*(1/p3[1])*np.exp(-X*X/(2*p3[1]**2))
    eps=2.2204*10**(-16)
    h=arg
    h[h<(eps*np.amax(h))]=0
    sumh=np.sum(h)
    if sumh!=0:
        h=h/sumh
    h=h/np.amax(h)
    return h

def exponcauchy(p2,p3):
    siz=(p2-1)/2*(np.ones(2))
    xxrange = np.arange(-siz[1],siz[1]+1)
    yyrange = np.arange(-siz[1],siz[1]+1)
    X,Y = np.meshgrid(xxrange,yyrange)
    #arg=((1/p3[0])/((1/p3[0])**2+(X*X)))*((1/p3[1])/((1/p3[1])**2+(Y*Y)))
    arg=(1/3.14159)*(1/p3[0])*np.exp(-X*X/(2*p3[0]**2))*((1/(p3[1]*3.14159))/((1/p3[1])**2+(Y*Y)))
    eps=2.2204*10**(-16)
    h=arg
    h[h<(eps*np.amax(h))]=0
    sumh=np.sum(h)
    if sumh!=0:
        h=h/sumh
    h=h/np.amax(h)
    return h
def cauchyexpon(p2,p3):
    siz=(p2-1)/2*(np.ones(2))
    xxrange = np.arange(-siz[1],siz[1]+1)
    yyrange = np.arange(-siz[1],siz[1]+1)
    X,Y = np.meshgrid(xxrange,yyrange)
    #arg=((1/p3[0])/((1/p3[0])**2+(X*X)))*((1/p3[1])/((1/p3[1])**2+(Y*Y)))
    arg=((1/(p3[0]*3.14159))/((1/p3[0])**2+(X*X)))*(1/3.14159)*(1/p3[1])*np.exp(-Y*Y/(2*p3[0]**2))
    eps=2.2204*10**(-16)
    h=arg
    h[h<(eps*np.amax(h))]=0
    sumh=np.sum(h)
    if sumh!=0:
        h=h/sumh
    h=h/np.amax(h)
    return h

def subpix2(z,ii,jj):
    trange = np.arange(7)
    ttrange = np.arange(7)
    X,Y = np.meshgrid(trange,ttrange)
    outgrid = interpolate.interp2d(X,Y,z,kind='quintic')
    xx=yy=np.arange(61)/10.
    l=outgrid(xx,yy)
    l=l[30-9:30+10,30-9:30+10]
    ind=find(l,l==np.amax(l))
    #print l
    #print ind[0][0],ind[1][0]
    ni=ii+(ind[0][0]-9.)/10
    nj=jj+(ind[1][0]-9.)/10
    #print ii,jj
    #print ni,nj
    return[ni,nj]

def dephcl(z):
    e = lambda v,z,: np.sum(np.abs(z-z[9,9]*cauchy(19,v)),1)
    vi=[0.3,0.3]    
    v, success = leastsq(e, vi, args=(z), maxfev=1000)
    if v[0]<0.001 or v[0]>2 or v[1]<0.001 or v[1]>2 :
        v[0]=v[1]=0.3+np.random.normal(0, 0.05, 1)
    return v

def dephcg(z):
    e = lambda v,z,: np.sum(np.abs(z-z[9,9]*expon(19,v)),1)
    vi=[1,1]
    #z[z<0]=0
    v, success = leastsq(e, vi, args=(z), maxfev=1000)
    if v[0]<0.1 or v[0]>4 or v[1]<0.1 or v[1]>4 :
        v[0]=v[1]=2+np.random.normal(0, 0.05, 1)

    return v

def dephcaprio(z,a,b,c):
    if c[0]=='g' and c[1]=='g':
            e = lambda v,z,: np.sum(np.abs(z-z[9,9]*expon(19,v)),1)
            vi=[a,b]
            #z[z<0]=0
            v, success = leastsq(e, vi, args=(z), maxfev=1000)
            if np.abs(float(v[0]-a))>1 or v[0]<0.5 or v[0]>6: 
                v[0]=a+np.random.normal(0, 0.05, 1)
            if np.abs(float(v[1]-b))>1 or v[0]<0.5 or v[0]>6:     
                v[1]=b+np.random.normal(0, 0.05, 1)
                
    if c[0]=='l' and c[1]=='l':
            e = lambda v,z,: np.sum(np.abs(z-z[9,9]*cauchy(19,v)),1)
            a=float(1/float(a))
            b=float(1/float(b))
            vi=[a,b]    
            v, success = leastsq(e, vi, args=(z), maxfev=1000)
            #print vi
            if np.abs(float(v[0]-float(1/float(a))))>0.5 or v[0]<0.08 or v[0]>4:
                v[0]=a+np.random.normal(0, 0.05, 1)
                #print float(1/float(a))
            v[0]=1/v[0]

            if  np.abs(float(v[1]-float(1/float(b))))>0.5 or v[1]<0.08 or v[1]>4:
                v[1]=b+np.random.normal(0, 0.05, 1)
            v[1]=1/v[1]
    if c[0]=='g' and c[1]=='l':
            e = lambda v,z,: np.sum(np.abs(z-z[9,9]*exponcauchy(19,v)),1)
            b=float(1/float(b)) 
            vi=[a,b]   
 
            v, success = leastsq(e, vi, args=(z), maxfev=1000)
            #print 'ham',v
            #print vi
            if np.abs(float(v[0]-a))>1 or v[0]<0.5 or v[0]>6:  
                v[0]=a+np.random.normal(0, 0.05, 1)
                #print float(1/float(a))
            if  np.abs(float(v[1]-float(b)))>0.5 or v[1]<0.08 or v[1]>4:
                v[1]=b+np.random.normal(0, 0.05, 1)
            v[1]=1/v[1]
    if c[0]=='l' and c[1]=='g':
            e = lambda v,z,: np.sum(np.abs(z-z[9,9]*cauchyexpon(19,v)),1)
            a=float(1/float(a)) 
            vi=[a,b]    
            v, success = leastsq(e, vi, args=(z), maxfev=1000)
            #print vi
            if np.abs(float(v[1]-b))>1 or v[0]<0.5 or v[0]>6:  
                v[1]=b+np.random.normal(0, 0.05, 1)
                #print float(1/float(a))
            if  np.abs(float(v[0]-float(1/float(a))))>0.5 or v[1]<0.08 or v[1]>4:
                v[0]=a+np.random.normal(0, 0.05, 1)
            v[0]=1/v[0]
    #print c 
    #print vi
    #print v                                       
    return v,c

def exp_hand(z,newn):
    c=0
    for i in range(len(z)):
        if z[i]==newn:
            c=1
            break
    return c
def exp_yn(amp,ampref,test):
    artest=np.array(test)
    ol=np.nonzero(artest==1)
    #print np.size(ol)
    #print len(amp)/2
    #print amp
    if np.size(ol)>len(amp)/2:
        ver=0
    else:
        o=np.nonzero(artest==0)
        vamp=np.array(amp)
        vampref=np.array(ampref)
        ivamref=vampref[o]
        ivam=vamp[o]
        ii=ivamref[ivamref==np.amax(ivamref)]
        jj=ivam[ivamref==np.amax(ivamref)]
        ver=1
        if np.size(ol)>1:
            ol=ol[0]
            #print np.size(ol)
        #print vampref
        #print vamp
        #print ol
        if np.size(ol)>1:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj))) 
                if (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))>50 or (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))<0.0200:
                    #print 'lela'
                    ver=0
        else:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0]) <0.001
                if (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])>50 or (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])<0.0200:
                    #print 'lela'
                    ver=0
            
    return ver
def exp_ync(amp,ampref,test):
    artest=np.array(test)
    ol=np.nonzero(artest==1)
    #print np.size(ol)
    #print len(amp)/2
    if np.size(ol)>len(amp)/2:
        ver=0
    else:
        o=np.nonzero(artest==0)
        vamp=np.array(amp)
        vampref=np.array(ampref)
        ivamref=vampref[o]
        ivam=vamp[o]
        ii=ivamref[ivamref==np.amax(ivamref)]
        jj=ivam[ivamref==np.amax(ivamref)]
        ver=1
        if np.size(ol)>1:
            ol=ol[0]
            #print np.size(ol)
        #print vampref
        #print vamp
        #print ol
        if np.size(ol)>1:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj))) 
                if (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))>10 or (((vampref[ol[kkk]]/ii)/(vamp[ol[kkk]]/jj)))<0.065:
                    #print 'lela'
                    ver=0
        else:
            for kkk in range(np.size(ol)):
                #print float(((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0]) 
                if (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])>10 or (((vampref[ol[kkk]][0]/ii)/(vamp[ol[kkk]][0]/jj))[0])<0.065:
                    #print 'lela'
                    ver=0
            
    return ver
        
def ser_mat(tabp):    
    boc=-1
    inp=[]
    inm=[]
    while boc<len(tabp)-1:
        boc+=1
        a=str(tabp[boc][2])
        newn=''
        for j in range(len(a)):
            if (a[j]=='_')==1:
                break
            newn+=a[j]
        
        r1=[]
        r1.append((boc))
        for jj in range(boc+1,len(tabp)):
            nomn=str(tabp[jj][2])
            #print nomn[0:j]
            try:
                if nomn[0:j+1]==newn+'_':
                    r1.append((jj))
                    #print 'ok'
                else:
                    break
            except:
                    break
        boc=boc+len(r1)-1
        inm.append(newn)
        inp.append(r1)
    return inm,inp

def tab2tab(X):
    Z=[]
    for i in range(len(X)):
        Z.append([X[i]])
    Z=np.array(Z)
    return Z

def norm(a,b,c):
    A=(float(a)/float(a+b))*float(c)
    B=(float(b)/float(a+b))*float(c)
    return A,B
    
    
