import random, copy, math,sys
import numpy as np
import fonc_util as fon

def k_means(t, nbclusters=2, nbiter=3, medoids=False, soft=False, beta=200.0,\
        distance=lambda x,y: math.sqrt(np.dot(x-y,(x-y).conj())),\
        responsability=lambda beta,d: math.exp(-1 * beta * d)):

    nbobs = t.shape[0]
    nbfeatures = t.shape[1]
    # find xranges for each features
    min_max = []
    for f in xrange(nbfeatures):
        min_max.append((t[:,f].min(), t[:,f].max()))

    ### Soft => Normalization, otherwise "beta" has no meaning!
    if soft:
        for f in xrange(nbfeatures):
            t[:,f] -= min_max[f][0]
            t[:,f] /= (min_max[f][1]-min_max[f][0])
    min_max = []
    for f in xrange(nbfeatures):
        min_max.append((t[:,f].min(), t[:,f].max()))
    ### /Normalization # ugly

    result = {}
    quality = 0.0 # sum of the means of the distances to centroids
    random.seed()
    tmpdist = np.ndarray([nbobs,nbclusters], np.float64) # distance obs<->clust
    tmpresp = np.ndarray([nbobs,nbclusters], np.float64) # responsability o<->c
    # iterate for the best quality
    for iteration in xrange(nbiter):
        clusters = [[] for c in xrange(nbclusters)]
        # Step 1: place nbclusters seeds for each features
        centroids = [np.array([random.uniform(min_max[f][0], min_max[f][1])\
                for f in xrange(nbfeatures)], np.float64)\
                for c in xrange(nbclusters)]
        old_centroids = [np.array([-1 for f in xrange(nbfeatures)], np.float64)\
                for c in xrange(nbclusters)] # should not be init, TODO
        new_sum = math.fsum([distance(centroids[c], old_centroids[c])\
                for c in xrange(nbclusters)])
        old_sum = sys.maxint
        np.seterr(invalid='raise')
        # iterate until convergence
        while new_sum < old_sum :
            old_centroids = copy.deepcopy(centroids)
            old_sum = new_sum
            for c in xrange(nbclusters):
                clusters[c] = []
            # precompute distance to all centroids/medoids for all observations
            for c in xrange(nbclusters):
                for o in xrange(nbobs):
                    tmpdist[o,c] = distance(centroids[c], t[o,:])
            if soft:
                # Step 2: compute the degree of assignment for each object
                for o in xrange(nbobs):
                    for c in xrange(nbclusters):
                        tmpresp[o,c] = responsability(beta, tmpdist[o,c])
                for o in xrange(nbobs):
                    tmpresp[o,:] /= math.fsum(tmpresp[o,:])
            else:
                # Step 2: assign each object to the closest centroid
                for o in xrange(nbobs):
                    clusters[tmpdist[o,:].argmin()].append(o)
            # Step 3: recalculate the positions of the nbclusters centroids
            for c in xrange(nbclusters):
                if medoids:
                    if soft:
                        print "ERROR: Soft medoids not implemented"
                        sys.exit(-1)
                    else:
                        tmpmin = sys.maxint
                        argmin = 0
                        for o in clusters[c]:
                            if tmpdist[o,c] < tmpmin:
                                tmpmin = tmpdist[o,c]
                                argmin = o
                        centroids[c] = t[argmin,:]
                else:
                    mean = np.array([0 for i in xrange(nbfeatures)], np.float64)
                    if soft:
                        for o in xrange(nbobs):
                            mean += tmpresp[o,c] * t[o,:]
                        mean /= math.fsum(tmpresp[:,c])
                    else:
                        for o in clusters[c]:
                            mean += t[o,:]
                        mean = map(lambda x: x/len(clusters[c]), mean)
                    centroids[c] = np.array(mean, np.float64)
            #print centroids
            new_sum = math.fsum([distance(centroids[c], old_centroids[c])\
                    for c in xrange(nbclusters)])
            #print "(k-means) old and new sum: ", old_sum, new_sum
        if soft:
            for o in xrange(nbobs):
                clusters[tmpdist[o,:].argmin()].append(o)
        quality = math.fsum([math.fsum([tmpdist[o][c] for o in clusters[c]])\
                /(len(clusters[c])+1) for c in xrange(nbclusters)])
        if not quality in result or quality > result['quality']:
            result['quality'] = quality
            result['centroids'] = centroids
            result['clusters'] = clusters
    return result
def flou(t,a,b,c):

    if t<=a:
        nt=0
    else:
        if t<=b :
            nt=(t-a)**2/((b-a)*(c-a))
        else:
            if t<=c:
                nt=1-((t-c)**2/((c-b)*(c-a)))
            else:
                nt=1
    
    
    return nt
    
    
def msatkin(tab):  
   # print tab
    ntab=np.array(tab)
    #print ntab
    ki5a=[]
    baya=[]
    a3lamb=[]
    a3lamk=[]
    yezzi=[]
    
    for h in range(len(ntab)):
        for h1 in range(len(ntab)):
            if h!=h1 and len(np.nonzero(np.array(yezzi)==h)[0])<1:
                if np.abs(float(ntab[h1,0])-float(ntab[h,0]))<0.0001 and np.abs(float(ntab[h1,1])-float(ntab[h,1]))<0.00001 :
                    #print h,yezzi,len(np.nonzero(np.array(yezzi)==h)[0])
                    yezzi.append(h1)
                    cond=0
                    if np.abs(float(ntab[h,11]))<3 and np.abs(float(ntab[h,12]))<5:
                        if np.abs(float(ntab[h1,11]))<3 and np.abs(float(ntab[h1,12]))<5:
                            cond=1
                            a3lamb.append(h1)
                            a3lamb.append(h)
                            t,tt=fon.norm(float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]))
                            ntab[h1,3]=str(t)
                            ntab[h,3]=str(tt)
                            #print t,tt,float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]),h,h1

                        else:
                            ki5a.append(h1)
                            baya.append(h)
                            t,tt=fon.norm(float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]))
                            ntab[h1,3]=str(t)
                            ntab[h,3]=str(tt)
                            #print t,tt,float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]),h,h1
                    else:
                        if np.abs(float(ntab[h1,11]))<3 and np.abs(float(ntab[h1,12]))<5:
                            ki5a.append(h)
                            baya.append(h1)
                            t,tt=fon.norm(float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]))
                            ntab[h1,3]=str(t)
                            ntab[h,3]=str(tt)
                            #print t,tt,float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]),h,h1
                        else:
                            cond=2
                            a3lamk.append(h)
                            a3lamk.append(h1)
                            t,tt=fon.norm(float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]))
                            ntab[h1,3]=str(t)
                            ntab[h,3]=str(tt)
                            #print t,tt,float(ntab[h1,4]),float(ntab[h,4]),float(ntab[h1,3]),h,h1


                    
                    #print ntab[h1][2],ntab[h1][12],ntab[h1][12]
                    #print ntab[h][2],ntab[h][12],ntab[h][12]
                    #break
    #print ki5a, ntab[ki5a,2],ntab[ki5a,11],ntab[ki5a,12]
    #print '%%%%%%%%%%'
    #print baya, ntab[baya,2],ntab[baya,11],ntab[baya,12]
    #print '%%%%%%%%%%'
    #if len(a3lamb)>0:
        #print a3lamb, ntab[a3lamb,2],ntab[a3lamb,11],ntab[a3lamb,12]
    #else:
       # print 'que dalle'      
    #print '%%%%%%%%%%'
    #if len(a3lamk)>0:
        #print a3lamk, ntab[a3lamk,2],ntab[a3lamk,11],ntab[a3lamk,12]
   # else:
        #print 'que dalle'
    
    return ntab
