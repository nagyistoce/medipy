"""
Peak and metabolite finding
"""
import numpy as np
import fonc_util as f
def alig(ic,jc,data2,amplitref,wa7da):
        amplitref=float(amplitref)*1.00000000000000000000000001
        peak=[]
        amp=[]
        pos=[]
        #r1.append((boc))
        #peak=np.array([[0,0]])
        #amp= np.array([0])
        #pos=np.array([[0,0]])
        indi=0
	
        z=data2[ic-3:ic+4,jc-5:jc+6]
	#print z
        m=np.size(z[:,1])
        n=np.size(z[1,:])
        for i in range(1,m-2):
            for j in range(1,n-1):
                zi=z[i-1:i+2,j-1:j+2]
                im=f.find(zi,zi==np.amax(zi))
                if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>wa7da:
                    #print i,j
                    #t1=np.array([[3-i,5-j]])*1.
                    peak.append((3-i,5-j))
                    #peak= np.concatenate((peak, t1), axis=0)
                    #t2=np.array([np.amax(zi)])*1.
                    #amp= np.concatenate((amp, t2), axis=0)
                    amp.append((np.amax(zi)))
        #print peak
        #raw_input()
        if len(peak)>0:
                
                #peak=np.delete(peak,0,0)
                #amp=np.delete(amp,0)
                if len(peak)==1:
                    #print 'd5al21'
                    pos=peak
                    pos=mod(pos[0])
                     #print 'oooook'
                else:
                    
                    avi=np.abs(peak)
                    #print avi
                    av=f.find(avi[:,0],avi[:,0]==np.amin(avi[:,0]))
                    #print av,'oooo',np.abs(peak)
                    av=av[0]
                    #print av,'ooooooooooooo'
                    if np.size(av)==1:
                        #print 'd5al22'
                        #print int(av[0][0])
                        pos=peak[int(av[0])]
                        pos=mod(pos)
                        #print 'oooook'
                    else:
                        
                        #print 'd5al'
                        #print av
                        iv=mod(av)
                        pook=aff(peak,iv)
                        #print pook,'kkkk'
                        bvi=np.abs(pook)
                        bv=f.find(bvi[:,1],bvi[:,1]==np.amin(bvi[:,1]))
                        #if np.size(bv)==1:
                            #print 'd5al23'
                        pos=pook[int(bv[0][0])]
                            #print pos
                        pos=mod(pos)
                        #else:
                            #print 'd5al24'
                            #print bv[0]
                            #bv=int(bv[0][0])
                            #print bv,pook[np.array(bv)[0],:]
                            #piik=pook[bv,:][0]
                            #amp1=amp[av]
                            #amp2=amp1[bv]
                            #dif=np.abs(amp[bv]-amplitref)
                            #difa=f.find(dif,dif==np.amin(dif))
                            #pos= piik[difa,:][0]
                            #pos=mod(pos)
        else:
            amplitref=float(amplitref)*1.00000000000000000000000001
            peak=[]
            amp=[]
            pos=[]
            #r1.append((boc))
            #peak=np.array([[0,0]])
            #amp= np.array([0])
            #pos=np.array([[0,0]])
            #print 'ena hna'
            indi=0
            z=data2[ic-6:ic+7,jc-6:jc+7]
            m=np.size(z[:,1])
            n=np.size(z[1,:])
            for i in range(1,m-2):
                for j in range(1,n-1):
                    zi=z[i-1:i+2,j-1:j+2]
                    im=f.find(zi,zi==np.amax(zi))
                    if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>wa7da:
                        #print i,j
                        #t1=np.array([[3-i,5-j]])*1.
                        peak.append((6-i,6-j))
                        #print 'ouuuf'
                        #peak= np.concatenate((peak, t1), axis=0)
                        #t2=np.array([np.amax(zi)])*1.
                        #amp= np.concatenate((amp, t2), axis=0)
                        amp.append((np.amax(zi)))
            #print peak
        #raw_input()
            if len(peak)>0:
                    #print 'd5al'
                    #peak=np.delete(peak,0,0)
                    #amp=np.delete(amp,0)
                    if len(peak)==1:
                        #print 'd5al1'
                        pos=peak
                        pos=mod(pos[0])
                        #print 'oooook'
                    else:
                        avi=np.abs(peak)
                        #print avi
                        av=f.find(avi[:,0],avi[:,0]==np.amin(avi[:,0]))
                        #print av,'oooo',np.abs(peak)
                        av=av[0]
                        #print av,'ooooooooooooo'
                        if np.size(av)==1:
                            #print 'd5al2'
                            #print int(av[0][0])
                            #print av,'hhhhhhhhhhhhhh'
                            pos=peak[int(av[0])]
                            pos=mod(pos)
                            #print 'oooook'
                        else:
                            
                            #print 'd5al'
                            #print av
                            iv=mod(av)
                            pook=aff(peak,iv)
                            #print pook,'kkkk'
                            bvi=np.abs(pook)
                            bv=f.find(bvi[:,1],bvi[:,1]==np.amin(bvi[:,1]))
                            if np.size(bv)==1:
                                #print 'd5al3'
                                pos=pook[int(bv[0][0])]
                                pos=mod(pos)
                            else:
                                #print 'd5al4'
                                bv=bv[0]
                                piik=aff(pook,bv)
                                #amp1=amp[av]
                                #amp2=amp1[bv]
                                cd=deff(aff(amp,bv),amplitref)
                                dif=np.abs(cd)
                                difa=f.find(dif,dif==np.amin(dif))
                                difa=difa[0]
                                #print difa,'jjjjjjjjjjjjjj'
                                pos= piik[difa[0]]
                                pos=mod(pos)
                    #pos=pos[0][:]         
            else:
                indi=1
                tt=data2[ic-1:ic+2,jc-2:jc+3]*1.
                #print tt
                if len(f.find(tt,tt>wa7da)[0]) >11:
                    #print 'im in first'
                    im=f.find(tt,tt==np.amax(tt))
                    #print im
                    pos.append((1-im[0][0],2-im[1][0]))
                    pos=mod(pos[0])
                else:
                        #print 'im in'
                        pos=np.array([100,100])
                        indi=1
                        
                        nz=data2[ic-2:ic+3,jc-3:jc+4]
                        
                        val=amplitref
                        nm=np.size(nz[:,1])
                        nn=np.size(nz[1,:])
                        zios=np.zeros([nm,nn])
                        for ni in range(1,nm-2):
                            for nj in range(1,nn-2):
                                nzi=nz[ni-1:ni+2,nj-1:nj+2]
                                #print nzi[1,1]
                                if (float(nzi[1,1])>wa7da):
                                        difo=f.find(nzi,nzi[1,1]>nzi)
                                        zios[ni,nj]=np.size(difo[1])
                        #print zios
                        #raw_input()
                       
                        if np.amax(zios)>4:
                            #print nz
                            #print zios
                            hv=f.find(zios,zios>4)
                            #print hv
                            #print hv[0],hv[1]
                            crit=np.abs(nz-val)
                            #print crit
                            pos=[hv[0][0],hv[1][0]]
                            for ll in range(1,len(hv[1])-1):
                                #print 'd5al!!!'
                                if crit[hv[0][ll],hv[1][ll]]<crit[pos[0],pos[1]]:
                                    pos=[hv[0][ll],hv[1][ll]]
                            pos=np.array([2,3])-pos
                            #print pos
                            #nic=np.array([ic,jc])-pos
                            #tt=data2[nic[0]-1:nic[0]+2,nic[1]-1:nic[1]+2]
                            #tt=nz[pos[0]-1:pos[0]+2,pos[1]-2:pos[1]+2]
                            #print tt
                            #im=f.find(tt,tt==np.amax(tt))
                            #pos[0]=pos[0]+(1-im[0][0])
                            #pos[1]=pos[1]+(1-im[1][0])
                            #tt=nz[pos[0]-1:pos[0]+2,pos[1]-1:pos[1]+2]
                            #im=f.find(tt,tt==np.amax(tt))
                            #pos[0]=pos[0]+im[0][0]-1
                            #pos[1]=pos[1]+im[1][0]-1
                            #print nz
                            #pos=np.array([2,4])-pos
                            #print pos[1],'hhhhhhhhhhhhhhhhh'
                        pos=[pos[0],pos[1]]
        #print pos
        if np.abs(pos[0])==4:
            pos=np.array([100,100])
            pos=[pos[0],pos[1]]
        if pos[0]==100 :
                pos=[]
                indi=1
                tt=data2[ic-1:ic+2,jc-2:jc+3]*1.
                #print 'one more time'
                #print tt
                if len(f.find(tt,tt>wa7da)[0]) >11:
                    im=f.find(tt,tt==np.amax(tt))
                    #print im
                    pos.append((1-im[0][0],2-im[1][0]))
                    pos=mod(pos[0])
                    #print pos
                else:
                    pos=[100,100]
        return pos,indi
    

def aligrec(ic,jc,data2,amplitref):
        amplitref=float(amplitref)*1.00000000000000000000000001
        peak=[]
        amp=[]
        pos=[]
        #r1.append((boc))
        #peak=np.array([[0,0]])
        #amp= np.array([0])
        #pos=np.array([[0,0]])
        indi=0
    
        z=data2[ic-19:ic+20,jc-11:jc+12]
    #print z
        m=np.size(z[:,1])
        n=np.size(z[1,:])
        for i in range(1,m-2):
            for j in range(1,n-1):
                zi=z[i-1:i+2,j-1:j+2]
                im=f.find(zi,zi==np.amax(zi))
                if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>300000:
                    #print i,j
                    #t1=np.array([[3-i,5-j]])*1.
                    peak.append((19-i,11-j))
                    #peak= np.concatenate((peak, t1), axis=0)
                    #t2=np.array([np.amax(zi)])*1.
                    #amp= np.concatenate((amp, t2), axis=0)
                    amp.append((np.amax(zi)))
        #print peak
        #raw_input()
        if len(peak)>0:
                
                #peak=np.delete(peak,0,0)
                #amp=np.delete(amp,0)
                if len(peak)==1:
                    #print 'd5al21'
                    pos=peak
                    pos=mod(pos[0])
                     #print 'oooook'
                else:
                    
                    avi=np.abs(peak)
                    #print avi
                    av=f.find(avi[:,0],avi[:,0]==np.amin(avi[:,0]))
                    #print av,'oooo',np.abs(peak)
                    av=av[0]
                    #print av,'ooooooooooooo'
                    if np.size(av)==1:
                        #print 'd5al22'
                        #print int(av[0][0])
                        pos=peak[int(av[0])]
                        pos=mod(pos)
                        #print 'oooook'
                    else:
                        
                        #print 'd5al'
                        #print av
                        iv=mod(av)
                        pook=aff(peak,iv)
                        #print pook,'kkkk'
                        bvi=np.abs(pook)
                        bv=f.find(bvi[:,1],bvi[:,1]==np.amin(bvi[:,1]))
                        #if np.size(bv)==1:
                            #print 'd5al23'
                        pos=pook[int(bv[0][0])]
                            #print pos
                        pos=mod(pos)
                        #else:
                            #print 'd5al24'
                            #print bv[0]
                            #bv=int(bv[0][0])
                            #print bv,pook[np.array(bv)[0],:]
                            #piik=pook[bv,:][0]
                            #amp1=amp[av]
                            #amp2=amp1[bv]
                            #dif=np.abs(amp[bv]-amplitref)
                            #difa=f.find(dif,dif==np.amin(dif))
                            #pos= piik[difa,:][0]
                            #pos=mod(pos)
        else:
            amplitref=float(amplitref)*1.00000000000000000000000001
            peak=[]
            amp=[]
            pos=[]
            #r1.append((boc))
            #peak=np.array([[0,0]])
            #amp= np.array([0])
            #pos=np.array([[0,0]])
            #print 'ena hna'
            indi=0
            z=data2[ic-6:ic+7,jc-6:jc+7]
            m=np.size(z[:,1])
            n=np.size(z[1,:])
            for i in range(1,m-2):
                for j in range(1,n-1):
                    zi=z[i-1:i+2,j-1:j+2]
                    im=f.find(zi,zi==np.amax(zi))
                    if im[0][0]==1 and im[1][0]==1 and zi[im[0][0],im[1][0]]>190000:
                        #print i,j
                        #t1=np.array([[3-i,5-j]])*1.
                        peak.append((6-i,6-j))
                        #print 'ouuuf'
                        #peak= np.concatenate((peak, t1), axis=0)
                        #t2=np.array([np.amax(zi)])*1.
                        #amp= np.concatenate((amp, t2), axis=0)
                        amp.append((np.amax(zi)))
            #print peak
        #raw_input()
            if len(peak)>0:
                    #print 'd5al'
                    #peak=np.delete(peak,0,0)
                    #amp=np.delete(amp,0)
                    if len(peak)==1:
                        #print 'd5al1'
                        pos=peak
                        pos=mod(pos[0])
                        #print 'oooook'
                    else:
                        avi=np.abs(peak)
                        #print avi
                        av=f.find(avi[:,0],avi[:,0]==np.amin(avi[:,0]))
                        #print av,'oooo',np.abs(peak)
                        av=av[0]
                        #print av,'ooooooooooooo'
                        if np.size(av)==1:
                            #print 'd5al2'
                            #print int(av[0][0])
                            #print av,'hhhhhhhhhhhhhh'
                            pos=peak[int(av[0])]
                            pos=mod(pos)
                            #print 'oooook'
                        else:
                            
                            #print 'd5al'
                            #print av
                            iv=mod(av)
                            pook=aff(peak,iv)
                            #print pook,'kkkk'
                            bvi=np.abs(pook)
                            bv=f.find(bvi[:,1],bvi[:,1]==np.amin(bvi[:,1]))
                            if np.size(bv)==1:
                                #print 'd5al3'
                                pos=pook[int(bv[0][0])]
                                pos=mod(pos)
                            else:
                                #print 'd5al4'
                                bv=bv[0]
                                piik=aff(pook,bv)
                                #amp1=amp[av]
                                #amp2=amp1[bv]
                                cd=deff(aff(amp,bv),amplitref)
                                dif=np.abs(cd)
                                difa=f.find(dif,dif==np.amin(dif))
                                difa=difa[0]
                                #print difa,'jjjjjjjjjjjjjj'
                                pos= piik[difa[0]]
                                pos=mod(pos)
                    #pos=pos[0][:]         
            else:
                indi=1
                tt=data2[ic-1:ic+2,jc-2:jc+3]*1.
                #print tt
                if len(f.find(tt,tt>190000)[0]) >11:
                    #print 'im in first'
                    im=f.find(tt,tt==np.amax(tt))
                    #print im
                    pos.append((1-im[0][0],2-im[1][0]))
                    pos=mod(pos[0])
                else:
                        #print 'im in'
                        pos=np.array([100,100])
                        indi=1
                        
                        nz=data2[ic-2:ic+3,jc-3:jc+4]
                        
                        val=amplitref
                        nm=np.size(nz[:,1])
                        nn=np.size(nz[1,:])
                        zios=np.zeros([nm,nn])
                        for ni in range(1,nm-2):
                            for nj in range(1,nn-2):
                                nzi=nz[ni-1:ni+2,nj-1:nj+2]
                                #print nzi[1,1]
                                if (float(nzi[1,1])>190000):
                                        difo=f.find(nzi,nzi[1,1]>nzi)
                                        zios[ni,nj]=np.size(difo[1])
                        #print zios
                        #raw_input()
                       
                        if np.amax(zios)>4:
                            #print nz
                            #print zios
                            hv=f.find(zios,zios>4)
                            #print hv
                            #print hv[0],hv[1]
                            crit=np.abs(nz-val)
                            #print crit
                            pos=[hv[0][0],hv[1][0]]
                            for ll in range(1,len(hv[1])-1):
                                #print 'd5al!!!'
                                if crit[hv[0][ll],hv[1][ll]]<crit[pos[0],pos[1]]:
                                    pos=[hv[0][ll],hv[1][ll]]
                            pos=np.array([2,3])-pos
                            #print pos
                            #nic=np.array([ic,jc])-pos
                            #tt=data2[nic[0]-1:nic[0]+2,nic[1]-1:nic[1]+2]
                            #tt=nz[pos[0]-1:pos[0]+2,pos[1]-2:pos[1]+2]
                            #print tt
                            #im=f.find(tt,tt==np.amax(tt))
                            #pos[0]=pos[0]+(1-im[0][0])
                            #pos[1]=pos[1]+(1-im[1][0])
                            #tt=nz[pos[0]-1:pos[0]+2,pos[1]-1:pos[1]+2]
                            #im=f.find(tt,tt==np.amax(tt))
                            #pos[0]=pos[0]+im[0][0]-1
                            #pos[1]=pos[1]+im[1][0]-1
                            #print nz
                            #pos=np.array([2,4])-pos
                            #print pos[1],'hhhhhhhhhhhhhhhhh'
                        pos=[pos[0],pos[1]]
        #print pos
        return pos,indi
def aff(p,ind):
    np=[]
    for i in range(len(ind)):
        np.append((p[ind[i]]))
    return np
def mod(ind):
    iv=[]
    for kkk in range(np.size(ind)):
        iv.append((ind[kkk]))
    return iv
def deff(vect,val):
    np=[]
    for i in range(len(vect)):
        np.append((vect[i]-val))
    return np
