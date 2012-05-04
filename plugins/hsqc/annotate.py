"""
Metabolite annotation
"""
import corps
import logging
import sys
import warnings
import numpy as np
from medipy.io import load
import os
import fonc_util as f
#import segmentation as seg
from medipy.base import Image
import gen_in
import lawij
import ser_xml
import MCMC as m
import xml_gen
import locale

def ann(root,wa7da):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    #root='/home/miv/belghith/Bureau/KAROM/Akram/nmr/RP/5/pdata/1'
    tabp,tabc=gen_in.ini1(root)
    
    boc=-1
    lolo= load(os.path.join(root,"2rr"))
    H=lolo.data
    D=H[0,:,:]
    H=D
    #print corps.interquartile_range(H)
    #print np.size(H,0),np.size(H,1)
    list_exp=ser_xml.exceptions(root)
    raxep=[]
    while boc<len(tabp)-1:
        boc+=1
        a=str(tabp[boc][2])
        #print a,boc
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
        #print tabp[r1[0]:r1[len(r1)-1]+1]
        #raw_input()
        #print len(r1)
        nt=tabp[r1[0]:r1[len(r1)-1]+1]
        #print nt
        #print 'Start'
        test=[]
        testc3=[]
        con=0
        ampref=[]
        amp=[]
        newrt=[]
        #print newn
        ham=0
        ed5il=0
        for jj in range(len(nt)):
    
            #print newn
            #print jj
            #print nt[jj][0],nt[jj][1]
            #print nt[jj][0],nt[jj][1]
            r,indi=lawij.alig(nt[jj][0],nt[jj][1],H,nt[jj][3],wa7da)
            
            #print r,indi,nt[jj][2],nt[jj][0]-r[0],nt[jj][1]-r[1],nt[jj][7]
            
            if nt[jj][7]=='n':
                ed5il=1
            
            #raw_input()
            if r[0]==100 and nt[jj][7]=='y':
                #print r,indi,nt[jj][2],nt[jj][0]-r[0],nt[jj][1]-r[1],nt[jj][7]
                ham=1
                break
            
            if indi==0 :
                con=con+1
                
                if np.abs(r[0])==4 or np.abs(r[1])==6:
                    testc3.append((1))
                    
                else:
                    testc3.append((0))
                test.append((0))
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                zayneb=H[(nt[jj][0]-r[0])-3:(nt[jj][0]-r[0])+4,(nt[jj][1]-r[1])-3:(nt[jj][1]-r[1])+4]*1.
                nr=f.subpix2(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                zayneb=H[(nt[jj][0]-r[0])-9:(nt[jj][0]-r[0])+10,(nt[jj][1]-r[1])-9:(nt[jj][1]-r[1])+10]*1.
                #nr=f.subpix(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                chl=f.dephcl(zayneb)
                chg=f.dephcg(zayneb)
                #ch,congl=f.dephc(zayneb)
                ch,congl=f.dephcaprio(zayneb,float(nt[jj][4]),float(nt[jj][5]),nt[jj][6])
                #print congl
                #cax = ax.imshow(zayneb, interpolation='nearest')
                #plt.show()
                #print ch
                #plt.show()
                #print nt[jj][3]
                ampref.append(float(nt[jj][3])*1.)
                amp.append( H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.)
                #print str(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])])
                newrt.append((nr[0],nr[1],nt[jj][2],str(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]),float(nt[jj][3]),chg[1],chl[0],chl[1],ch[0],ch[1],congl,r[0],r[1]))
                #print newrt
            else:
                if r[0]<100 :
                    if np.abs(r[0])==3 or np.abs(r[1])==5:
                        testc3.append((1))
                        
                    else:
                        testc3.append((0))
                    con=con+1
                    test.append((1))
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111)
                    zayneb=H[(nt[jj][0]-r[0])-3:(nt[jj][0]-r[0])+4,(nt[jj][1]-r[1])-3:(nt[jj][1]-r[1])+4]*1.
                    nr=f.subpix2(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                    #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                    zayneb=H[(nt[jj][0]-r[0])-9:(nt[jj][0]-r[0])+10,(nt[jj][1]-r[1])-9:(nt[jj][1]-r[1])+10]*1.
                    #nr=f.subpix(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                    #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                    chl=f.dephcl(zayneb)
                    chg=f.dephcg(zayneb)
                    #ch,congl=f.dephc(zayneb)
                    ch,congl=f.dephcaprio(zayneb,float(nt[jj][4]),float(nt[jj][5]),nt[jj][6])
                    ampref.append(float(nt[jj][3])*1.)
                    amp.append(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.)
                    newrt.append((nr[0],nr[1],nt[jj][2],H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.,float(nt[jj][3]),chg[1],chl[0],chl[1],ch[0],ch[1],congl,r[0],r[1]))
                else:
                    test.append((2))
                    testc3.append((2))
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111)
                    zayneb=H[nt[jj][0]-7:(nt[jj][0])+8,(nt[jj][1])-7:(nt[jj][1])+8]*1.
                    ampref.append(float(nt[jj][3])*1.)
                    amp.append(0)
                    #cax = ax.imshow(zayneb, interpolation='nearest')
                    #plt.show()
                    #raw_input()
                    newrt.append((0,0,0,0,0,0,0,0,0,0,0,0,0))
                    

        #raw_input()
        #print newn
        #print ampref
        #print amp,newn,testc3
        #print nt
        #print newrt
        #raw_input()
        #amptest=np.nonzero(amp>0)
        o=np.nonzero(testc3==0)
        vamp=np.array(amp)
        ivamref=vamp[o]
        o=np.nonzero(ivamref>210000)
        #print newn,'test'
        #print 'ham',ham
        #if (float(len(o[0]))*1.000000001)/float(len(ivamref)*1.00000001)>0.4 or f.exp_hand(list_exp,newn)==1 or ham==0:
        if ham==0:
            #print 'd5aal'
            if len(nt)==con:
                if len(nt)==1:
                    raxep.append(newrt[0])
                    #print 'accepted'
                else:
                    artestc3=np.array(testc3)
                    olc3=np.nonzero(artestc3==1)
                    if np.size(olc3)>0:
                        if f.exp_ync(amp,ampref,testc3)==1:
                            #print 'ouffffff'
                            artest=np.array(test)
                            ol=np.nonzero(artest==1)
                            o=np.nonzero(artest==0)
                            if np.size(ol)>0:
                                #print 'accepted with some conditions'
                                #print f.exp_ync(amp,ampref,testc3)
                                #if f.exp_ync(amp,ampref,testc3)==0:
                                    #print 'llllllaaaaaaaaa'
                                if f.exp_yn(amp,ampref,test)==1:
                                    for no in range(len(newrt)):
                                        raxep.append(newrt[no])
                                elif f.exp_hand(list_exp,newn)==1:
                                    artest=np.array(test)
                                    rnt=np.array(newrt)
                                    ol=np.nonzero(artest==1)
                                    o=np.nonzero(artest==0)
                                    vo=rnt[o]
                                    #raw_input()
                                    for no in range(len(vo)):
                                        #print '%%%%%%%%%'
                                        zi=lawij.mod(vo[no])
                                        #print zi
                                        #print artest
                                        #raw_input()
                                        raxep.append((zi))
                            else:
                                #print f.exp_ync(amp,ampref,testc3)
                                for no in range(len(nt)):
                                    raxep.append(newrt[no])
                                #print 'accepted'                                    
                        else:
                            #print 'ouuuuut'
                            #print f.exp_hand(list_exp,newn)
                            if f.exp_hand(list_exp,newn)==1 :
                                artest=np.array(test)
                                artestc3=np.array(testc3)
                                rnt=np.array(newrt)
                                #print nt
                                #raw_input()
                                ol=np.nonzero(artest==1)
                                condlist=[artest==0,artestc3==0]
                                g=condlist[0]*condlist[1]
                                #print g,artest,artestc3
                                o=np.nonzero(artest==0)
                                #print '%%%hhh%%%%%%'
                                #print rnt
                                #print test
                                vo=rnt[g]
                                #print vo,'%%%hhh%%%%%%'
                                #raw_input()
                                for no in range(len(vo)):
                                            #print '%%%%%%%%%'
                                            zi=lawij.mod(vo[no])
                                            #print zi
                                            #print artest
                                            #raw_input()
                                            raxep.append((zi))
                    else:
                                #print f.exp_ync(amp,ampref,testc3)
                                for no in range(len(nt)):
                                    raxep.append(newrt[no])
                                #print 'accepted'
    
            else:
                if ham==0:
                    if f.exp_hand(list_exp,newn)==1 or ed5il==1:
                        artest=np.array(test)
                        #print newrt
                        rnt=xml_gen.cell2tab8(newrt)
                        #print nt
                        #raw_input()
                        ol=np.nonzero(artest==1)
                        o=np.nonzero(artest==0)
                        #print '%%%hhh%%%%%%'
                        rnt=np.array(rnt)
                        #print test
                        vo=rnt[o]
                        #print vo
                        #print len(vo)
                        #raw_input()
                        for no in range(len(vo)):
                                    #print '%%%%%%%%%'
                                    zi=lawij.mod(vo[no])
                                    #print zi
                                    #print artest
                                    #raw_input()
                                    raxep.append((zi))
                    #print 'may be...'
                #else:
                    #print 'refused without discussion'    
            #print test
            
                
            
            #print 'DONE'
            #print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            #print '        '
            #raw_input()
    locale.setlocale(locale.LC_ALL, old_locale)
    return raxep,tabc
    #for lll in range(len(raxep)):
        #print  raxep[lll]  




def ann_old(root,wa7da):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    #root='/home/miv/belghith/Bureau/KAROM/Akram/nmr/RP/5/pdata/1'
    tabp,tabc=gen_in.ini1(root)
    
    boc=-1
    lolo= load(os.path.join(root,"2rr"))
    H=lolo.data
    D=H[0,:,:]
    H=D
    #print corps.interquartile_range(H)
    #print np.size(H,0),np.size(H,1)
    list_exp=ser_xml.exceptions(root)
    raxep=[]
    while boc<len(tabp)-1:
        boc+=1
        a=str(tabp[boc][2])
        #print a,boc
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
        #print tabp[r1[0]:r1[len(r1)-1]+1]
        #raw_input()
        #print len(r1)
        nt=tabp[r1[0]:r1[len(r1)-1]+1]
        #print nt
        #print 'Start'
        test=[]
        testc3=[]
        con=0
        ampref=[]
        amp=[]
        newrt=[]
        #print newn
        for jj in range(len(nt)):
            #print newn
            #print jj
            #print nt[jj][0],nt[jj][1]
            #print nt[jj][0],nt[jj][1]
            #print wa7da
            r,indi=lawij.alig(nt[jj][0],nt[jj][1],H,nt[jj][3],wa7da)
            
            print r,indi,nt[jj][2],nt[jj][0]-r[0],nt[jj][1]-r[1]
            #raw_input()
            if indi==0 :
                con=con+1
                
                if np.abs(r[0])==4 or np.abs(r[1])==6:
                    testc3.append((1))
                    
                else:
                    testc3.append((0))
                test.append((0))
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                zayneb=H[(nt[jj][0]-r[0])-3:(nt[jj][0]-r[0])+4,(nt[jj][1]-r[1])-3:(nt[jj][1]-r[1])+4]*1.
                nr=f.subpix2(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                zayneb=H[(nt[jj][0]-r[0])-9:(nt[jj][0]-r[0])+10,(nt[jj][1]-r[1])-9:(nt[jj][1]-r[1])+10]*1.
                #nr=f.subpix(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                chl=f.dephcl(zayneb)
                chg=f.dephcg(zayneb)
                #ch,congl=f.dephc(zayneb)
                ch,congl=f.dephcaprio(zayneb,float(nt[jj][4]),float(nt[jj][5]),nt[jj][6])
                #print congl
                #cax = ax.imshow(zayneb, interpolation='nearest')
                #plt.show()
                #print ch
                #plt.show()
                ampref.append(float(nt[jj][3])*1.)
                amp.append( H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.)
                #print str(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])])
                newrt.append((nr[0],nr[1],nt[jj][2],str(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]),chg[0],chg[1],chl[0],chl[1],ch[0],ch[1],congl,r[0],r[1]))
                #print newrt
            else:
                if r[0]<100 :
                    if np.abs(r[0])==3 or np.abs(r[1])==5:
                        testc3.append((1))
                        
                    else:
                        testc3.append((0))
                    con=con+1
                    test.append((1))
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111)
                    zayneb=H[(nt[jj][0]-r[0])-3:(nt[jj][0]-r[0])+4,(nt[jj][1]-r[1])-3:(nt[jj][1]-r[1])+4]*1.
                    nr=f.subpix2(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                    #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                    zayneb=H[(nt[jj][0]-r[0])-9:(nt[jj][0]-r[0])+10,(nt[jj][1]-r[1])-9:(nt[jj][1]-r[1])+10]*1.
                    #nr=f.subpix(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                    #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                    chl=f.dephcl(zayneb)
                    chg=f.dephcg(zayneb)
                    #ch,congl=f.dephc(zayneb)
                    ch,congl=f.dephcaprio(zayneb,float(nt[jj][4]),float(nt[jj][5]),nt[jj][6])
                    ampref.append(float(nt[jj][3])*1.)
                    amp.append(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.)
                    newrt.append((nr[0],nr[1],nt[jj][2],H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.,chg[0],chg[1],chl[0],chl[1],ch[0],ch[1],congl,r[0],r[1]))
                else:
                    test.append((2))
                    testc3.append((2))
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111)
                    zayneb=H[nt[jj][0]-7:(nt[jj][0])+8,(nt[jj][1])-7:(nt[jj][1])+8]*1.
                    ampref.append(float(nt[jj][3])*1.)
                    amp.append(0)
                    #cax = ax.imshow(zayneb, interpolation='nearest')
                    #plt.show()
                    #raw_input()
                    newrt.append((0,0,0,0,0,0,0,0,0,0,0,0,0))

        #raw_input()
        #print newn
        #print ampref
        #print amp,newn,testc3
        #print nt
        #print newrt
        #raw_input()
        #amptest=np.nonzero(amp>0)
        o=np.nonzero(testc3==0)
        vamp=np.array(amp)
        ivamref=vamp[o]
        o=np.nonzero(ivamref>210000)
        #print newn,'test'
        if (float(len(o[0]))*1.000000001)/float(len(ivamref)*1.00000001)>0.4 or f.exp_hand(list_exp,newn)==1:
            #print newn,'d5aal'
            if len(nt)==con:
                if len(nt)==1:
                    raxep.append(newrt[0])
                    #print 'accepted'
                else:
                    artestc3=np.array(testc3)
                    olc3=np.nonzero(artestc3==1)
                    if np.size(olc3)>0:
                        if f.exp_ync(amp,ampref,testc3)==1:
                            #print 'ouffffff'
                            artest=np.array(test)
                            ol=np.nonzero(artest==1)
                            o=np.nonzero(artest==0)
                            if np.size(ol)>0:
                                #print 'accepted with some conditions'
                                #print f.exp_ync(amp,ampref,testc3)
                                #if f.exp_ync(amp,ampref,testc3)==0:
                                    #print 'llllllaaaaaaaaa'
                                if f.exp_yn(amp,ampref,test)==1:
                                    for no in range(len(newrt)):
                                        raxep.append(newrt[no])
                                elif f.exp_hand(list_exp,newn)==1:
                                    artest=np.array(test)
                                    rnt=np.array(newrt)
                                    ol=np.nonzero(artest==1)
                                    o=np.nonzero(artest==0)
                                    vo=rnt[o]
                                    #raw_input()
                                    for no in range(len(vo)):
                                        #print '%%%%%%%%%'
                                        zi=lawij.mod(vo[no])
                                        #print zi
                                        #print artest
                                        #raw_input()
                                        raxep.append((zi))
                            else:
                                #print f.exp_ync(amp,ampref,testc3)
                                for no in range(len(nt)):
                                    raxep.append(newrt[no])
                                #print 'accepted'                                    
                        else:
                            #print 'ouuuuut'
                            #print f.exp_hand(list_exp,newn)
                            if f.exp_hand(list_exp,newn)==1 :
                                artest=np.array(test)
                                artestc3=np.array(testc3)
                                rnt=np.array(newrt)
                                #print nt
                                #raw_input()
                                ol=np.nonzero(artest==1)
                                condlist=[artest==0,artestc3==0]
                                g=condlist[0]*condlist[1]
                                #print g,artest,artestc3
                                o=np.nonzero(artest==0)
                                #print '%%%hhh%%%%%%'
                                #print rnt
                                #print test
                                vo=rnt[g]
                                #print vo,'%%%hhh%%%%%%'
                                #raw_input()
                                for no in range(len(vo)):
                                            #print '%%%%%%%%%'
                                            zi=lawij.mod(vo[no])
                                            #print zi
                                            #print artest
                                            #raw_input()
                                            raxep.append((zi))
                    else:
                                #print f.exp_ync(amp,ampref,testc3)
                                for no in range(len(nt)):
                                    raxep.append(newrt[no])
                                #print 'accepted'
    
            else:
                if f.exp_hand(list_exp,newn)==1 :
                    artest=np.array(test)
                    #print newrt
                    rnt=xml_gen.cell2tab8(newrt)
                    #print nt
                    #raw_input()
                    ol=np.nonzero(artest==1)
                    o=np.nonzero(artest==0)
                    #print '%%%hhh%%%%%%'
                    rnt=np.array(rnt)
                    #print test
                    vo=rnt[o]
                    #print vo
                    #print len(vo)
                    #raw_input()
                    for no in range(len(vo)):
                                #print '%%%%%%%%%'
                                zi=lawij.mod(vo[no])
                                #print zi
                                #print artest
                                #raw_input()
                                raxep.append((zi))
                    #print 'may be...'
                #else:
                    #print 'refused without discussion'    
            #print test
            
                
            
            #print 'DONE'
            #print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            #print '        '
            #raw_input()
    locale.setlocale(locale.LC_ALL, old_locale)
    return raxep,tabc
    #for lll in range(len(raxep)):
        #print  raxep[lll]  






def rec(root):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    #root='/home/miv/belghith/Bureau/KAROM/Akram/nmr/RP/5/pdata/1'
    tabp,tabc=gen_in.ini2(root)
    
    
    boc=-1
    lolo= load(os.path.join(root,"2rr"))
    H=lolo.data
    D=H[0,:,:]
    H=D
    #print corps.interquartile_range(H)
    #print np.size(H,0),np.size(H,1)
    list_exp=ser_xml.exceptions(root)
    raxep=[]
    while boc<len(tabp)-1:
        boc+=1
        a=str(tabp[boc][2])
        #print a,boc
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
        #print tabp[r1[0]:r1[len(r1)-1]+1]
        #raw_input()
        #print len(r1)
        nt=tabp[r1[0]:r1[len(r1)-1]+1]
        #print nt
        #print 'Start'
        test=[]
        testc3=[]
        con=0
        ampref=[]
        amp=[]
        newrt=[]
        #print newn
        for jj in range(len(nt)):
            #print newn
            #print jj
            #print nt[jj][0],nt[jj][1]
            #print nt[jj][0],nt[jj][1]
            r,indi=lawij.aligrec(nt[jj][0],nt[jj][1],H,nt[jj][3])
            
            #print r,indi,nt[jj][2],nt[jj][0]-r[0],nt[jj][1]-r[1]
            #raw_input()
            if indi==0 :
                con=con+1
                
                #print 'ok'
                testc3.append((0))
                test.append((0))
                #fig = plt.figure()
                #ax = fig.add_subplot(111)
                zayneb=H[(nt[jj][0]-r[0])-3:(nt[jj][0]-r[0])+4,(nt[jj][1]-r[1])-3:(nt[jj][1]-r[1])+4]*1.
                nr=f.subpix2(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                zayneb=H[(nt[jj][0]-r[0])-9:(nt[jj][0]-r[0])+10,(nt[jj][1]-r[1])-9:(nt[jj][1]-r[1])+10]*1.
                #nr=f.subpix(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                chl=f.dephcl(zayneb)
                chg=f.dephcg(zayneb)
                #ch,congl=f.dephc(zayneb)
                ch,congl=f.dephcaprio(zayneb,float(nt[jj][4]),float(nt[jj][5]),nt[jj][6])
                #print congl
                #cax = ax.imshow(zayneb, interpolation='nearest')
                #plt.show()
                #print ch
                #plt.show()
                ampref.append(float(nt[jj][3])*1.)
                amp.append( H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.)
                #print str(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])])
                newrt.append((nr[0],nr[1],nt[jj][2],str(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]),chg[0],chg[1],chl[0],chl[1],ch[0],ch[1],congl,r[0],r[1]))
                #print newrt
            else:
                if r[0]<100 :
                    if np.abs(r[0])==3 or np.abs(r[1])==5:
                        testc3.append((1))
                        
                    else:
                        testc3.append((0))
                    con=con+1
                    test.append((1))
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111)
                    zayneb=H[(nt[jj][0]-r[0])-3:(nt[jj][0]-r[0])+4,(nt[jj][1]-r[1])-3:(nt[jj][1]-r[1])+4]*1.
                    nr=f.subpix2(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                    #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                    zayneb=H[(nt[jj][0]-r[0])-9:(nt[jj][0]-r[0])+10,(nt[jj][1]-r[1])-9:(nt[jj][1]-r[1])+10]*1.
                    #nr=f.subpix(zayneb,(nt[jj][0]-r[0]),(nt[jj][1]-r[1]))
                    #print (nt[jj][0]-r[0]),(nt[jj][1]-r[1]),nr
                    chl=f.dephcl(zayneb)
                    chg=f.dephcg(zayneb)
                    #ch,congl=f.dephc(zayneb)
                    ch,congl=f.dephcaprio(zayneb,float(nt[jj][4]),float(nt[jj][5]),nt[jj][6])
                    ampref.append(float(nt[jj][3])*1.)
                    amp.append(H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.)
                    newrt.append((nr[0],nr[1],nt[jj][2],H[(nt[jj][0]-r[0]),(nt[jj][1]-r[1])]*1.,chg[0],chg[1],chl[0],chl[1],ch[0],ch[1],congl,r[0],r[1]))
                else:
                    test.append((2))
                    testc3.append((2))
                    #fig = plt.figure()
                    #ax = fig.add_subplot(111)
                    zayneb=H[nt[jj][0]-7:(nt[jj][0])+8,(nt[jj][1])-7:(nt[jj][1])+8]*1.
                    ampref.append(float(nt[jj][3])*1.)
                    amp.append(0)
                    #cax = ax.imshow(zayneb, interpolation='nearest')
                    #plt.show()
                    #raw_input()
                    newrt.append((0,0,0,0,0,0,0,0,0,0,0,0,0))

        #raw_input()
        #print newn
        #print ampref
        #print amp,newn,testc3
        #print nt
        #print newrt
        #raw_input()
        #amptest=np.nonzero(amp>0)
        o=np.nonzero(testc3==0)
        vamp=np.array(amp)
        ivamref=vamp[o]
        o=np.nonzero(ivamref>100)
        #print newn,'test'
        if (float(len(o[0]))*1.000000001)/float(len(ivamref)*1.00000001)>0.4 or f.exp_hand(list_exp,newn)==1:
            #print newn,'d5aal'
            if len(nt)==con:
                if len(nt)==1:
                    raxep.append(newrt[0])
                    #print 'accepted'
                else:
                    artestc3=np.array(testc3)
                    olc3=np.nonzero(artestc3==1)
                    if np.size(olc3)>0:
                        if f.exp_ync(amp,ampref,testc3)==1:
                            #print 'ouffffff'
                            artest=np.array(test)
                            ol=np.nonzero(artest==1)
                            o=np.nonzero(artest==0)
                            if np.size(ol)>0:
                                #print 'accepted with some conditions'
                                #print f.exp_ync(amp,ampref,testc3)
                                #if f.exp_ync(amp,ampref,testc3)==0:
                                    #print 'llllllaaaaaaaaa'
                                if f.exp_yn(amp,ampref,test)==1:
                                    for no in range(len(newrt)):
                                        raxep.append(newrt[no])
                                elif f.exp_hand(list_exp,newn)==1:
                                    artest=np.array(test)
                                    rnt=np.array(newrt)
                                    ol=np.nonzero(artest==1)
                                    o=np.nonzero(artest==0)
                                    vo=rnt[o]
                                    #raw_input()
                                    for no in range(len(vo)):
                                        #print '%%%%%%%%%'
                                        zi=lawij.mod(vo[no])
                                        #print zi
                                        #print artest
                                        #raw_input()
                                        raxep.append((zi))
                            else:
                                #print f.exp_ync(amp,ampref,testc3)
                                for no in range(len(nt)):
                                    raxep.append(newrt[no])
                                #print 'accepted'                                    
                        else:
                            #print 'ouuuuut'
                            #print f.exp_hand(list_exp,newn)
                            if f.exp_hand(list_exp,newn)==1 :
                                artest=np.array(test)
                                artestc3=np.array(testc3)
                                rnt=np.array(newrt)
                                #print nt
                                #raw_input()
                                ol=np.nonzero(artest==1)
                                condlist=[artest==0,artestc3==0]
                                g=condlist[0]*condlist[1]
                                #print g,artest,artestc3
                                o=np.nonzero(artest==0)
                                #print '%%%hhh%%%%%%'
                                #print rnt
                                #print test
                                vo=rnt[g]
                                #print vo,'%%%hhh%%%%%%'
                                #raw_input()
                                for no in range(len(vo)):
                                            #print '%%%%%%%%%'
                                            zi=lawij.mod(vo[no])
                                            #print zi
                                            #print artest
                                            #raw_input()
                                            raxep.append((zi))
                    else:
                                #print f.exp_ync(amp,ampref,testc3)
                                for no in range(len(nt)):
                                    raxep.append(newrt[no])
                                #print 'accepted'
    
            else:
                if f.exp_hand(list_exp,newn)==1 :
                    artest=np.array(test)
                    #print newrt
                    rnt=xml_gen.cell2tab8(newrt)
                    #print nt
                    #raw_input()
                    ol=np.nonzero(artest==1)
                    o=np.nonzero(artest==0)
                    #print '%%%hhh%%%%%%'
                    rnt=np.array(rnt)
                    #print test
                    vo=rnt[o]
                    #print vo
                    #print len(vo)
                    #raw_input()
                    for no in range(len(vo)):
                                #print '%%%%%%%%%'
                                zi=lawij.mod(vo[no])
                                #print zi
                                #print artest
                                #raw_input()
                                raxep.append((zi))
                    #print 'may be...'
                #else:
                    #print 'refused without discussion'    
            #print test
            
                
            
            #print 'DONE'
            #print '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
            #print '        '
            #raw_input()
    locale.setlocale(locale.LC_ALL, old_locale)
    return raxep,tabc
    #for lll in range(len(raxep)):
        #print  raxep[lll]  
