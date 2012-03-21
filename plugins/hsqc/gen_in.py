"""
Converting steps
"""
import ser_xml as ser
import numpy as np
import xml_gen as xml
import conv
import ordon
import os
def ini(root):
    CA,CB= conv.convc(root)
    HA,HB= conv.convh(root)
    rooti="peak.xml"
    nom=os.path.join(root,rooti)
    z=ser.read(nom)
    tabc=ordon.gen_ord(z)
    K=len(tabc)
    tabp=[]
    #print CA,CB,HA,HB
    #print (1/float(CA))*(float(tabc[0][0])-float(CB))
    #print tabc[0][0]
    for k in range(K):
        ta0=np.floor((1/float(CA))*(float(tabc[k][0])-float(CB)))
        ta1=np.floor((1/float(HA))*(float(tabc[k][1])-float(HB)))
        tabp.append((ta0,ta1,tabc[k][2],tabc[k][3]))
    
    return tabp,tabc

def haya(root):
    CA,CB= conv.convc(root)
    HA,HB= conv.convh(root)
    rooto="peaklist.xml"
    rooti="peak.xml"
    nom=os.path.join(root,rooti)
    z=ser.read(nom)
    tabc=ordon.gen_ord(z)
    nom=os.path.join(root,rooto)
    z=ser.read(nom)
    tabd=ordon.gen_ord(z)
    K=len(tabc)
    tabamp=[]
    #print CA,CB,HA,HB
    #print (1/float(CA))*(float(tabc[0][0])-float(CB))
    #print tabc[0][0]
    i=0
    for k in range(K):
        if tabc[k][2]==tabd[i][2]:            
            tabamp.append(float(tabd[i][3]))
            i+=1
        else:
            tabamp.append(0)
    
    return tabamp
def fini(tabc,root):
    #root='/home/miv/belghith/Bureau/KAROM/Akram/nmr/PH/5/pdata/1'
    CA,CB= conv.convc(root)
    HA,HB= conv.convh(root)
    K=len(tabc)
    tabp=[]
    #print CA,CB,HA,HB
    #print (1/float(CA))*(float(tabc[0][0])-float(CB))
    #print tabc[0][0]
    for k in range(K):
        #ta0=np.floor((1/float(CA))*(float(tabc[k][0])-float(CB)))
        ta0=((float(CA))*(float(tabc[k][0])+1)+float(CB))
        #ta1=np.floor((1/float(HA))*(float(tabc[k][1])-float(HB)))
        ta1=((float(HA))*(float(tabc[k][1])+1)+float(HB))
        tabp.append((ta0,ta1,tabc[k][2],tabc[k][3]))
    
    return tabp