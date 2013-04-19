import numpy as np
import os
import fonc_util as f
#import segmentation as seg
import gen_in
import xml_gen
import fark
def runmainc(path1,path2):
#if __name__=="__main__":
    #base_directory='/home/miv/belghith/Bureau/group1'
    base_directory=path1
    i=0
    k={}
    for root,dirs,files in os.walk(base_directory):
        if '2rr'in files:
            k[i]=str(root)
            i=i+1
    #base_directory='/home/miv/belghith/Bureau/group2'
    base_directory=path2
    i=0
    k1={}
    for root,dirs,files in os.walk(base_directory):
        if '2rr'in files:
            k1[i]=str(root)
            i=i+1
    #print k,k1
    #raw_input()
    Xt=[]
    for ll in range(len(k)):
        root=k[ll]
        tabc=gen_in.haya(root)
        #tabc=f.tab2tab(tabc)
        #tabc=np.array(tabc)
        #print tabc.T#np.random.rand(4, 50)
        #raw_input()
        Xt.append(tabc)  
        A=np.array(Xt) 
    Xt=[]
    for ll in range(len(k1)):
        root=k[ll]
        tabc=gen_in.haya(root)
        #tabc=f.tab2tab(tabc)
        #tabc=np.array(tabc)
        #print tabc.T#np.random.rand(4, 50)
        #raw_input()
        Xt.append(tabc)  
        B=np.array(Xt) 
    #print A, B 
    #raw_input()
    tabp,tabc=gen_in.ini(k1[0])
    inm,inp=f.ser_mat(tabp)
    #Xt = np.random.rand(10, np.size(tabp,0))
    #Xo = np.random.rand(10, np.size(tabp,0))+5
    #print k,k1
    #raw_input()
    a=fark.detfark(A,B,inm,inp)
    xml_gen.gen_change(a,base_directory,'res_change.txt')
    #print a,a[1][0],len(a)
