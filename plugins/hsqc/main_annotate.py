"""
Main function
"""
import ser_xml
import xml_gen
import annotate
from shutil import copytree
from shutil import rmtree
import os
import gen_in
import fun
import numpy as np
def mainanno(path):
    #z=ser_xml.exceptions('except.txt')
    #print len(z[0])
    #if str(z[0])=="akram":
        #print 'good'
    #root='/home/miv/belghith/Bureau/KAROM/Akram/nmr/PP/5/pdata/1'
    for ll in range(len(path)):
        root=path[ll]
        nroot=root[0:len(root)-10]
        nroot=nroot+'new'
        ext= root[len(root)-9:len(root)]
        #raw_input()
        if os.path.isdir(nroot)==False:
            copytree(root[0:len(root)-10],nroot)
        else:
             rmtree(nroot)
             copytree(root[0:len(root)-10],nroot)
        #print nroot
        nnroot=os.path.join(nroot,ext)
        raxep=annotate.ann(nnroot)
        xml_gen.gen_s_ch(raxep,nnroot,'peaklistwithch.xml')
        tab=ser_xml.read(os.path.join(nnroot,"peaklistwithch.xml"))
        #print tab
        tabc=xml_gen.cell2tab(tab)
        ntab=gen_in.fini(tabc,nnroot)
        rntab=np.array(ntab)
        #condlist=[np.abs(rntab[:,0]-22.7)<0.3,np.abs(rntab[:,0]-1.33)<0.03]
        #g=condlist[0]*condlist[1]
        #print g
        for h in range(np.size(rntab[:,0])):
            if np.abs(float(rntab[h,0])-22.7)<0.3 and np.abs(float(rntab[h,1])-1.33)<0.03 :
                dc=float(rntab[h,0])-22.7000000000
                dh=float(rntab[h,1])-1.330000000
                break
        #print 22.7*np.ones([1,np.size(rntab[:,0])])[0]-100
        #print rntab[:,0]-1
        fun.update(nnroot,'procs',dh)
        fun.update(nnroot,'proc2s',dc)
        ntab=gen_in.fini(tabc,nnroot)
        xml_gen.gen_s(ntab,nnroot,'peaklist.xml')