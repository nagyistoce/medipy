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
import func
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
        if nroot[len(nroot)-1]=='/':
            nroot=root[0:len(root)-11]
            nroot=nroot+'new'
            ext= root[len(root)-10:len(root)]
            #raw_input()
            if os.path.isdir(nroot)==False:
                copytree(root[0:len(root)-11],nroot)
            else:
                rmtree(nroot)
                copytree(root[0:len(root)-11],nroot)
        else:
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
        raxep,ttt=annotate.rec(nnroot)
        tab=ser_xml.read(os.path.join(nnroot,"ref.xml"))
        tabc=xml_gen.cell2tab(tab)
        #print tabc[0][3]
        ntab=gen_in.fini(raxep,nnroot)
        rntab=np.array(ntab)
        #print rntab
        dc=float(rntab[0,0])-float(tabc[0][0])
        dh=float(rntab[0,1])-float(tabc[0][1])
        
        
        
        
        fun.update(nnroot,'procs',dh)
        fun.update(nnroot,'proc2s',dc)
        
        raxep,ttt=annotate.ann(nnroot,float(tabc[0][3]))
        
        raxep=func.msatkin(raxep)
        #xml_gen.gen_s_chggen(ttt,'/home/miv/belghith/Bureau','peaknewuuuu.xml')
        #xml_gen.gen_s_ch(raxep,nnroot,'peaklistwithch.xml')
        #tab=ser_xml.read(os.path.join(nnroot,"peaklistwithch.xml"))
        #print tab
        #tabc=xml_gen.cell2tab(tab)
        ntab=gen_in.fini(raxep,nnroot)
        rntab=np.array(ntab)
        #condlist=[np.abs(rntab[:,0]-22.7)<0.3,np.abs(rntab[:,0]-1.33)<0.03]
        #g=condlist[0]*condlist[1]
        ntab=gen_in.fini(raxep,nnroot)
        #rntab=np.array(raxep)
        #nrex=raxep
        #nrex[0:1,:]=rntab[0:1,:]
        #xml_gen.gen_s_ch(nrex,nnroot,'peaklistwithchppm.xml')
        xml_gen.gen_s_chg(xml_gen.tabp2tabg(ntab,raxep),nnroot,'peaklistwithchg.xml')
        xml_gen.gen_s_chl(xml_gen.tabp2tabl(ntab,raxep),nnroot,'peaklistwithchl.xml')
        xml_gen.gen_s_chgl(xml_gen.tabp2tabgl(ntab,raxep),nnroot,'peaklistwithchlgl.xml')
        xml_gen.gen_s(ntab,nnroot,'peaklist.xml')