"""
Corpus generating
"""
import ser_xml as ser
import numpy as np
import os
import locale

from medipy.base import find_resource

def gen(nom):
    
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    tab=ser.read_s()
    k=len(tab)/4
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open('nom', 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h,0])
        m3='" F2="'
        m4=str(tab[h,1])
        m5='" annotation="'
        m6=str(tab[h,2])
        m7='" intensity="'
        m8=str(tab[h,3])
        m9='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()
    
    locale.setlocale(locale.LC_ALL, old_locale)
    
def gen_s_ch(tab,root,nom):
    
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    #tab=ser.read_s()
    k=len(tab)
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open(os.path.join(root,nom), 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h][0])
        m3='" F2="'
        m4=str(tab[h][1])
        m5='" annotation="'
        m6=str(tab[h][2])
        m7='" intensity="'
        m8=str(tab[h][3])
        m9='" depHg="'
        m10=str(tab[h][4])
        m11='" depCg="'
        m12=str(tab[h][5])
        m13='" depHl="'
        m14=str(tab[h][6])
        m15='" depCl="'
        m16=str(tab[h][7])
        m17='" depHgl="'
        m18=str(tab[h][8])
        m19='" depCgl="'
        m20=str(tab[h][9])
        m21='" gl="'
        m22=str(tab[h][10])
        m23='" depc="'
        m24=str(tab[h][11])
        m25='" deph="'
        m26=str(tab[h][12])
        m27='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17+m18+m19+m20+m21+m22+m23+m24+m25+m26+m27
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()  
    locale.setlocale(locale.LC_ALL, old_locale)
    
    
def gen_s_chgl(tab,root,nom):
    
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    #tab=ser.read_s()
    k=len(tab)
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open(os.path.join(root,nom), 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h][0])
        m3='" F2="'
        m4=str(tab[h][1])
        m5='" annotation="'
        m6=str(tab[h][2])
        m7='" intensity="'
        m8=str(tab[h][3])
        m9='" depH="'
        m10=str(tab[h][4])
        m11='" depC="'
        m12=str(tab[h][5])
        m13='" gl="'
        m14=str(tab[h][6])
        m15='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m13+m14+m15
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()  
    locale.setlocale(locale.LC_ALL, old_locale)    
def gen_s_chg(tab,root,nom):
    
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    #tab=ser.read_s()
    k=len(tab)
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open(os.path.join(root,nom), 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h][0])
        m3='" F2="'
        m4=str(tab[h][1])
        m5='" annotation="'
        m6=str(tab[h][2])
        m7='" intensity="'
        m8=str(tab[h][3])
        m9='" depH="'
        m10=str(tab[h][4])
        m11='" depC="'
        m12=str(tab[h][5])
        m13='" gl="'
        m14=str('g')
        m15='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m13+m14+m15
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()  
    locale.setlocale(locale.LC_ALL, old_locale)    

def gen_s_chggen(tab,root,nom):
    
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    print "d5allllll"
    #tab=ser.read_s()
    k=len(tab)
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open(os.path.join(root,nom), 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h][0])
        m3='" F2="'
        m4=str(tab[h][1])
        m5='" annotation="'
        m6=str(tab[h][2])
        m7='" intensity="'
        m8=str(tab[h][3])
        m9='" etalH="'
        m10="0.3"
        m11='" etalC="'
        m12="0.3"
        m13='" gl="'
        m14=str('l')
        m15='" important="'
        m16=str('y')
        m17='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m13+m14+m15+m16+m17
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()  
    locale.setlocale(locale.LC_ALL, old_locale)    
   
def gen_s_chl(tab,root,nom):
    
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    #tab=ser.read_s()
    k=len(tab)
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open(os.path.join(root,nom), 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h][0])
        m3='" F2="'
        m4=str(tab[h][1])
        m5='" annotation="'
        m6=str(tab[h][2])
        m7='" intensity="'
        m8=str(tab[h][3])
        m9='" depH="'
        m10=str(tab[h][4])
        m11='" depC="'
        m12=str(tab[h][5])
        m13='" gl="'
        m14=str('l')
        m15='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m13+m14+m15
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()  
    locale.setlocale(locale.LC_ALL, old_locale)        
def cell2tab(tabc):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    tabb=[]
    k=len(tabc)/4
    for h in range(k):
        tabb.append((tabc[h,0],tabc[h,1],tabc[h,2],tabc[h,3]))
    locale.setlocale(locale.LC_ALL, old_locale)
    return tabb

def cell2tab1(tabc):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    tabb=[]
    k=len(tabc)/8

    for h in range(k):
        tabb.append((tabc[h,0],tabc[h,1],tabc[h,2],tabc[h,3],tabc[h,4],tabc[h,5],tabc[h,6],tabc[h,7]))
    locale.setlocale(locale.LC_ALL, old_locale)
    return tabb
def cell2tab8(tabc):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    tabb=[]
    k=len(tabc)/4
    #print range(len(tabc)),tabc[1][0]
    for h in range(len(tabc)):
        tabb.append((tabc[h][0],tabc[h][1],tabc[h][2],tabc[h][3],tabc[h][4],tabc[h][5],tabc[h][6],tabc[h][7],tabc[h][8],tabc[h][9],tabc[h][10],tabc[h][11],tabc[h][12]))
    locale.setlocale(locale.LC_ALL, old_locale)
    return tabb

def tabp2tabg(tab1,tab2):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    ntab=[]
    for h in range(len(tab1)):
        ntab.append((tab2[h][0],tab2[h][1],tab1[h][2],tab1[h][3],tab2[h][4],tab2[h][5]))

    locale.setlocale(locale.LC_ALL, old_locale)
    return ntab

def tabp2tabl(tab1,tab2):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    ntab=[]
    for h in range(len(tab1)):
        ntab.append((tab2[h][0],tab2[h][1],tab1[h][2],tab1[h][3],tab2[h][6],tab2[h][7]))

    locale.setlocale(locale.LC_ALL, old_locale)
    return ntab

def tabp2tabgl(tab1,tab2):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    ntab=[]
    for h in range(len(tab1)):
        ntab.append((tab2[h][0],tab2[h][1],tab1[h][2],tab1[h][3],tab2[h][8],tab2[h][9],tab2[h][10]))

    locale.setlocale(locale.LC_ALL, old_locale)
    return ntab


def gen_s(tab,root,nom):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    #tab=ser.read_s()
    k=len(tab)
    f=open(find_resource(os.path.join("resources", "components", "belghith", 'entete.txt')), 'r')
    l=f.readlines()
    f.close()
    f=open(os.path.join(root,nom), 'w')
    l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='<Peak2D F1="'
        m2=str(tab[h][0])
        m3='" F2="'
        m4=str(tab[h][1])
        m5='" annotation="'
        m6=str(tab[h][2])
        m7='" intensity="'
        m8=str(tab[h][3])
        m9='" type="0"/>\n'
        m=m1+m2+m3+m4+m5+m6+m7+m8+m9
        l1=f.writelines(m)
    l1=f.writelines('</PeakList2D> \n')
    l1=f.writelines('</PeakList> \n')
    f.close()
    
    locale.setlocale(locale.LC_ALL, old_locale)
    
    
    
def gen_change(tab,root,nom):
    old_locale = locale.getlocale()
    locale.setlocale(locale.LC_ALL, "C")
    
    #tab=ser.read_s()
    k=len(tab)
    f=open(os.path.join(root,nom), 'w')
    #l1=f.writelines(l)
    #l1=f.writelines('\n')
    for h in range(k):
        m1='"metabolite"= '
        m2=str(tab[h][0])
        m3=', "change degree"= '
        m4=str(tab[h][1])
        m5=', "concentration"= '
        m6=str(tab[h][2])
        m=m1+m2+m3+m4+m5+m6
        l1=f.writelines(m)
        l1=f.writelines('\n')
    f.close()
    
    locale.setlocale(locale.LC_ALL, old_locale)