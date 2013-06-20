"""
Functions used in ppm2pixel converting
"""
import numpy as np
import os
from decimal import Decimal, ROUND_HALF_UP
import shutil
def offset(path,type):
        C=os.path.join(path,type)
        f=open(C, 'r')
        l=f.readlines()
        C=l[63][11:]
        f.close()
        C = Decimal(C)
        C=C.quantize(Decimal('.0001'), rounding=ROUND_HALF_UP)
        return C
    
def SW(path,type):
        C=os.path.join(path,type)
        f=open(C, 'r')
        l=f.readlines()
        #print l[87][3:5]
        if l[87][3:5]=='SW':
            C=l[87][9:]
        else:
            C=l[89][9:]
        f.close()
        C = Decimal(C)
        C=C.quantize(Decimal('.0001'), rounding=ROUND_HALF_UP)
        return C

def SI(path,type):
        C=os.path.join(path,type)
        f=open(C, 'r')
        l=f.readlines()
        if l[78][3:5]=='SI':
            C=l[78][7:]
        else:
            C=l[79][7:]
        f.close()
        C = Decimal(C)
        C=C.quantize(Decimal('.01'), rounding=ROUND_HALF_UP)
        return C
def SF(path,type):
        C=os.path.join(path,type)
        f=open(C, 'r')
        l=f.readlines()
        if l[77][3:5]=='SF':
            C=l[77][7:]
        else:
            C=l[78][7:]
        f.close()
        C = Decimal(C)
        C=C.quantize(Decimal('.0001'), rounding=ROUND_HALF_UP)
        return C   
def update(root,type,dc):
    f=open(os.path.join(root,type), 'r')
    g=open(os.path.join(root,'temp'), 'w')
    l=f.readlines()
    for h in range(np.size(l[:])):
            if h==63:
                    m=l[63][0:11]+str(float(l[63][11:])-dc)+'\n'
                    l1=g.writelines(m)
            else:
                    l1=g.writelines(l[h])
    f.close()
    g.close()
    shutil.move(os.path.join(root,'temp'),os.path.join(root,type))
