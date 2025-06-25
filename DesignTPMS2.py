'''
    DesignTPMS.py
    
    Create Solid TPMS Structure
    Build STL
    
    Adapted from Matlab to Python by TMR on 23 January 2024
    Original Author: Fabian Guenther et al.
'''
#%% Imports
import numpy as np
from math import *
from itertools import product
#%% Params
from params import *
#%% Transformation
from transform import T
#%% Functions
def TPMS(form):
    if form=='P':
        return lambda x,y,z: -1*(cos(2*pi*x)+cos(2*pi*y)+cos(2*pi*z))
        #return lambda x,y,z: -1*(cos(2*pi*x/labda[0]+theta[0])+cos(2*pi*y/labda[1]+theta[1])+cos(2*pi*z/labda[2]+theta[2]))
    if form=='G':
        return lambda x,y,z: cos(2*pi*x)*sin(2*pi*y)+cos(2*pi*y)*sin(2*pi*z)+cos(2*pi*z)*sin(2*pi*x)
        #return lambda x,y,z: cos(2*pi*x/labda[0]+theta[0])*sin(2*pi*y/labda[1]+theta[1])+cos(2*pi*y/labda[1]+theta[1])*sin(2*pi*z/labda[2]+theta[2])+cos(2*pi*z/labda[2]+theta[2])*sin(2*pi*x/labda[0]+theta[0])
    if form=='D':
        #return lambda x,y,z: sin(labda[0]*2*pi*x+theta[0])*sin(labda[1]*2*pi*y+theta[1])*sin(labda[2]*2*pi*z+theta[2]) + sin(labda[0]*2*pi*x+theta[0])*cos(labda[1]*2*pi*y+theta[1])*cos(labda[2]*2*pi*z+theta[2]) + cos(labda[0]*2*pi*x+theta[0])*sin(labda[1]*2*pi*y+theta[1])*cos(labda[2]*2*pi*z+theta[2]) + cos(labda[0]*2*pi*x+theta[0])*cos(labda[1]*2*pi*y+theta[1])*sin(labda[2]*2*pi*z+theta[2])
        return lambda x,y,z: sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)+sin(2*pi*x)*cos(2*pi*y)*cos(2*pi*z)+cos(2*pi*x)*sin(2*pi*y)*cos(2*pi*z)+cos(2*pi*x)*cos(2*pi*y)*sin(2*pi*z)
    if form=='I':
        return lambda x,y,z: 2*(cos(2*pi*x)*cos(2*pi*y)+cos(2*pi*y)*cos(2*pi*z)+cos(2*pi*z)*cos(2*pi*z))-(cos(4*pi*x)+cos(4*pi*y)+cos(4*pi*z))
    print('ERROR')
def calcVF(F,r,X,V0,dV):
    V=0
    # xs=[];ys=[];zs=[];
    n=len(X)
    for i in range(n):
        x,y,z=X[i]
        if abs(F(x,y,z))<r:
            V+=dV
            # print(V)
            # xs.append(x)
            # ys.append(y)
            # zs.append(z)
    
    return V/V0
def NR(F,r0,X,V0,dV,goal,tol=None):
    if tol==None: tol=1*dV/V0
    dr=0.01#10/len(X)**(1/3)
    vf=calcVF(F,r0,X,V0,dV)
    err=vf-goal
    #print(err)
    #print(tol)
    if abs(err)<tol:
        return r0
    else:
        #g=0
        #while g==0:
        g=(calcVF(F,r0+dr,X,V0,dV)-calcVF(F,r0-dr,X,V0,dV))/(2*dr)
        if g==0: 
        	g=2/pi
        	#tol=tol*2
        r=r0-err/g
        print("r="+str(r))
        return NR(F,r,X,V0,dV,goal,tol)
#%% Vfrac Calibration
# this is independant of the transformations

rstart = -0.5 # initial guess
volTol = 0.0001 # solver tolerance

N=25
xi = np.linspace(0,2*pi,N)
yi = np.linspace(0,2*pi,N)
zi = np.linspace(0,2*pi,N)

X=list(product(xi,yi,zi))

F = TPMS(form)

V0=pi**3
dV=(pi**3)/(N**3)

r=NR(F,rstart,X,V0,dV,vf,dV) # run solver
print('r='+str(r))
vf=calcVF(F,r,X,V0,dV)
print('vfrac='+str(vf))

## write setup

text=""
#line = lambda x,y,z: '\n\tins box\n\t\tp1='+str(x-dx/2)+','+str(y-dy/2)+','+str(z-dz/2)+'\n\t\tp2='+str(x+dx/2)+','+str(y+dy/2)+','+str(z+dz/2)+'\n\tendi'
xi = np.linspace(X0[0],X0[1],int(1+(X0[1]-X0[0])/dx))
yi = np.linspace(Y0[0],Y0[1],int(1+(Y0[1]-Y0[0])/dy))
zi = np.linspace(Z0[0],Z0[1],int(1+(Z0[1]-Z0[0])/dz))

# x,y,z=np.meshgrid(xi,yi,zi)
X=list(product(xi,yi,zi))
count=0
line = lambda x,y,z: '\n\tins box\n\t\tp1='+str(x)+','+str(y)+','+str(z)+'\n\t\tp2='+str(x+dx)+','+str(y+dy)+','+str(z+dz)+'\n\tendi'
for xx in X:
	x,y,z=xx
	X,Y,Z=T(x,y,z)
	#print(F(X,Y,Z))
	if abs(F(X,Y,Z))<r:
		text+=line(x,y,z)
		count+=1
print("N_Boxels: "+str(count))
with open(fname_in,'r') as IN:
	s=IN.read()
	with open(fname_out,'w') as OUT:
		OUT.write(s.replace(tag,text))#.replace('{L}',str(L)))
