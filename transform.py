from params import labda,theta as th
import numpy as np

from math import cos as c
from math import sin as s

from random import random as rand
from math import pi

delta=np.array([rand() for _ in range(3)])

def T(x,y,z):
	r=[x,y,z]
	r=[(r[i]+delta[i])*labda[i] for i in range(3)]
	
	R1=[	[1,0,0],
		[0,c(th[0]),s(th[0])],
		[0,-s(th[0]),c(th[0])]
	]
	R2=[	[c(th[1]),0,s(th[1])],
		[0,1,0],
		[s(th[1]),0,c(th[1])]
	]
	R3=[	[c(th[2]),s(th[2]),0],
		[-s(th[2]),c(th[2]),0],
		[0,0,1]
	]
	R=np.matmul(np.matmul(R3,R2),R1)

	return list(np.matmul(R,r))
