"""
All quantities are in SI units.
This is the code for the temperature profile of a cylinder.
Thermal conductivity is assumed to be constant here.
Temp at r=0 > Temp at r=R
b = sub-diagonal vector
a = diagonal vector
c = super_diagonal vector 
"""
import numpy as np
from matplotlib import pyplot as plt

a=[]
b=[]
c=[]

n=50# number of nodes
r=1
del_r=r/(n-1)
Tr=300 #temperature at r=R
q=1000000
k=205

h=20000
T_inf=273

switch=1 # 0 for inverse and 1 for TDMA
neumann=0
dirichlet=1

def TDMAsolver(a, b, c, d):
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]
    # print(bc)
    # print(dc)
    # print(cc)
    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

x=np.zeros((n,n),dtype="float")

x[0][0]=-1
a.append(-1)#TDMA f
x[0][1]=1
c.append(1) # TDMA f

for i in range(1,n-1):
	for j in range(0,n):
		if i==j:
			x[i][j] = -4*del_r*i
			a.append(-4*del_r*i)#TDMA

x[n-1][n-1]=1
a.append(1.0)#TDMA l

for i in range(1,n-1):
	for j in range(0,n):
		if i-j==1:
			x[i][j]= 2*i*del_r-del_r
			b.append(2*i*del_r-del_r)#TDMA

b.append(0.0)#TDMA l

for i in range(1,n-1):
	for j in range(0,n):
		if j-i==1:
			x[i][j]=2*del_r*i+del_r
			c.append(2*i*del_r+del_r)#TDMA

y=[0]*n

for i in range(1,n-1):
	y[i]=-2*del_r*del_r*del_r*i*q/k

y[0]=-q*del_r*del_r/k
y[0]=y[0]/4
y[n-1]=Tr

_y=[0]*n

for i in range(0,n):
	_y[i]=i*del_r

if neumann==1:
	a[-1]=(2*del_r*h*(2*r+del_r)/k) - 4*r
	b[-1]=4*r
	x[n-1][n-1]=(2*del_r*h*(2*r+del_r)/k) - 4*r
	x[n-1][n-2]=4*r
	y[-1]=-2*r*del_r*del_r*q/k + 2*(2*r+del_r)*del_r*h*T_inf/k
	if switch ==1:
		f=TDMAsolver(b,a,c,y)
		print(f)
		plt.plot(_y,f)
		plt.show()
		print(x)
	elif switch == 0:
		x_inv=np.linalg.inv(x)
		t=np.dot(x_inv,y)
		print(t)
		plt.plot(_y,t)
		plt.show()
if dirichlet==1:
	if switch==0:
		x_inv=np.linalg.inv(x)
		t=np.dot(x_inv,y)
		print(t)
		plt.plot(_y,t)
		plt.show()

	elif switch==1:
		f=TDMAsolver(b,a,c,y)
		print(f)
		plt.plot(_y,f)
		plt.show()

# print(x)
# print(y)
