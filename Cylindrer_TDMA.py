"""
All quantities are in SI units.
This is the code for the temperature profile of a cylinder.
Thermal conductivity is assumed to be constant here.
Temp at r=0 > Temp at r=R
b = sub-diagonal vector
a = diagonal vector
c = super_diagonal vector 
"""
#############################################################
Tr=100 # temp at r=R
q_v=100 # heat gen
k=164 # thermal conductivity
r=4 # radius
n=50 # number of nodes
##############################################################
import numpy as np
import pandas as pd
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
import matplotlib.pyplot as plt

del_r=r/n

m=n-1 #number of divisions.

a=[]
b=[]
c=[]

# b=np.array([0.004, 0.012, 0.02, 0.028, 0.036000000000000004, 0.044, 0.052000000000000005, 0.06])
# a=np.array([1, -0.016, -0.032, -0.048, -0.064, -0.08, -0.096, -0.112, -0.06])
# c=np.array([0, 0.012, 0.02, 0.028, 0.036000000000000004, 0.044, 0.052000000000000005, 0.06])
# y=np.array([300, 0, 0, 0, 0, 0, 0, 0, 0])
# t=TDMAsolver(b,a,c,y)
# print(t)

if q_v>0:
	y=[0]*n
	y[0]=Tr
	a.append(1)
	for i in range(1,n-1): 
		a.append(-4*del_r*i)
	a.append(-1)
	for i in range(1,n-1):
		b.append((2*del_r*i)-del_r)
	b.append(1)
	c.append(0)
	for i in range(1,n-1):
		c.append((2*del_r*i)+del_r)

	y[n-1] = -q_v/(4*k)

	for i in range(1,n-1):
		y[i]=-(q_v*del_r*del_r*del_r*2*i)/k

	t=TDMAsolver(b,a,c,y)
	print(t)

	r_x=[0]*n

	for i in range(0,n):
		r_x[i]=(n-i)*del_r
	plt.plot(r_x,t)
	plt.show()
	_x=np.array(r_x)
	_y=np.array(t)
	_y=_y[::-1]
	arr= {'x':_x,'y':_y}
	df=pd.DataFrame(arr)
	df.to_excel("plot.xlsx",index=True)
else:
	a.append(1)
	for i in range(1,n-2):
		a.append(-4*del_r*i)
	a.append(-4*del_r*(m-1)+2*del_r*(m-1)+del_r)
	for i in range(1,n-1):
		b.append((2*del_r*i)-del_r)
	c.append(0)
	for i in range(1,n-2):
		c.append((2*del_r*i)+del_r)

	y=[0]*(n-1)
	y[0]=Tr

	t=TDMAsolver(b,a,c,y)
	print(t)

	r_x=[0]*(n-1)

	for i in range(0,n-1):
		r_x[i]=(n-1-i)*del_r
	plt.plot(r_x,t)
	plt.show()

# print((a))
# print((b))
# print((c))
# print((y))

