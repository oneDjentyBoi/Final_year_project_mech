"""
All quantities are in SI units.
This is the code for the temperature profile of a cylinder.
Thermal conductivity is assumed to be constant here.
Temp at r=0 > Temp at r=R
e = sub-diagonal vector
f = diagonal vector
g = super_diagonal vector 
"""
#############################################################
Tr=300 # temp at r=R
q_v=1000000 # heat gen
k=100 # thermal conductivity
r=0.04 # radius
n=50 # number of divisions
Tf = Tr + (q_v*r*r)/4 # temp at r=0; solved from analytical formula.
##############################################################
import matplotlib.pyplot as plt
del_r=r/n 
e=[0]*(n-2)
f=[0]*(n-2)
g=[0]*(n-3)
x=[0]*(n-2)
y=[0]*(n-2)
y[0]=-((2*del_r-del_r)*Tr + (q_v/k)*2*del_r*del_r*del_r)
y[n-3]=-((2*del_r*(n-2)+del_r)*Tf + (q_v/k)*2*del_r*(n-2)*del_r*del_r) 
for i in range(1,n-3):
	y[i] = -(q_v/k)*2*del_r*del_r*del_r*(i+1)
for i in range(1,n-2):
	e[i]=2*del_r*(i+1)-del_r
for i in range(0,n-2):
	f[i]=-4*del_r*(i+1)
for i in range(0,n-3):
	g[i]=2*del_r*(i+1) + del_r
for i in range(1,n-2):
	factor=e[i]/f[i-1]
	f[i]=f[i]-factor*g[i-1]
	y[i]=y[i]-factor*y[i-1]
x[n-3]=y[n-3]/f[n-3]
for i in range(n-4,-1,-1):
	x[i] = (y[i] - g[i]*x[i+1])/f[i]
# print(e)
# print(f)
# print(g)
# print(y)
# print(x)
y1=[0]*n
y1[0]=Tf
y1[n-1]=Tr
x1=[0]*n
for i in range(n-3,-1,-1):
	y1[n-i-2]=x[i]
for i in range(1,n):
	x1[i]= del_r*i
plt.plot(x1,y1)
plt.show()