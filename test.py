import numpy as np
from matplotlib import pyplot as plt

n=50# number of divs
x=np.zeros((n,n),dtype="float")

q_v=100
k=100

# x[0][0]=-1
# x[0][0]=1
# x[n-1][n-1]=1

r=1
del_r=r/n

for i in range(0,n):
	for j in range(0,n):
		if i==j:
			x[i][j]=-(4*del_r*(i+1))

for i in range(0,n):
	for j in range(0,n):
		if i-j==1:
			x[i][j]=(2*del_r*(i+1))-del_r

for i in range(0,n):
	for j in range(1,n):
		if j-i==1:
			x[i][j]=(2*del_r*(i+1))+del_r

x[n-1][n-2]=1
x[n-1][n-1]=-1
y=np.zeros((n,1),dtype="float")
y[0]=-((2*del_r)-del_r)*300
# y[n-1]=-((2*del_r*(n-2))+del_r)*200
for i in range(1,n-1):
	y[i]=-(q_v/k)*2*del_r*del_r*del_r*(i)
y[n-1]=-q_v/(4*k)


y[n-1]=-q_v/(4*k)
# y[n-1]=300

for i in range(1,n-1):
	y[i]=-i*2*q_v*del_r*del_r*del_r/k
op=np.linalg.inv(x)
t=np.dot(op,y)
print (x)
print(y)
_x=[]
for i in range(1,n+1):
	_x.append((n-i)*del_r)
t=list(t)
plt.plot(_x,t)
plt.show()
print(t)