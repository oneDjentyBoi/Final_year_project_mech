import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n=1500
r=1
del_r=r/(n-1)
Tr=300
q=10000
k=105.6
x=np.zeros((n,n),dtype="float")

x[0][0]=-1
for i in range(1,n):
	for j in range(0,n):
		if i==j:
			x[i][j]=-2
x[n-1][n-1]=1
for i in range(1,n-1):
	for j in range(0,n):
		if i-j==1:
			x[i][j]=1
for i in range(0,n-1):
	for j in range(0,n):
		if j-i==1:
			x[i][j]=1

y=[-q*del_r*del_r/k]*n

y[0]=y[0]/2
y[n-1]=Tr

x_inv=np.linalg.inv(x)
t=np.dot(x_inv,y)
print(t)

_y=[0]*n

for i in range(0,n):
	_y[i]=i*del_r

plt.plot(_y,t)
plt.show()

# print(x)
# print(y)