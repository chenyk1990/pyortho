#  DEMO script (python version) for calculating local orthogonalization and its application in random noise attenuation 
#  
#  Copyright (C) 2022 Yangkang Chen
#  
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details: http://www.gnu.org/licenses/
#  
#  Reference:   1. Random noise attenuation using local signal-and-noise orthogonalization
#               Chen and Fomel, 2015, Geophysics
#               2. Ground-Roll Noise Attenuation Using a Simple and Effective Approach Based on 
#               Local Band-Limited Orthogonalization, Chen et al., 2015, IEEE Geoscience and Remote Sensing Letters
#               3. Iterative deblending with multiple constraints based on shaping regularization,
#               Chen, 2015, IEEE Geoscience and Remote Sensing Letters
#               4. Orthogonalized morphological reconstruction for weak signal detection in micro-seismic monitoring:
#               Methodology, Huang et al., 2018, GJI
#               5. Surface-related multiple leakage extraction using local primary-and-multiple 
#               orthogonalization, Zhang et al., 2020, Geophysics
#               6. Non-stationary local signal-and-noise orthogonalization, Chen et al.,
#               2021, Geophysics
#               7. Local primary-and-multiple orthogonalization for leaked internal multiple crosstalk estimation and attenuation on full-wavefield migrated images
#               Zhang, et al., 2021, Geophysics

## generate synthetic data
#This synthetic data was used in Huang et al., 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
import numpy as np
import matplotlib.pyplot as plt

## generate the synthetic data
a1=np.zeros([300,20])
[n,m]=a1.shape
a3=np.zeros([300,20])
a4=np.zeros([300,20])

k=-1;
a=0.1;
b=1;
pi=np.pi

ts=np.arange(-0.055,0.055+0.002,0.002)
b1=np.zeros([len(ts)])
b2=np.zeros([len(ts)])
b3=np.zeros([len(ts)])
b4=np.zeros([len(ts)])

for t in ts:
    k=k+1;
    b1[k]=(1-2*(pi*30*t)*(pi*30*t))*np.exp(-(pi*30*t)*(pi*30*t));
    b2[k]=(1-2*(pi*40*t)*(pi*40*t))*np.exp(-(pi*40*t)*(pi*40*t));
    b3[k]=(1-2*(pi*40*t)*(pi*40*t))*np.exp(-(pi*40*t)*(pi*40*t));
    b4[k]=(1-2*(pi*30*t)*(pi*30*t))*np.exp(-(pi*30*t)*(pi*30*t));

t1=np.zeros([m],dtype='int')
t3=np.zeros([m],dtype='int')
t4=np.zeros([m],dtype='int')
for i in range(m):
  t1[i]=np.round(140);
  t3[i]=np.round(-6*i+180);
  t4[i]=np.round(6*i+10);
  a1[t1[i]:t1[i]+k+1,i]=b1; 
  a3[t3[i]:t3[i]+k+1,i]=b1; 
  a4[t4[i]:t4[i]+k+1,i]=b1; 

temp=a1[0:300,:]+a3[0:300,:]+a4[0:300,:];

shot=np.zeros([300,20,20])
for j in range(20):
    a4=np.zeros([300,20]);
    for i in range(m):
    	t4[i]=np.round(6*i+10+3*j); 
    	a4[t4[i]:t4[i]+k+1,i]=b1;
  
    	t1[i]=np.round(140-2*j);
    	a1[t1[i]:t1[i]+k+1,i]=b1;

    shot[:,:,j]=a1[0:300,:]+a3[0:300,:]+a4[0:300,:];

d0=shot

## add noise
[n1,n2,n3]=d0.shape
np.random.seed(201415)
n=0.1*np.random.randn(n1,n2,n3);
dn=d0+n;
print(np.std(dn))

from fxydmssa import fxydmssa
d1=fxydmssa(dn,0,120,0.004,3,3);	#DMSSA (when damping factor =1, there are heavy damages)
noi1=dn-d1;

## compare with matlab
import scipy
from scipy import io
datas = {"d0":d0,"dn": dn, "d1": d1, "noi1": noi1}
scipy.io.savemat("datas.mat", datas)

## plotting
fig = plt.figure(figsize=(5, 6))
fig.add_subplot(3, 1, 1)
plt.imshow(dn.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.2, 0.2))

fig.add_subplot(3, 1, 2)
plt.imshow(d1.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.2, 0.2))

fig.add_subplot(3, 1, 3)
plt.imshow(noi1.transpose(0,2,1).reshape(n1,n2*n3),cmap='jet',clim=(-0.2, 0.2))
plt.show()





