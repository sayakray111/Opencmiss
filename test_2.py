import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv
X = []
Cp = []
Cx = []
Cy = []
Cz = []
C = []
#f = open("/hpc/sray036/VirtualPregnancy/repro-examples/uteroplacental/darcy-models/02-darcy-ell-full/output/StaticDarcy.part0.exnode",'r')
g = open("StaticDarcy.part0.exnode")
currentline  = 1
Pdrop = 100*0.1
Pbase = 0
porosity = 0.6
vol_frac = 1.00 - porosity

piperadius = .5
density = 1060.0
pipelength = 3.0

delP = (Pdrop-Pbase)/(pipelength*1e-2)

d = .500
perm = d ** 2. * (1. - vol_frac) ** 3. / (180. * vol_frac ** 2.)
perm = 0.0052*1e-4
viscosity = 0.112

perm_vis = perm/viscosity
print('Permeability is =',perm)

gamma = 1/(1+(2.5*(1-porosity)))
print('Beta = ',gamma)
print('Permeability = ',perm_vis)
Darcy_Number = perm/(piperadius*piperadius*1e-4)
omega = np.sqrt(gamma/Darcy_Number)

print('Omega = ',omega)
v_z = (perm_vis)*(delP)*100

print(v_z)
print('Darcy Number = ',Darcy_Number)
C1 = np.sqrt(gamma/perm)
ST = g.readlines()
r = []
f = []
fd = []
sum1 = 0.0
avg = 0.0
print(omega)
Re = density*v_z*(2*piperadius*1e-4)/viscosity
print('Reynolds number is = ',Re)
with open("StaticDarcy.part0.exnode") as file:
    for i,l in enumerate(file):
        pass
N = i+1
count = int((N-35)/15)
z_lim = round(float(ST[23+(count-1)*15].strip()),2)
print(z_lim)
Sum1 = 0.0
Sum2 = 0.0
Sum3 = 0.0
Nodes = 0
for n in range(1,count):
    #x_val = float(ST[9*n+6].strip())
    x_val = round(float(ST[21+(n-1)*15].strip()),2)
    #y_val = float(ST[9*n+7].strip())
    y_val = round(float(ST[22+(n-1)*15].strip()),2)
    #z_val = float(ST[9*n+8].strip())
    z_val = round(float(ST[23+(n-1)*15].strip()),2)
    #print(z_val)
    if(z_val == z_lim):
        f2 = round(float(ST[26+(n-1)*15].strip()),2)
        r1 = np.sqrt((x_val**2)+(y_val**2))
        dim_r1 = r1/(piperadius)
        ratio = iv(0,((omega)*dim_r1))/iv(0,(omega))
        f1 = perm_vis*delP*(1-ratio)*100
        r.append(r1)
        f.append(f1)
        fd.append(f2)
        Sum1 = Sum1+f2
        Sum2 = Sum2+f1
        error = (f1-f2)/(f1+1e-40)
        Sum3 = Sum3+abs(error)
        Nodes+=1

dict1 = {r[i]: fd[i] for i in range(0, len(fd))}
sorted_dic = sorted(dict1.items(), key=lambda x: x[0])
print(sorted_dic)
x1 = []
x2 = []
fd1 = []
Sum_Err = 0
for i in sorted_dic:
    x1.append(i[0])
    x2.append(i[1])
    ratio = (iv(0, (C1 * i[0] * 1e-2))) / (iv(0, (C1 * piperadius * 1e-2)))
    f11 = perm_vis * delP * (1 - ratio) * 100
    fd1.append(f11)
    err1 = abs(i[1]-f11)/(f11+1e-40)
    Sum_Err = Sum_Err+err1

avg_error = Sum3/Nodes
plt.xlabel('Radius(cm)')
plt.ylabel('Velocity(cm/sec)')
plt.plot(x1, fd1,'r')
plt.plot(x1, x2,'b')
plt.show()


avg1 = Sum3/Nodes
print('Average value (numerical) = ',Sum1/Nodes)
print('Average value (analytical) = ',Sum2/Nodes)
l1 = Sum1/Nodes
l2 = Sum2/Nodes
print('the average error is = ',(l1-l2)/l1*100,' ',avg1)
avg = 0
for k1 in range(0,len(x1)):
    if(x1[k1]==min(x1)):
        val = x2[k1]

print('the maximum value of velocity is  = ',val)
print('The analytical maximum velocity is = ',max(fd1))
# M = 6250 N = 598.54
# M = 9900 N = 599.07
# M = 10800 N = 599.07
# M = 14700 N = 599.44