from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('nor_n10_r0rt100.dat')
data2 = np.loadtxt('nor_n10_r0rt1000.dat')
data3 = np.loadtxt('nor_n10_r0rt10000.dat')
data4 = np.loadtxt('nor_n10_r0rt100000.dat')
data5 = np.loadtxt('nor_n1000_r0rt100.dat')
data6 = np.loadtxt('nor_n1000_r0rt1000.dat')
data7 = np.loadtxt('nor_n1000_r0rt10000.dat')
data8 = np.loadtxt('nor_n1000_r0rt100000.dat')

def fit (x,T,A):
    return np.exp(-x/T)+A

x1 = data1[:,0]
y1 = data1[:,1]

x2 = data2[:,0]
y2 = data2[:,1]

x3 = data3[:,0]
y3 = data3[:,1]

x4 = data4[:,0]
y4 = data4[:,1]

x5 = data5[:,0]
y5 = data5[:,1]

x6 = data6[:,0]
y6 = data6[:,1]

x7 = data7[:,0]
y7 = data7[:,1]

x8 = data8[:,0]
y8 = data8[:,1]

param1, cov1 = curve_fit(fit, x1, y1)
param2, cov2 = curve_fit(fit, x2, y2)
param3, cov3 = curve_fit(fit, x3, y3)
param4, cov4 = curve_fit(fit, x4, y4)
param5, cov5 = curve_fit(fit, x5, y5)
param6, cov6 = curve_fit(fit, x6, y6)
param7, cov7 = curve_fit(fit, x7, y7)
param8, cov8 = curve_fit(fit, x8, y8)

#print(param1[0],param2[0],param3[0],param4[0])
y_T1 = [param1[0],param2[0],param3[0],param4[0]]
y_T2 = [param5[0],param6[0],param7[0],param8[0]]

x = [2,3,4,5]

fig, ax = plt.subplots(1, 1)
fig.subplots_adjust(bottom=0.15)
ax.plot(x, y_T1, 'o', c='red', label = 'N=10')
ax.plot(x, y_T2, 'o', c='blue', label = 'N=1000')
plt.xlabel(r'log(|r|)', fontsize = 18)
plt.ylabel(r'parameter T', fontsize = 18)
ax.legend(loc='best')
fig.savefig("n10n1000_r_param_T.png")
plt.show()

