import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('nor_n10_r0rt100000.dat')
data2 = np.loadtxt('nor_n100_r0rt100000.dat')
data3 = np.loadtxt('nor_n1000_r0rt100000.dat')
data4 = np.loadtxt('nor_n10000_r0rt100000.dat')

x1 = data1[:,0]
y1 = data1[:,1]

x2 = data2[:,0]
y2 = data2[:,1]

x3 = data3[:,0]
y3 = data3[:,1]

x4 = data4[:,0]
y4 = data4[:,1]

fig, ax = plt.subplots(1, 1)

ax.plot(x1, y1, 'o', markersize=1, c='black', label = r'$N=10$')
ax.plot(x2, y2, 'o', markersize=1, c='red', label = r'$N=100$')
ax.plot(x3, y3, 'o', markersize=1, c='blue', label = r'$N=1000$')
ax.plot(x4, y4, 'o', markersize=1, c='green', label = r'$N=10000$')

plt.xlabel("Time", fontsize = 18)
plt.ylabel("S.D.", fontsize = 18)

#plt.xlim(0,1)

ax.legend(loc='best')
fig.savefig("tmax_nor_n_all_r100000.png")
plt.show()
