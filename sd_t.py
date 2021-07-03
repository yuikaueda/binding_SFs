import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt('nor_n10_r0rt1000000.dat')
data2 = np.loadtxt('nor_n1000_r0rt1000000.dat')

x1 = data1[:,0]
y1 = data1[:,1]

x2 = data2[:,0]
y2 = data2[:,1]

fig, ax = plt.subplots(1, 1)

ax.plot(x1, y1, 'o', markersize=1, c='black', label = r'$N=10$')
ax.plot(x2, y2, 'o', markersize=1, c='red', label = r'$N=1000$')

plt.xlabel("Time", fontsize = 18)
plt.ylabel("S.D.", fontsize = 18)

plt.xlim(0,1)

ax.legend(loc='best')
fig.savefig("nor_n10n1000_r1000000.png")
plt.show()
