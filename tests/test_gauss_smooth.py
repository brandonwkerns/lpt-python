import numpy as np
from context import lpt
import matplotlib.pylab as plt

kernel = lpt.helpers.gauss_smooth_kernel(121,121,20,20)

#print(kernel)
print(np.sum(kernel))

fig=plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax1.pcolormesh(kernel)
ax2 = fig.add_subplot(1,2,2)
ax2.plot(kernel[:,61])
plt.show()


bogus_data = np.zeros([100,100])
bogus_data[40:61, 40:61] = 1
fig=plt.figure()
ax1 = fig.add_subplot(1,2,1)
H1 = ax1.pcolormesh(bogus_data, vmin=0, vmax=1)
plt.colorbar(H1)

bogus_data_filtered = lpt.helpers.gauss_smooth(bogus_data,19,19,3,3)

ax2 = fig.add_subplot(1,2,2)
H2 = ax2.pcolormesh(bogus_data_filtered, vmin=0, vmax=1)
plt.colorbar(H2)

plt.show()
