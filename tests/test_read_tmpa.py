import numpy as np
from context import lpt
import matplotlib.pylab as plt

fn = '/home/orca/data/satellite/trmm_global_rainfall/2011/11/20111101/3B42.20111101.06.7.HDF'
DATA = lpt.readdata.read_tmpa_hdf(fn)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
DATA['precip'][DATA['precip'] < 0.001] = np.NaN
H1 = ax1.pcolormesh(DATA['lon'], DATA['lat'], DATA['precip'], vmin=0.0, vmax=5.0, cmap='jet')
plt.colorbar(H1)
plt.show()


fnrt = '/home/orca/data/satellite/trmm_global_rainfall/rt/2019/01/20190128/3B42RT.2019012800.7.bin'
RT = lpt.readdata.read_tmpa_rt_bin(fnrt)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
RT['precip'][RT['precip'] < 0.001] = np.NaN
H2 = ax1.pcolormesh(RT['lon'], RT['lat'], RT['precip'], vmin=0.0, vmax=5.0, cmap='jet')
plt.colorbar(H2)
plt.show()
