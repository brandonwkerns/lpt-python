import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt

"""
fn = '/home/orca/data/satellite/trmm_global_rainfall/2011/11/20111101/3B42.20111101.06.7.HDF'
DATA = lpt.readdata.read_tmpa_hdf(fn)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
DATA['precip'][DATA['precip'] < 0.001] = np.NaN
H1 = ax1.pcolormesh(DATA['lon'], DATA['lat'], DATA['precip'], vmin=0.0, vmax=5.0, cmap='jet')
plt.colorbar(H1)
plt.show()
"""

#fnrt = '/home/orca/data/satellite/cmorph/rt/2019/02/20190217/CMORPH_V0.x_RT_8km-30min_2019021700'
#RT = lpt.readdata.read_cmorph_rt_bin(fnrt)
RT = lpt.readdata.read_cmorph_at_datetime(dt.datetime(2019,2,17,0,0,0), verbose=True)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
RT['precip'][RT['precip'] < 0.001] = np.NaN
H2 = ax1.pcolormesh(RT['lon'], RT['lat'], RT['precip'][0,:,:], vmin=0.0, vmax=5.0, cmap='jet')
plt.colorbar(H2)
plt.show()
