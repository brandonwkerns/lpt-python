import numpy as np
from context import lpt
import matplotlib.pylab as plt
import datetime as dt

objdir='/home/orca/bkerns/public_html/realtime_mjo_tracking/lpt/data/tmpa/objects'
obj1 = 20190327000002
obj2 = 20190328000001

overlapping_points = lpt.helpers.calc_overlapping_points(obj1, obj2, objdir)
print(overlapping_points)

lon1, lat1, mask1 = lpt.helpers.get_lpo_mask(obj1, objdir)
lon2, lat2, mask2 = lpt.helpers.get_lpo_mask(obj2, objdir)

plt.close('all')
plt.contour(lon1, lat1, mask1, levels=(0.5,), colors='b')
plt.contour(lon2, lat2, mask2, levels=(0.5,), colors='r')
plt.text(0.1, 0.9, ('(N1(blue), N2(red), Noverlap) = '+str(overlapping_points)),transform=plt.gca().transAxes)
plt.savefig('test_lpo_overlap.png')
print('test_lpo_overlap.png')
