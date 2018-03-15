import numpy as np
from utm_wgs84 import LatLon
import pandas as pd

print('Reading data file')
lld = pd.read_csv('B2.xyz', sep=';', usecols=[
                  0, 1, 2], dtype=np.float32).values


def ll_to_utm(ll):
    utm = LatLon(ll[0], ll[1]).to_utm()
    return np.array([utm.easting, utm.northing, ll[2]])


print('Converting to utm')
utmd = np.apply_along_axis(ll_to_utm, 1, lld)

print(utmd)
