import sys
import numpy
import pandas


ts = pandas.read_table(sys.argv[1])
ts['wc'] = ts.z - ts.zbottom
ts = ts.sort('id').reset_index()
print ts.head()
ts['diffx'] = 0
ts['diffy'] = 0
ts['diffx'][1:] = numpy.diff(ts.x)
ts['diffy'][1:] = numpy.diff(ts.y)

cols = ['Time', 'x','y','z','state','wc','diffx','diffy']
print ts[ts.id==1].sort('Time')[ts.Time >= 731586][ts.Time < 731587][cols[:-1]]