from blb import *
import math
import time, csv
import random

def mean(sample):
    return sum(sample)*1.0/len(sample)

def stddev(sample):
    mn = mean(sample)
    return math.sqrt(sum([ (x-mn)*(x-mn) for x in sample ] )*1.0/(len(sample)-1))


class PyMeanStdev_BLB(BLB):
    def compute_estimate(self, sample):
        return stddev(sample)
    
    def reduce_bootstraps(self, sample):
        return mean(sample)

    def average(self, sample):
        return mean(sample)

class CMeanStdev_BLB(BLB):
    def __init__(self, **kargs):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean'
        BLB.__init__(self, **kargs)

class CilkMeanStdev_BLB(BLB):
    def __init__(self, **kargs):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean'
        kargs['with_cilk'] = True
        BLB.__init__(self, **kargs)

simple_data = range( 10000 )

f= open('trainingData.csv')
trainingdata = csv.reader(f)
data = []
first = True
for event in trainingdata:
    if first:
        first = False
        continue
    else:
        data.append(float(event[1]))


def time_impl(impl, data, n_iters):
    #force compile before timing
    impl.run(data)
    #time it
    start = time.clock()
    for i in xrange(n_iters):
        impl.run(data)
    return (time.clock() - start )/float(n_iters)

def make_data( n ):
    return [ random.random() for i in xrange(n) ]

if __name__ == '__main__':
    py_blb = PyMeanStdev_BLB()
    openMP_blb = CMeanStdev_BLB(with_openMP=True)
    n_iters = 100
    start = time.clock()
    for i in xrange(n_iters):
        pass
    bias = (time.clock() - start)/float(n_iters)
    
    for size in [ 4, 10, 16, 25, 50, 64, 100, 256 ]:
        
    
        datas = make_data( size * 1000 )
        cdatas = numpy.array( datas, dtype='float32' )
        print 'Running pure python version with size %d' % ( size * 1000 )
        start = time.clock()
        for i in xrange( n_iters ):
            py_blb.run( datas )
        py_time = (time.clock() - start)/float(n_iters) - bias
#        print 'py_time was %f' % py_time
        openMP_blb.run( cdatas )
        print 'Running OpenMP version with size %d' % (size * 1000)
        start = time.clock()
        for i in xrange( n_iters ):
            openMP_blb.run(cdatas)
        openMP_time = (time.clock() - start)/float(n_iters) - bias
#        print 'openMP time was %f' % openMP_time
        print 'OpenMP test completed for size %d. Speedup was %fx' % ((size * 1000), (py_time/openMP_time))

    print 'all done!'
