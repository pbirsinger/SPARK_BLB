from blb import *
import math
import time, csv

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
    for i in xrange(1):
        impl.run(data)
    return (time.clock() - start )/float(1)

if __name__ == '__main__':

    py_blb = PyMeanStdev_BLB()
    c_blb = CMeanStdev_BLB()
    cilk_blb = CilkMeanStdev_BLB()
    omp_blb = CMeanStdev_BLB(with_openMP=True)

    py_simple_time = time_impl(py_blb, simple_data, 10)
    c_simple_time = time_impl(c_blb, simple_data, 10)
    cilk_simple_time = time_impl(cilk_blb, simple_data, 10)
    omp_simple_time = time_impl( omp_blb, simple_data, 10)
    
    py_time = time_impl(py_blb, data, 10)
    c_time = time_impl(c_blb, data, 10)
    cilk_time = time_impl(cilk_blb, data, 10)
    omp_time = time_impl( omp_blb, data, 10)

    print 'speedup on simple data:'
    print 'c speedup: %fx' % (py_simple_time/c_simple_time)
    print 'cilk speedup: %fx' % (py_simple_time/cilk_simple_time)
    print 'omp speedup: %fx' % (py_simple_time/omp_simple_time)

    print 'speedup on real data:'
    print 'c speedup: %fx' % (py_time/c_time)
    print 'cilk speedup: %fx' % (py_time/cilk_time)
    print 'omp speedup: %fx' % (py_time/omp_time)
    
