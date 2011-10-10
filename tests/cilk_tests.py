from blb import *
import math
import unittest2 as unittest
import time

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

def percent_error( s1, s2 ):
    return abs( (s1 - s2)/float(s2) )

class CilkTest(unittest.TestCase):
    
    def setUp(self):
        self.py_blb = PyMeanStdev_BLB()
        self.c_blb = CMeanStdev_BLB()
        self.cilk_blb = CilkMeanStdev_BLB()

    def test_error(self):
        data = range(10000)
       
        py_res = self.py_blb.run(data)
        c_res = self.c_blb.run(data)
        cilk_res = self.cilk_blb.run(data)
        
        self.assertTrue( .02 > percent_error( py_res, c_res ) )
        self.assertTrue( .02 > percent_error( c_res, cilk_res ) )
        self.assertTrue( .02 > percent_error( py_res, cilk_res ) )

    def test_speedup(self):
        data = range( 10 * 1000 )
        n_iters = 100
        #force compilation
        self.c_blb.run(data)
        self.cilk_blb.run(data)

        start = time.clock()
        for i in xrange(n_iters):
            self.py_blb.run(data)
        py_time = (time.clock() - start)/float(n_iters)

        start = time.clock()
        for i in xrange(n_iters):
            self.c_blb.run(data)
        c_time = (time.clock() - start)/float(n_iters)

        start = time.clock()
        for i in xrange(n_iters):
            self.cilk_blb.run(data)
        cilk_time = (time.clock() - start)/float(n_iters)

        print 'c speedup: %fx' % (py_time/c_time)
        print 'cilk speedup: %fx' % (py_time/cilk_time)
if __name__ == '__main__':
    unittest.main()
