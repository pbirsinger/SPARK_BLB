#!/usr/bin/env python

import math, random, unittest2 as unittest
from blb import *

def mean(sample):
  return sum(sample)*1.0/len(sample)

def stddev(sample):
  mn = mean(sample)
  return math.sqrt(sum([(x-mn)*(x-mn) for x in sample])*1.0/(len(sample)-1))

def variance(sample):
  mn = mean(sample)
  return stddev(sample)**2

c_variance = """
    float mean = 0.0;
    for( unsigned int i=0; i< size; i++ ){
        mean += weights[i]*data[i];
    }
    mean /= DATA_SIZE;
    float var = 0.0;
    for( unsigned int i=0; i<size; i++ ){
        float datum = (data[i] - mean);
        var += weights[i]*datum*datum;
    }
    return var / DATA_SIZE;
"""

class MeanMean_BLB(BLB):
  def compute_estimate(self, sample):
    return mean(sample)

  def reduce_bootstraps(self, sample):
    return mean(sample)

  def average(self, sample):
    return mean(sample)

class SDMean_BLB(BLB):
  def compute_estimate(self, sample):
    return mean(sample)

  def reduce_bootstraps(self, sample):
    return stddev(sample)

  def average(self, sample):
    return mean(sample)

class MeanSD_BLB(BLB):
  def compute_estimate(self, sample):
    return stddev(sample)

  def reduce_bootstraps(self, sample):
    return mean(sample)

  def average(self, sample):
    return mean(sample)

class SDSD_BLB(BLB):
  def compute_estimate(self, sample):
    return stddev(sample)

  def reduce_bootstraps(self, sample):
    return stddev(sample)

  def average(self, sample):
    return mean(sample)

class CMeanSD_BLB(BLB):
    def __init__( self, **kargs ):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean'
        BLB.__init__( self, **kargs )

class CMeanMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'mean'
    self.average= 'mean'
    BLB.__init__( self )

class CSDMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'stdev'
    self.average = 'mean'
    BLB.__init__( self )

class CSDSD_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'stdev'
    self.reduce_bootstraps = 'stdev'
    self.average = 'mean'
    BLB.__init__( self )

class MeanVariance_BLB(BLB):
  def compute_estimate(self, sample):
    return variance( sample )

  def reduce_bootstraps( self, sample ):
    return mean(sample)

  def average(self, sample):
    return mean(sample)

class CMeanVariance_BLB(BLB):
  def __init__(self):
    self.compute_estimate = c_variance
    self.reduce_bootstraps = 'mean'
    self.average = 'mean'
    BLB.__init__( self )

def percent_error( true, measured ):
  return abs( float( measured - true )/true )

# mean = 50.0, variance = 312.5, stdev = 17.67
data  = [ 25.0 for i in xrange(2500) ] + [ 50.0 for i in xrange(5000) ] + [ 75.0 for i in xrange(2500) ]

BLB_FAIL_MSG = 'Test Case: %s; Python result: %f; C result: %f; Percent Error: %f'
class BLBTest(unittest.TestCase):
  def setUp(self):
    # max allowable percent error to be considered correct.
    self.threshold = .1
  
  def test_error(self):
    """
    Ensure that basic computations are accurate.
    """
    # mean of mean
    py_blb = MeanMean_BLB()
    c_blb = CMeanMean_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Mean of Mean', py_res, c_res, p) )
    # mean of stdev
    py_blb = MeanSD_BLB()
    c_blb = CMeanSD_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Mean of SD', py_res, c_res, p) )

  def test_numnpy(self):
    """
    Ensure the numpy c version works.
    """
    data = [ random.random() for i in xrange(10000) ]
    py_blb = MeanMean_BLB()
    py_res = py_blb.run(data)
    c_blb = CMeanMean_BLB()
    npy_data = numpy.array(data, dtype='float32')
    c_res = c_blb.run(npy_data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Numpy Test', py_res, c_res, p) )


  def test_custom(self):
    """
    Ensure custom c functions as classifiers work.
    (Assuming they are correctly written!)
    """
    py_blb = MeanVariance_BLB()
    py_res = py_blb.run(data)
    c_blb = CMeanVariance_BLB()
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Custom Classifier Test', py_res, c_res, p) )


if __name__ == '__main__':
  unittest.main()
