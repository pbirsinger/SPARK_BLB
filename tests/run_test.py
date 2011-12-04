#!/usr/bin/env python

import math, random, unittest2 as unittest
from blb import *

DIM = 8

def mean(sample, dim):
  if dim == 1:
    print sample
    avg = 0.0
    for i in xrange(len(sample)):
      avg += sample[i]
#    print "Mean computed (python): " + str(avg)
    return avg / len(sample)

#make sure data is homogeneous
#  for i in xrange(dim):
#    print str(sample[i]) + " "
#  print '\n'
  avg = [0.0]*dim
  for i in xrange(len(sample)):
    avg[i % dim] += sample[i]
  for j in xrange(dim):
    avg[j] /= (len(sample) / dim)
#  print "Mean computed (python): " + str(avg)
  return avg

def stddev(sample, dim):
  if dim == 1:
    avg = mean(sample, dim);
    stdev = 0.0
    for i in xrange(len(sample)):
      stdev += (sample[i] - avg) ** 2
    stdev /= len(sample)
#    print "STDEV computed (python): " + str(math.sqrt(stdev))
    return math.sqrt(stdev)
  dev = [0.0] * dim
  mean_vec = mean(sample, dim)
  for vec in [ sample[i*dim:(i+1)*dim] for i in xrange(len(sample)/dim)] :
    for j in xrange(dim):
      dev[j] += (vec[j] - mean_vec[j])**2
  for j in xrange(dim):
    dev[j] /= len(sample) / dim
    dev[j] = math.sqrt(dev[j])
  return dev
    
def squaredNorm(vec, dim):
  norm = 0.0
  for i in xrange(dim):
    norm += vec[i] * vec[i]
  return norm

def variance(sample):
  mn = mean(sample, DIM)
  return stddev(sample)**2

c_variance = """
    float mean = 0.0;
    for( int i=0; i< size; i++ ){
        mean += data[ indicies[i] ];
    }
    mean /= size;
    float var = 0.0;
    for( int i=0; i<size; i++ ){
        float datum = data[ indicies[i] ];
        var += pow( datum - mean, 2 );
    }
    return var / size;
"""

def flatten(sample):
  flat = []
  for item in sample:
    flat.extend(item)
  return flat

def norm(mean_vec):
  norm = 0.0
  for i in xrange(len(mean_vec)):
    norm += mean_vec[i] ** 2
  return math.sqrt(norm)


class MeanMean_BLB(BLB):
  def compute_estimate(self, sample):
    return mean(sample, DIM)

  def reduce_bootstraps(self, sample):
    sample = flatten(sample)
    return mean(sample, DIM)

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM)
    return norm(mean_vec)

class SDMean_BLB(BLB):
  def compute_estimate(self, sample):
    return mean(sample, DIM)

  def reduce_bootstraps(self, sample):
    return stddev(flatten(sample), DIM)

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM)
    print "Average vector for python "
    print mean_vec
    print "\n"
    return norm(mean_vec)

class MeanSD_BLB(BLB):
  def compute_estimate(self, sample):
    return stddev(sample, DIM)

  def reduce_bootstraps(self, sample):
    sample = flatten(sample)
    return mean(sample, DIM)

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM)
    return norm(mean_vec)

class SDSD_BLB(BLB):
  def compute_estimate(self, sample):
    return stddev(sample, DIM)

  def reduce_bootstraps(self, sample):
    return stddev(flatten(sample), DIM)

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM)
    return norm(mean_vec)

class CMeanSD_BLB(BLB):
    def __init__( self ):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean_norm'
        BLB.__init__( self )

class CMeanMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'mean'
    self.average= 'mean_norm'
    BLB.__init__( self )

class CSDMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'stdev'
    self.average = 'mean_norm'
    BLB.__init__( self )

class CSDSD_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'stdev'
    self.reduce_bootstraps = 'stdev'
    self.average = 'mean_norm'
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
    self.average = 'mean_norm'
    BLB.__init__( self )

def percent_error( true, measured ):
  return abs( float( measured - true )/true )

def generate_data(dim):
  data = [0.0]*10000
  for i in xrange(10000 / dim):
    if (i % 2 == 0):
      val = 2000.0
    else:
      val = 1000.0
    for j in xrange(dim):
      data[i*dim + j] = val
  return data

BLB_FAIL_MSG = 'Test Case: %s; Python result: %f; C result: %f; Percent Error: %f'
data = generate_data(8)


class BLBTest(unittest.TestCase):
  def setUp(self):
    # max allowable percent error to be considered correct.
    self.threshold = .1
  
  def test_error(self):
    """
    Ensure that basic computations are accurate.
    """
    ## biased to prevent small-argument errors
    #data = [ 100 + random.random() for i in xrange(10000) ]

    # mean of mean
    py_blb = MeanMean_BLB()
    c_blb = CMeanMean_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    if (p<= self.threshold): 
      print "Mean of Mean passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Mean of Mean', py_res, c_res, p) )

    # mean of stdev
    print "Computing MEAN of STDEV"
    py_blb = MeanSD_BLB()
    c_blb = CMeanSD_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    if (p<= self.threshold):
      print "Mean of SD passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Mean of SD', py_res, c_res, p) )

    # stdev of mean
    print "Computing STDEV of MEAN"
    py_blb = SDMean_BLB()
    c_blb = CSDMean_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    if (p <= self.threshold):
      print "SD of mean passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('SD of Mean', py_res, c_res, p) )

    # stdev of stdev
    print "Computing STDEV of STDEV"
    py_blb = SDSD_BLB()
    c_blb = CSDSD_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    if (p <= self.threshold):
      print "SD of SD passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('SD of SD', py_res, c_res, p) )


  def test_numnpy(self):
    """
    Ensure the numpy c version works.
    """
    #data = [ random.random() for i in xrange(10000) ]
    print "Computing SD of MEAN"
    py_blb = SDMean_BLB()
    py_res = py_blb.run(data)
    c_blb = CSDMean_BLB()
    npy_data = numpy.array(data, dtype='float32')
    c_res = c_blb.run(npy_data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Numpy Test', py_res, c_res, p) )


  def test_custom(self):
    """
    Ensure custom c functions as classifiers work.
    (Assuming they are correctly written!)
    """
    """
    #data = [ random.random() for i in xrange(10000) ]
    py_blb = MeanVariance_BLB()
    py_res = py_blb.run(data)
    c_blb = CMeanVariance_BLB()
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Custom Classifier Test', py_res, c_res, p) )
    """


if __name__ == '__main__':
  unittest.main()
