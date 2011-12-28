#!/usr/bin/env python

import math, random, unittest2 as unittest
from blb import *

DIM = 10

def mean(sample, dim, verify=False ):
    if dim == 1:
        avg = 0.0
        for i in xrange(len(sample)):
            avg += sample[i]
        if verify:
	    print "Mean computed (python): " + str(avg) + " for dim " + str(dim) + " mean is " + str(avg / len(sample))
        return avg / len(sample)
    else:
        print "Calculating mean of " + str(sample)
        avg = [0.0]*dim
        print str(avg)
        print "number of items is " + str(len(sample)) + " and dim is " + str(dim)
        for i in xrange(len(sample)):
            avg[i % dim] += sample[i]
        print str(avg)
	for j in xrange(dim):
	    avg[j] /= (len(sample) / dim)
	if verify:
	    print "Mean computed (python): " + str(avg) + " for dim " + str(dim) + " mean is " + str(avg)
        return avg

def stddev(sample, dim):
    if dim == 1:
        avg = mean(sample, dim);
        stdev = 0.0
        for i in xrange(len(sample)):
            stdev += (sample[i] - avg) ** 2
        stdev /= len(sample)
        return math.sqrt(stdev)
    else:
        dev = [0.0] * dim
        mean_vec = mean(sample, dim, verify=True)
        for vec in [ sample[i*dim:(i+1)*dim] for i in xrange(len(sample)/dim)] :
            for j in xrange(dim):
	        dev[j] += (vec[j] - mean_vec[j])**2
	for j in xrange(dim):
	    dev[j] /= len(sample) / dim
	    dev[j] = math.sqrt(dev[j])
        return dev
    
def variance(sample):
  mn = mean(sample, DIM)
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

def flatten(sample):
  flat = []
  for item in sample:
    flat.extend(item)
  return flat

def norm(mean_vec):
  norm = reduce( float.__add__, [ x*x for x in mean_vec ] )
  return math.sqrt(norm)


class MeanMean_BLB(BLB):
  def compute_estimate(self, sample):
    return mean(sample, DIM)

  def reduce_bootstraps(self, sample):
    sample = flatten(sample)
    return mean(sample, DIM)

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM, True)
    #print "mean_vec is " + str(mean_vec)
    #print "Norm of mean_vec is " + str(norm(mean_vec))
    x = norm( mean_vec )
    print 'Python Mean of Mean (Norm): %f' % x
    return x

class SDMean_BLB(BLB):
  def compute_estimate(self, sample):
    return mean(sample, DIM)

  def reduce_bootstraps(self, sample):
    print stddev([ vec[0] for vec in sample ], 1)
    return stddev(flatten(sample), DIM)

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM)
    return norm(mean_vec)

class MeanSD_BLB(BLB):
  def compute_estimate(self, sample):
    est = stddev(sample, DIM)
    print "Estimate computed for sample size " + str(len(sample)) + " is " + str(est)
    return est

  def reduce_bootstraps(self, sample):
    sample = flatten(sample)
    est = mean(sample, DIM)
    print "Reduced bootstrap for sample size " + str(len(sample)) + " : " + str(est)
    return est

  def average(self, sample):
    sample = flatten(sample)
    mean_vec = mean(sample, DIM)
    est = norm(mean_vec)
    print "Average for sample size " + str(len(sample)) + " is " + str(est)
    return est

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
    def __init__( self, **kargs ):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean_norm'
        BLB.__init__( self, **kargs )

class CMeanMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'mean'
    self.average= 'mean_norm'
    BLB.__init__( self, dimension=DIM )

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
  if measured == true:
      return 0
  return abs( float( measured - true )/true )

def generate_data(dim):
  data = [0.0]*10000
  for i in xrange(10000 / dim):
    if (i % 2 == 0):
      val = 1000.0
    else:
      val = 2000.0
    for j in xrange(dim):
      data[i*dim + j] = val
  return data

BLB_FAIL_MSG = 'Test Case: %s; Python result: %f; C result: %f; Percent Error: %f'
data = generate_data(DIM)


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
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Mean of Mean', py_res, c_res, p) )
    print "C Mean of Mean test passed: C result: " + str(c_res) + " and Python result: " + str(py_res)

    # mean of stdev
#    print "Computing MEAN of STDEV"
    py_blb = MeanSD_BLB()
    c_blb = CMeanSD_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
#    if (p<= self.threshold):
#      print "Mean of SD passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Mean of SD', py_res, c_res, p) )

    # stdev of mean
#    print "Computing STDEV of MEAN"
    py_blb = SDMean_BLB()
    c_blb = CSDMean_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
#    if (p <= self.threshold):
#      print "SD of mean passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('SD of Mean', py_res, c_res, p) )

    # stdev of stdev
#    print "Computing STDEV of STDEV"
    py_blb = SDSD_BLB()
    c_blb = CSDSD_BLB()
    py_res = py_blb.run(data)
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
#    if (p <= self.threshold):
#      print "SD of SD passed"
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('SD of SD', py_res, c_res, p) )


#  def test_numnpy(self):
#    """
#    Ensure the numpy c version works.
#    """
#    #data = [ random.random() for i in xrange(10000) ]
#    print "Computing SD of MEAN"
#    py_blb = SDMean_BLB()
#    py_res = py_blb.run(data)
#    c_blb = CSDMean_BLB()
#    npy_data = numpy.array(data, dtype='float32')
#    c_res = c_blb.run(npy_data)
#    p = percent_error(py_res, c_res)
#    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Numpy Test', py_res, c_res, p) )


  def test_custom(self):
    """
    Ensure custom c functions as classifiers work.
    (Assuming they are correctly written!)
    """
    """
    py_blb = MeanVariance_BLB()
    py_res = py_blb.run(data)
    c_blb = CMeanVariance_BLB()
    c_res = c_blb.run(data)
    p = percent_error(py_res, c_res)
    self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Custom Classifier Test', py_res, c_res, p) )
    """


if __name__ == '__main__':
  unittest.main()
