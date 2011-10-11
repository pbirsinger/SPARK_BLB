import csv
import math
import numpy

from blb import *

def mean(sample):
  return sum(sample) * 1.0 / len(sample)

def stddev(sample):
  mn = mean(sample)
  return math.sqrt(sum([(x - mn) * (x - mn) for x in sample]) * 1.0 / (len(sample) - 1))

def variance(sample):
  mn = mean(sample)
  return stddev(sample) ** 2

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
    def __init__(self):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean'
        BLB.__init__(self)

class CMeanMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'mean'
    self.average = 'mean'
    BLB.__init__(self)

class CSDMean_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'mean'
    self.reduce_bootstraps = 'stdev'
    self.average = 'mean'
    BLB.__init__(self)

class CSDSD_BLB(BLB):
  def __init__(self):
    self.compute_estimate = 'stdev'
    self.reduce_bootstraps = 'stdev'
    self.average = 'mean'
    BLB.__init__(self)

class MeanVariance_BLB(BLB):
  def compute_estimate(self, sample):
    return variance(sample)

  def reduce_bootstraps(self, sample):
    return mean(sample)

  def average(self, sample):
    return mean(sample)

class CMeanVariance_BLB(BLB):
  def __init__(self):
    self.compute_estimate = c_variance
    self.reduce_bootstraps = 'mean'
    self.average = 'mean'
    BLB.__init__(self)

f = open('trainingData.csv')
trainingdata = csv.reader(f)
data = []
first = True
for event in trainingdata:
    if first:
        first = False
        continue
    else:
        data.append(float(event[1]))
 
    
blb = MeanMean_BLB()
result = blb.run(data)
print ("Mean of Mean: ", result)

blb = SDMean_BLB()
result = blb.run(data)
print ("SD of Mean: ", result)

blb = MeanSD_BLB()
result = blb.run(data)
print ("Mean of SD: ", result)

blb = SDSD_BLB()
result = blb.run(data)
print ("SD of SD: ", result)

blb = CMeanSD_BLB()
result = blb.run(data)
print ("Mean of SD... in C: ", result)

blb = CMeanMean_BLB()
result = blb.run(data)
print ("Mean of Mean... in C: ", result)

blb = CSDSD_BLB()
result = blb.run(data)
print ("SD of SD... in C: ", result)

blb = MeanVariance_BLB()
result = blb.run(data)
print ("Mean of Variance: ", result)

blb = CMeanVariance_BLB()
result = blb.run(data)
print ("Mean of Variance... in C: ", result)
