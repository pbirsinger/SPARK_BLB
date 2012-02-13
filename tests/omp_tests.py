import unittest2 as unittest
from run_test import CMeanSD_BLB, percent_error
from bootstrap import Bootstrap
import numpy, random

def generate_data( n, d ):
    #generate 8 numbers, weights adding up to n
    m = min(8, n)
    weights = numpy.random.multinomial( n, [1.0/m]*m )
    vals = range( 100, 1000, 100 ) 
    vals = [ random.choice( vals ) for i in range( m ) ]   
    data = numpy.empty( n*d, dtype='float64' )
    offset = 0
    mean = numpy.zeros( d )
    for i in range( m ):
	step = weights[i]*d
	numpy.put( data, xrange( offset, offset+step ), [ vals[i] ]*step )
	mean += weights[i]*numpy.array( [vals[i]]* d )
	offset += step
    for i in range( len(data )):
	if data[i] not in vals:
	    raise ValueError( "%d not in vals, index %d" % ( data[i], i ) )
    mean /= n
    var = numpy.zeros( d )
    for i in range( m ):
	vec = numpy.array( [vals[i]]*d )
	residual = vec - mean
	var += weights[i]*residual*residual
    var /= n
    return data, mean, numpy.sqrt( var ) , vals, weights

def norm( vec ):
    return numpy.linalg.norm( vec ) 

class  CMeanSD_BOOTSTRAP( Bootstrap ):
    def compute_estimate( self ):
	std()
    def average( self ):
	mean()

FAIL_MSG = 'Test "%s" failed: actual %s, measured %s, percent error %s'
class OMPTest( unittest.TestCase ):
    def setUp(self):
	self.threshold = .05
	data, mean, std, vals, weights = generate_data( 1250, 8 )
	self.data = data
	self.mean = mean
	self.std = std
	self.vals = vals
	self.weights = weights

    def test_blb(self):
	c_blb = CMeanSD_BLB(dimension = 8)
	two_blb = CMeanSD_BLB(with_openMP=True, dimension = 8)
	two_blb.omp_n_threads=2
	four_blb = CMeanSD_BLB(with_openMP=True, dimension = 8)
	four_blb.omp_n_threads=4
	res_serial = c_blb.run(self.data)
	res_two = two_blb.run(self.data)
	res_four = four_blb.run(self.data)
	actual = norm( self.std )
	p = percent_error( res_serial, actual)
	self.assertTrue(  p <= self.threshold , msg = FAIL_MSG % ('serial', actual, res_serial, p ) )
	p = percent_error(res_two, actual)
	self.assertTrue( p <= self.threshold , msg = FAIL_MSG % ('two threads', actual, res_two, p) )
	p = percent_error(res_four, actual)
	self.assertTrue( p <= self.threshold , msg = FAIL_MSG % ('four threads', actual, res_four, p) )

    def test_bootstrap( self ):
	c_bootstrap = CMeanSD_BOOTSTRAP(dimension = 8)
	two_bootstrap = CMeanSD_BOOTSTRAP(with_openMP=True, dimension = 8)
	two_bootstrap.omp_n_threads=2
	four_bootstrap = CMeanSD_BOOTSTRAP(with_openMP=True, dimension = 8)
	four_bootstrap.omp_n_threads=4
	res_serial = c_bootstrap.run(self.data)
	res_two = two_bootstrap.run(self.data)
	res_four = four_bootstrap.run(self.data)
	actual =  self.std 
	p = percent_error( res_serial, actual)
	self.assertTrue( (p <= self.threshold).all() , msg = FAIL_MSG % ('serial', actual, res_serial, p ) )
	p = percent_error(res_two, actual)
	self.assertTrue( (p <= self.threshold).all() , msg = FAIL_MSG % ('two threads', actual, res_two, p) )
	p = percent_error(res_four, actual)
	self.assertTrue( (p <= self.threshold).all() , msg = FAIL_MSG % ('four threads', actual, res_four, p) )

if __name__ == '__main__':
    unittest.main()