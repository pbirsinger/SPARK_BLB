import blb, math
import unittest2 as unittest
import numpy, random

class MeanSDBLB( blb.BLB ):

    def compute_estimate( X ):
	mean = vector( dim( X, 1), dtype(X) )
	square_sum = vector( dim( X, 1), dtype(X) )
	for x in X:
	    mean += x
	    square_sum += x**2
	mean /= len( X )
	square_sum /= len( X )
	return  sqrt(square_sum - mean**2)


    def reduce_bootstraps( bootstraps ):
	mean = vector( dim(bootstraps, 1), dtype(bootstraps) )
	for bootstrap in bootstraps:
	    mean += bootstrap
	mean /= len(bootstraps)
	return mean
	
    def average( subsamples ):
	mean = vector( dim(subsamples, 1), dtype(subsamples) )
	for x in subsamples:
	    mean += x
	mean /= len(subsamples)
	return mean    


class LinregBLB( blb.BLB ):

    def compute_estimate( X, y ):
	XtX = matrix( dim( X, 1 ), dim( X, 1 ), dtype(X) )
	Xy = vector( dim( X, 1 ), dtype(X) )
	for x in X:
	    XtX += outer_product( x, x )
	for i in dim( y, 0 ):
	    Xy += y[i]*X[i]
	return mv_solve( XtX, Xy )

    def reduce_bootstraps( bootstraps ):
	mean = vector( dim( bootstraps, 1 ), dtype(bootstraps) )
	for x in bootstraps:
	    mean += x
	mean /= len( bootstraps )
	return mean    

    def average(  subsamples ):
	mean = vector( dim(subsamples, 1), dtype( subsamples) )
	for x in subsamples:
	    mean += x
	mean /= len( subsamples )
	return mean    

def percent_error( actual, measured ):
    return abs( ( measured - actual )/actual )

def generate_sd_data():
    N = 100	
    d = 100
    sigma = math.sqrt( (N**2 - 1) / 12 )
    choices = range(1, N+1)
    data = []	
    while len( choices) > 0:
	x = random.choice( choices )
	choices.remove( x )
	data.append( [ float(x) ]* d )
    arr = numpy.array( data )
    arr.shape = (N, d)
    return arr, numpy.array( [sigma]*d  )

def generate_linreg_data():
    num_observations = 3000
    dimension = 90 
    X = numpy.empty( num_observations*dimension, dtype='float64' )
    y_eps = 0.0001
    y = numpy.empty( num_observations, dtype='float64' )
    mapping = numpy.random.randn( dimension )
    for i in xrange(num_observations):
	x = numpy.random.randn( dimension )
	numpy.put( y, i, (numpy.random.rand() * y_eps) + reduce( float.__add__, map( float.__mul__, x, mapping ) ) )
	numpy.put( X, xrange( i*dimension, (i+1)*dimension) , x )
    X.shape = ( num_observations, dimension )
    return X, y, mapping

class GenerationTest( unittest.TestCase ):

    def test_MeanSD( self ):
	tester = MeanSDBLB(subsample_len_exp=0.7, num_subsamples=64)
	data, sigma = generate_sd_data()
	sigma_est = tester.run( data )
#	print 'sigma', sigma
#	print 'sigma_est', sigma_est  
	err = percent_error( sigma, sigma_est )
#	print 'err', err
	self.assertTrue( ( err < 0.1 ).all(), 'Error out of bounds' )

    def _test_LinReg( self ):
	tester = LinregBLB(subsample_len_exp=0.7)
	X, y, mapping = generate_linreg_data()
	mapping_est = tester.run( X, y )
#	print 'mapping', mapping
#	print 'mapping_est', mapping_est
	x = [ a < 0.1 for a in percent_error( mapping, mapping_est ) ]
	self.assertTrue( reduce( lambda a,b: a and b, x ), "LinReg Generatrion test: error out of bounds" )

if __name__ == '__main__':
    unittest.main()