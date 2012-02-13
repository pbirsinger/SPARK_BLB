"""
Functional linear regression test.
"""

import blb
import numpy
from run_test import percent_error

def create_gaussian_data( num_observations, dimension, mapping=None, seed=None ):
    results = numpy.empty( num_observations*(dimension+1), dtype='float64' )
    y_eps = 0.01
    if seed is not None:
	numpy.random.seed( seed )
    if mapping is None:
	mapping = numpy.random.randn( dimension )
    offset = dimension + 1
    for i in xrange(num_observations):
	x = numpy.abs( numpy.random.randn( dimension ) )
	#y = (numpy.random.rand() * y_eps) + reduce( float.__add__, map( float.__mul__, x, mapping ) )
	y = reduce( float.__add__, x*mapping )
	numpy.put( results, i*offset, y )
	numpy.put( results, xrange( i*offset+1, (i+1)*offset) , x )
    return results, mapping

class LinRegBLB( blb.BLB ):
    def compute_estimate( self, sample ):
	linreg()

    def reduce_bootstraps( self, sample ):
	mean()
	
    def average( self, sample ):
	mean()



if __name__ == '__main__':
    data, mapping = create_gaussian_data( 3000, 50 )
    print 'mapping is, ', mapping
    cblb = LinRegBLB(dimension=51, num_subsamples=30)
    estimated_mapping = cblb.run( data )
    assert len(mapping) == len(estimated_mapping), "Mapping dimension are inequal!" 
    for i in xrange( len( mapping ) ):
	print 'checking covariate ', i
	p = percent_error( mapping[i], estimated_mapping[i] )
	if p > 0.05:
	    print 'Error: covariate %d is %f, estimated as %f, percent error %f' % (i, mapping[i], estimated_mapping[i], p )
    print estimated_mapping
    print 'done'