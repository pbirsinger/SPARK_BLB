"""
Functional linear regression test.
"""

import blb
import numpy
from run_test import percent_error

def create_gaussian_data( num_observations, dimension, mapping=None, seed=None ):
    results = numpy.zeros( num_observations*(dimension+1), dtype='float32' )
    if seed is not None:
	numpy.random.seed( seed )
    if mapping is None:
	mapping = numpy.random.randn( dimension )
    for i in xrange(num_observations):
	x = numpy.random.randn( dimension )
	y = reduce( float.__add__, map( float.__mul__, x, mapping ) )
	numpy.put( results, i*dimension, y )
	numpy.put( results, [ j+i*dimension  for j in xrange(1,dimension+1) ], x )
    return results, mapping

class LinRegBLB( blb.BLB ):
    def compute_estimate( self, sample ):
	linreg()

    def reduce_bootstraps( self, sample ):
	mean()
	
    def average( self, sample ):
	mean()



if __name__ == '__main__':
    data, mapping = create_gaussian_data( 200, 50 )
    print 'mapping is, ', mapping
    cblb = LinRegBLB(dimension=51, num_subsamples=30)
    estimated_mapping = cblb.run( data )
    assert len(mapping) == len(estimated_mapping), "Mapping dimension are inequal!" 
    for i in xrange( len( mapping ) ):
	print 'checking covariate ', i
	p = percent_error( mapping[i], estimated_mapping[i] )
	if p > 0.05:
	    print 'Error: covariate %d is %f, estimated as %f, percent error %f' % (i, mapping[i], estimated_mapping[i], p )
    print 'done'