import numpy, time
import blb, bootstrap

def create_gaussian_data( num_observations, dimension, mapping=None, seed=None ):
    results = numpy.empty( num_observations*(dimension+1), dtype='float64' )
    print 'memory allocated'
    y_eps = 0.01
    if seed is not None:
	numpy.random.seed( seed )
    if mapping is None:
	mapping = numpy.random.randn( dimension )
    checkins = [ int(num_observations*0.1*r) for r in range( 10 ) ]
    offset = dimension + 1
    for i in xrange(num_observations):
	if i in checkins:
	    print '%f percent' % (100*float( i ) / num_observations)
	x = numpy.random.randn( dimension )
	noise = y_eps * numpy.random.rand()
	y = noise + reduce( float.__add__, map( float.__mul__, x, mapping ) )
	numpy.put( results, i*offset, y )
	numpy.put( results, xrange( i*offset+1, (i+1)*offset), x )
    return results, mapping

class LinRegBLB( blb.BLB ):
    def compute_estimate( self, sample ):
	linreg()

    def reduce_bootstraps( self, sample ):
	std()

    def average( self, sample ):
	mean()

class LinRegBootstrap( bootstrap.Bootstrap ):
    def compute_estimate( self, sample ):
	linreg()

    def average( self, sample ):
	std()
   
class Regressor( blb.BLB ):
    def compute_estimate( self, sample ):
	linreg()

    def reduce_bootstraps( self, sample ):
	mean()

    def average( self, sample ):
	mean()

if __name__ == '__main__':
    #generate an estimation of std using many data sets from one mapping
    dim = 64
    n = 10 ** 6
    mapping = numpy.random.randn( dim )
    for i in [ 1, 2, 4, 8, 16, 32, 64 ]:
	estimates = numpy.empty( 2000, dim )
        regressor = Regressor( with_openMP = True, dimension = 65 )
        regressor.omp_n_threads = i
	start = time.time()  
        for i in xrange( 2000 ):
	    data, blah = create_gaussian_data( n, dim, mapping = mapping )
	    estimates[i] =  regressor.run( data )
        print 'Calculating std_true...'
        std_true = numpy.std( estimates, axis = 0 )
	ellapsed = time.time()
	print '%d threads, manual calculation, time = %f' % (i, ellapsed - start)  
        #estimate std with BLB, bootstrap using various thread counts
        data, blah = create_gaussian_data( n, dim, mapping = mapping )
        print 'data formed'
        bootstrapper = LinRegBootstrap(with_openMP=True, dimension = 65)
        blber = LinRegBLB(with_openMP=True, dimension = 65)
	bootstrapper.omp_n_threads = n
	blber.omp_n_threads = n
	start = time.time()
	bootstrapper.run( data )
	ellapsed = time.time()
	print '%d threads, bootstrap time %f' % (i, ellapsed - start )
	blber.run(data)
	b_ellapsed = time.time()
	print '%d threads, blb time %f' % (i, b_ellapsed - ellapsed )
