from linreg_experiment import create_gaussian_data, LinRegBLB, LinRegBootstrap
import numpy, time
import blb, bootstrap

def percent_error( measured, actual ):
    return numpy.abs( (measured - actual) / actual )
GB = 2 ** 21
if __name__ == '__main__':
    """
    test_time = 0
    spec = ( GB*32*64, numpy.ndarray )
    bootstrapper = LinRegBootstrap( with_openMP=True, dimension=65 )
    bootstrapper.omp_n_threads = 32
    bootstrapper.compile_for( None, key=spec )
    blber = LinRegBLB( with_openMP=True, dimension=65, subsample_len_exp=0.7, num_subsamples=64, num_bootstraps=40 )
    blber.omp_n_threads = n
    blber.compile_for( None, key=spec )
    data, mapping = create_gaussian_data( 32*GB, 64 )
    start = time.time()
    bootstrapper.run( data )
    ellapsed = time.time()
    test_time = ellapsed - start
    print 'Bootstrap size %d threads %d time %f' % ( 32, 32, test_time )
    start = time.time()
    blber.run( data )
    ellapsed = time.time()
    test_time = ellapsed - start
    print 'BLB size %d threads %d time %f' % (32, 32, test_time )
    """
    for size in [ 32 ]:
	data, mapping = create_gaussian_data( size*GB, 64 )
	print 'Size %d, data generated' % size
	for n in [  32  ]:
	    test_time = 0
	    tester = LinRegBootstrap( with_openMP=True, dimension=65 )
	    tester.omp_n_threads = n
	    start = time.time()
#	    estimated_mapping = tester.run( data )
	    ellapsed = time.time()
	    test_time = ellapsed - start
	    print 'Bootstrap size %d threads %d time %f' % (size, n, test_time)	
	    tester = LinRegBLB( with_openMP=True, dimension=65, subsample_len_exp=0.7, num_subsamples=64, num_bootstraps=40 )
	    tester.omp_n_threads = n
	    start = time.time()
	    estimated_mapping = tester.run( data )
	    ellapsed = time.time()
	    test_time = ellapsed - start
	    print 'BLB size %d threads %d time %f' % (size, n, test_time)
