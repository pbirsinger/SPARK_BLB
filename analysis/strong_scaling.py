import blb
import numpy, time
from cref import *


if __name__ == '__main__':

    dsize = long(8*GB)
    sizes = [ 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64 ]
    testers = {}
    for i in sizes:
	tester = CBLB(num_subsamples=64, subsample_len_exp=0.7, with_openMP=True, dimension=(2**18))
	tester.omp_n_threads = i
	tester.compile_for( [], key=( dsize, numpy.ndarray ) )
	testers[ i ] = tester

    data = generate_data((dsize,))
    print 'Data unpacked.'
    for i in sizes:
        tester = testers[i]
	start = time.time()
	tester.run(data)
	ellapsed = time.time()
	print "Time ellapsed for %d threads: %f" % ( i, (ellapsed-start) )