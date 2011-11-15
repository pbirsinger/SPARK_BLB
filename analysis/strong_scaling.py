import blb
import numpy, time
from cref import *


if __name__ == '__main__':

    data = generate_data((GB,))
    print 'Data unpacked.'
    for i in [ 1, 2, 4, 8, 16, 32, 64 ]:
        tester = CBLB(num_subsamples=64, subsample_len_exp=0.7, with_openMP=True)
	tester.omp_n_threads=i
	# force compile
	tester.run(data)
	print 'compilation complete. About to time...'
	start = time.time()
	tester.run(data)
	ellapsed = time.time()
	print "Time ellapsed for %d threads: %f" % ( i, (ellapsed-start) )