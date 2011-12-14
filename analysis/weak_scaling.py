import numpy
import time
import blb
from cref import *


class FakeArray( object ):
    def __init__(self, size, gen=True):
	self.size = size
	if gen:
	    self.arr = generate_data((size,))
	else:
	    self.arr = numpy.empty((1,))
    def __len__(self):
	return self.size

    def __array__(self):
	return 	self.arr

if __name__ == '__main__':
    numpy.random.seed(0x0adb)
    for i in [1, 2, 4, 8, 16, 32, 64]:
	tester = CBLB( subsample_len_exp=0.7, num_subsamples=64, with_openMP=True)
	#force compile
	tester.omp_n_threads = i
	tester.compile_for([], key=(int(i*GB), numpy.ndarray))
	print 'Size %d generating data' % i
	data = generate_data((i*GB,))
	#time
	print 'Size %d about to time' % i
	start = time.time()
	tester.run(data)
	ellapsed = time.time()
	print 'Time taken for size %d was %f' % ( i, (ellapsed-start) )