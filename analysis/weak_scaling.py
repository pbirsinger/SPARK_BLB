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
    data = []
    for i in [ 32 ]:
	print 'Size %d generating data' % i
	#data = numpy.append(data, generate_data((long((i/2)*GB),)))
	print 'Size %d, number of floats = %d' % (i, len(data))
	tester = CBLB( subsample_len_exp=0.7, num_subsamples=64, with_openMP=True)
	tester.omp_n_threads = i
	#force compile
	tester.compile_for([], key=(int(32*GB), numpy.ndarray))
	data = generate_data((32*GB,))
	#time
	print 'Size %d about to time' % i
	start = time.time()
	tester.run(data)
	ellapsed = time.time()
	print 'Time taken for size %d was %f' % ( i, (ellapsed-start) )