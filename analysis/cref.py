#Reference C time for weak & strong scaling analyses
import numpy
import time
import blb

class CBLB( blb.BLB ):
    def __init__(self, **kargs):
	self.compute_estimate = 'stdev'
	self.reduce_bootstraps = 'mean'
	self.average = 'mean'
	blb.BLB.__init__(self, **kargs)

def generate_data( size, dt='float32' ):
    if reduce( int.__mul__, size) == 0:
	return []
    try:
	numpy.dtype(dt)
    except TypeError:
	numpy.dtype('float32')
    return numpy.random.rand(*size)

#one gigabyte for a four-byte data type
GB = 25e7
 
if __name__ == '__main__':
    data = generate_data((GB,))
    print 'size of data: %f' % len(data)
    tester = CBLB(num_subsamples=64, subsample_len_exp=0.7)
    #fore compile
    tester.run(data)
    print 'compilation complete, about to time...'
    start = time.time()
    tester.run(data)
    ellapsed = time.time()
    print 'Serial time on one GB: %f' % (ellapsed-start)