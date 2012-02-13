import blb
import numpy, time

class CBLB ( blb.BLB ):
    def __init__(self, **kargs):
	self.compute_estimate = 'stdev'
	self.reduce_bootstraps = 'mean'
	self.average = 'mean'
	blb.BLB.__init__(self, **kargs)


def generate_data():
    numpy.random.seed(0xadb)
    return numpy.random.rand(*( 25 * 10**7,))

if __name__ == '__main__':

    data = generate_data()
    f.close()
    print 'Data unpacked.'
    niters = 5
    bstart = time.time()
    for i in xrange(niters):
	print 'Bias print'
    bias = time.time() - bstart
    for i in [ 2, 4, 8, 16, 32, 64 ]:
        tester = CBLB(subsample_len_exp=0.7)
	tmp_data = data[:len(data)/i]
	# force compile
	tester.run(tmp_data)
	print 'compilation complete. About to time...'
	start = time.time()
	for j in xrange(niters):
	    print 'Timing - iteration %d' % j
            tester.run(tmp_data)
	ellapsed = ((time.time() - start)/float(niters)) - bias
	print "Time ellapsed for %d threads: %f" % ( i, ellapsed )