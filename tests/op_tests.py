import blb, math
import unittest2 as unittest
import numpy, random

class TestBLBBase( blb.BLB ):
    def reduce_bootstraps( bootstraps ):
	mean = vector( dim( bootstraps, 1 ), dtype( bootstraps ) )
	for bootstrap in bootstraps:
	    mean += bootstrap
	mean /= len(bootstraps)
	return mean

    def average( subsamples ):
	mean = vector( dim( subsamples, 1), dtype( subsamples ) )
	for subsample in subsamples:
	    mean += subsample
	return mean / len(subsamples)


class OpsTester( unittest.TestCase ):
    
    def setUp( self ):
	self.data = numpy.ones( (100, 100) )
	self.scalar_data = numpy.ones( (100*100,) )
	
    def test_add( self ):
	sumres = numpy.array( [100.0]*100 )
 
	class VectorAdd( TestBLBBase ):
	    def compute_estimate( X ):
		sum = vector( dim( X, 1 ), dtype( X ) )
		for x in X:
		    sum = sum + x
		return sum

	#tester = VectorAdd()
	#res = tester.run( self.data )
	#print 'add res', res
	#self.assertTrue( numpy.equal(res, sumres).all(), "VectorAdd failed!" )
 
	class VectorAugAdd( TestBLBBase ):
	    def compute_estimate( X ):
		sum = vector( dim( X, 1 ), dtype( X ) )
		for x in X:
		    sum += x
		return sum

	tester = VectorAugAdd(with_openMP = True)
	res = tester.run( self.data )
	print 'augadd res', res
	self.assertTrue( numpy.equal(res, sumres).all(), "VectorAugAdd failed!" )

    def test_sub( self ):
	sumres = numpy.array( [-100.0]*100 )
	
	class VectorSub( TestBLBBase ):
	    def compute_estimate( X ):
		sum = vector( dim( X, 1 ), dtype( X ) )
		for x in X:
		    sum = sum - x
		return sum
	#tester = VectorSub()
	#res = tester.run( self.data )
	#print 'sub res', res
	#self.assertTrue( numpy.equal( res, sumres ).all(), "VectorSub failed!" )

	class VectorAugSub( TestBLBBase ):
	    def compute_estimate( X ):
		sum = vector( dim( X, 1 ), dtype( X ) )
		for x in X:
		    sum -= x
		return sum
	tester = VectorAugSub(with_openMP=True)
	res = tester.run( self.data )
	print 'augsub res', res
	self.assertTrue( numpy.equal( res, sumres ).all(), "VectorAugSub failed!" )

if __name__ == '__main__':
    unittest.main()