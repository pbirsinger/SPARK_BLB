from run_test import BLB, generate_data
import unittest2 as unittest
class NoBLB( BLB ):
    def compute_estimate( self, sample ):
	noop()
    
    def reduce_bootstraps( self, sample ):
	noop()

    def average( self, sample ):
	noop() 

class NoOpTest( unittest.TestCase ):
    def test( self ):
	data = generate_data(8)
	test_blb = NoBLB()
	res = test_blb.run(data)
	self.assertEqual( res, 0 )


if __name__ == '__main__':
    unittest.main()