import unittest2 as unittest
from run_test import CMeanSD_BLB, BLB_FAIL_MSG, percent_error, generate_data
import numpy

class OMPTest( unittest.TestCase ):
    def setUp(self):
	self.threshold = .1

    def test(self):
	ndatas = generate_data(8)
	c_blb= CMeanSD_BLB()
	two_blb = CMeanSD_BLB(with_openMP=True)
	two_blb.omp_n_threads=2
	four_blb = CMeanSD_BLB(with_openMP=True)
	four_blb.omp_n_threads=4
	res_serial = c_blb.run(ndatas)
	res_two = two_blb.run(ndatas)
	res_four = four_blb.run(ndatas)
	p = percent_error(res_serial, res_two)
	self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Serial/two', res_serial, res_two, p) )
	p = percent_error(res_serial, res_four)
	self.assertTrue( p <= self.threshold, msg = BLB_FAIL_MSG % ('Serial/four', res_serial, res_four, p) )
	
if __name__ == '__main__':
	unittest.main()