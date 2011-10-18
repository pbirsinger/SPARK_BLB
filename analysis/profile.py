import blb, yep

class ProfBLB(blb.BLB):
    def __init__(self, **kwargs):
        self.compute_estimate = 'stdev'
        self.reduce_bootstraps = 'mean'
        self.average = 'mean'
        blb.BLB.__init__(self, **kwargs)

if __name__ == '__main__':
    data1 = range(10000)
    data2 = range(50000)
    data3 = range(100000)

    tester = ProfBLB()
    tester.run(data3)

    yep.start('cilk.prof')
    for i in xrange(500):
        tester.run(data3)
    yep.stop()
    
    
