"""Main class representing the BLB algorithm.

"""

import random

class BLB:
    known_reducers= ['mean', 'stdev']
    def __init__(self, num_subsamples=100, num_bootstraps=25, 
                 subsample_len_exp=0.5):

        self.pure_python = False
        for method in [ 'compute_estimate', 'reduce_bootstraps', 'average' ]:
            if not hasattr(self, method):
                raise ValueError("No method %s defined" % method)
            else:
                method_f = getattr(self, method)
                if hasattr( method_f, '__call__' ):
                    print "Warning: using pure-python mode for %s" % method
                    self.pure_python = True
                        
        if self.pure_python and str in map( type, [ self.compute_estimate, self.average, 
                            self.reduce_bootstraps ] ):
            raise ValueError('Cannot use specialised methods in pure python mode')
        self.num_subsamples = num_subsamples
        self.num_bootstraps = num_bootstraps
        self.subsample_len_exp = subsample_len_exp

    def run(self, data):
        pure_python = self.pure_python
        if pure_python:
            subsample_estimates = []
            for i in range(self.num_subsamples):
                subsample = self.__subsample(data, self.subsample_len_exp)
                bootstrap_estimates = [] 
                for j in range(self.num_bootstraps):
                    bootstrap = self.__bootstrap(subsample)
                    estimate = self.compute_estimate(bootstrap)
                    bootstrap_estimates.append(estimate)
                subsample_estimates.append(self.reduce_bootstraps(bootstrap_estimates))
            return self.average(subsample_estimates)
        else:
            sub_n = int( pow( len(data), self.subsample_len_exp ) )
            import asp.codegen.templating.template as template
            blb_template = template.Template(filename="templates/blb_template.mako", disable_unicode=True)
            impl_template = template.Template(filename="templates/blb_impl.mako", disable_unicode=True)
            rendered = blb_template.render( sub_n=sub_n, n_data=len(data), n_subsamples=self.num_subsamples,
                                            n_bootstraps=self.num_bootstraps, compute_estimate=self.compute_estimate,
                                            reduce_bootstraps=self.reduce_bootstraps, average=self.average )

            
            impl_args ={}
            if self.compute_estimate in BLB.known_reducers:
                impl_args['use_classifier'] = self.compute_estimate
            else:
                impl_args['classifier'] = self.compute_estimate

            if self.reduce_bootstraps in BLB.known_reducers:
                impl_args['use_bootstrap_reducer'] = self.reduce_bootstraps
            else:
                impl_args['bootstrap_reducer'] = self.reduce_bootstraps

            if self.average in BLB.known_reducers:
                impl_args['use_subsample_reducer'] = self.average
            else:
                impl_args['subsample_reducer'] = self.average

            rendered_impl = impl_template.render( **impl_args )

            import asp.jit.asp_module as asp_module
            mod = asp_module.ASPModule()
            mod.add_function('compute_estimate', rendered_impl)
            mod.add_header('stdlib.h')
            mod.add_header('math.h')
            mod.add_header('time.h')
            mod.add_header('numpy/ndarrayobject.h')
            mod.add_function("compute_blb", rendered)
            mod.add_to_init('import_array();')
            f = open('/tmp/test.cpp', 'w+')
            f.write( str(mod.backends['c++'].module.generate()) )
            f.close()
            return mod.compute_blb(data)

    def __subsample(self, data, subsample_len_exp):
        subsample_len = int(len(data) ** subsample_len_exp)
        subsample = random.sample(data, subsample_len)
        return subsample

    def __bootstrap(self, data):
        bootstrap = [random.choice(data) for i in range(len(data))]
        return bootstrap

    # These three methods are to be implemented by subclasses
    #def compute_estimate(self, sample):
    #    '''The actual statistic being computed. E.g. mean, standard deviation,
    #    etc. This is run on just the bootstrapped samples in the inner loops.
    #    '''
    #    TypeError('compute_estimate not defined')
    #
    #def reduce_bootstraps(self, sample):
    #    TypeError('reduce_bootstraps not defined')
    #
    #def average(self, sample):
    #    TypeError('average not defined')

