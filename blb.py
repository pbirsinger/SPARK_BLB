"""
Main class representing the BLB algorithm.

"""

import random
import numpy

class BLB:
    known_reducers= ['mean', 'stdev']
    def __init__(self, num_subsamples=100, num_bootstraps=25, 
                 subsample_len_exp=0.5, with_cilk=False, with_openMP=False):

        self.with_cilk=with_cilk
	self.with_openMP = with_openMP
        self.pure_python = False
        for method in [ 'compute_estimate', 'reduce_bootstraps', 'average' ]:
            if not hasattr(self, method):
                raise ValueError("No method %s defined" % method)
            else:
                method_f = getattr(self, method)
                if hasattr( method_f, '__call__' ):
                    self.pure_python = True
                        
        if self.pure_python and str in map( type, [ self.compute_estimate, self.average, 
                            self.reduce_bootstraps ] ):
            raise ValueError('Cannot use specialised methods in pure python mode')
        self.num_subsamples = num_subsamples
        self.num_bootstraps = num_bootstraps
        self.subsample_len_exp = subsample_len_exp
        self.cached_mods = {}


    
    def fingerprint(self, data):
        """
        Return a tuple of problem information sufficient to
        determine compilation-equivalence.
        """
        return (len(data),type(data))

    def run(self, data):
        if self.pure_python:
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
            f = self.fingerprint(data)
            mod = None
            if f in self.cached_mods:
                mod = self.cached_mods[f]
            else:
                mod = self.build_mod(f)
                self.cached_mods[f] = mod
            return mod.compute_blb(data)

    def compile_for( self, data, key=None ):
        f = key if key else self.fingerprint(data)
	mod = None
	if f in self.cached_mods:
	    mod = self.cached_mods[f]
	else:
	    mod = self.build_mod(f)
	    mod.backends["c++"].compile()
	    self.cached_mods[f] = mod

    def build_mod(self, key):
        template_name = ''
        if self.with_openMP:
            template_name = 'blb_omp.mako'
        elif self.with_cilk:
            template_name = 'blb_cilk.mako'
        else:
            template_name = 'blb_template.mako'
            
        fwk_args = self.set_framework_args(key)
        import asp.codegen.templating.template as template
        blb_template = template.Template(filename="templates/%s" % template_name, disable_unicode=True)
        impl_template = template.Template(filename="templates/blb_impl.mako", disable_unicode=True)
        rendered = blb_template.render( **fwk_args )
        
        
        impl_args ={}
        impl_attributes={}
        impl_args['attributes'] = impl_attributes
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

        impl_attributes['with_cilk'] = self.with_cilk

        rendered_impl = impl_template.render( **impl_args )
        
        import asp.jit.asp_module as asp_module
        mod = asp_module.ASPModule()
        mod.add_function('compute_estimate', rendered_impl)
        mod.add_function("compute_blb", rendered)

        self.set_compiler_flags(mod)
        self.set_includes(mod)
        f = open('blbout.cpp','w+')
        f.write( str(mod.backends['c++'].module.generate()) )
        f.close()
        return mod

    def __subsample(self, data, subsample_len_exp):
        subsample_len = int(len(data) ** subsample_len_exp)
        subsample = random.sample(data, subsample_len)
        return subsample

    def __bootstrap(self, data):
        bootstrap = [random.choice(data) for i in range(len(data))]
        return bootstrap

    def set_includes(self, mod):
            mod.add_header('stdlib.h')
            mod.add_header('math.h')
            mod.add_header('time.h')
            mod.add_header('numpy/ndarrayobject.h')
            
            if self.with_cilk:
                mod.add_header('cilk/cilk.h')
		mod.add_header('cilk/cilk_api.h')
            elif self.with_openMP:
                mod.add_header('omp.h')
	
            mod.add_to_init('import_array();')

    def set_compiler_flags(self, mod):
        import asp.config
        
#        mod.backends["c++"].toolchain.cflags += ['-Llibprofiler.so.0']

        if self.with_cilk: # or asp.config.CompilerDetector().detect("icc"):
            mod.backends["c++"].toolchain.cc = "icc"
            mod.backends["c++"].toolchain.cflags += ["-intel-extensions", "-fast", "-restrict"]
            mod.backends["c++"].toolchain.cflags += ["-openmp", "-fno-fnalias", "-fno-alias"]
            mod.backends["c++"].toolchain.cflags += ["-I/usr/include/x86_64-linux-gnu"]
            mod.backends["c++"].toolchain.cflags.remove('-fwrapv')
            mod.backends["c++"].toolchain.cflags.remove('-O2')
            mod.backends["c++"].toolchain.cflags.remove('-g')
            mod.backends["c++"].toolchain.cflags.remove('-fno-strict-aliasing')
        else:
            mod.backends["c++"].toolchain.cflags += ["-fopenmp", "-O3", "-msse3"]

        if mod.backends["c++"].toolchain.cflags.count('-Os') > 0:
            mod.backends["c++"].toolchain.cflags.remove('-Os')
        if mod.backends["c++"].toolchain.cflags.count('-O2') > 0:
            mod.backends["c++"].toolchain.cflags.remove('-O2')

    def set_framework_args(self, key):
        '''
        Return a dictionary containing the appropriate kwargs for redering
        the framework template.

	the key argument is a fingerprint key for this specialiser
        '''
        # estimate cache line size
        ret = {}
        if key[1] is list:
            ret['seq_type'] = 'list'
        elif key[1] is numpy.ndarray:
            ret['seq_type'] = 'ndarray'
        ret['bootstrap_unroll'] = 1
        ret['sub_n'] = int( pow( key[0], self.subsample_len_exp ) )
        ret['n_data'] = key[0]
        ret['n_subsamples'] = self.num_subsamples
        ret['n_bootstraps'] = self.num_bootstraps
        if self.with_openMP:
            # specialise this somehow.
	    ret['omp_n_threads'] =  getattr(self, 'omp_n_threads', 1)
	    print 'DEBUG: omp_n_threads = %d' % ret['omp_n_threads']
        elif self.with_cilk:
            ret['cilk_n_workers'] = getattr(self, 'cilk_n_workers', 1)
	    print 'DEBUG: cilk_nworkers = %d' % ret['cilk_n_workers']
	    ret['parallel_loop'] = 'outer'
            
        return ret
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

