
import asp.jit.asp_module as asp_module
import avroInter.avro_backend as avro_backend    
import asp.codegen.ast_tools as ast_tools
from asp.codegen.codegenScala import *
import ast

import random
import numpy
import asp.config
import inspect, ast

#from blb_convert import BLBConverter
from blb_setup import gslroot, cache_dir

def combine(blb_funcs):
    parent = (open('blb_core_parallel.scala')).read()  
    return parent + blb_funcs 

class BLB:
    known_reducers= ['mean', 'stdev', 'mean_norm', 'noop']
    def __init__(self, num_subsamples=25, num_bootstraps=100, 
                 subsample_len_exp=0.6, with_cilk=False, with_openMP=False,
                 dimension=1, pure_python=False, use_scala=False):

        self.dim = dimension
        self.with_cilk=with_cilk
        self.with_openMP = with_openMP
        self.pure_python= pure_python
        self.use_scala = use_scala

        
        self.estimate_src = inspect.getsource(self.compute_estimate)
        self.estimate_ast = ast.parse(self.estimate_src.lstrip())

        self.reduce_src = inspect.getsource(self.reduce_bootstraps)
        self.reduce_ast = ast.parse(self.reduce_src.lstrip())

        self.average_src = inspect.getsource(self.average)
        self.average_ast = ast.parse(self.average_src.lstrip())

        self.num_subsamples = num_subsamples
        self.num_bootstraps = num_bootstraps
        self.subsample_len_exp = subsample_len_exp
        self.cached_mods = {}
    
    def fingerprint(self, data):
        """
        Return a tuple of problem information sufficient to
        determine compilation-equivalence.
        """
        return tuple( [ (datum.shape,datum.dtype) for datum in data ] )
          
    def run(self, *data):
        if self.pure_python:
            self.run_python(data)
        elif self.use_scala: 
            self.run_cloud(data)
        else: #run local c++ 
            f = self.fingerprint(data)
            mod = None
            if f in self.cached_mods:
                mod = self.cached_mods[f]
            else:   
                mod = self.build_mod(f)
                self.cached_mods[f] = mod
            return mod.compute_blb(data)
        
    def run_python(self,data):
        subsample_estimates = []
        for i in range(self.num_subsamples):
            subsample = self.__subsample(data, self.subsample_len_exp)
            bootstrap_estimates = [] 
            for j in range(self.num_bootstraps):
                bootstrap = self.__bootstrap(subsample)
                estimate = self.compute_estimate(bootstrap)
                bootstrap_estimates.append(estimate)
            subsample_est = self.reduce_bootstraps(bootstrap_estimates)
            subsample_estimates.append(subsample_est)
        return self.average(subsample_estimates)
    
    def run_cloud(self,data): 
        mod = asp_module.ASPModule(cache_dir = "/root/spark/examples/target/scala-2.9.1.final/classes/", use_scala=True)
                                                                  
        scala_estimate= ast_tools.ConvertPyAST_ScalaAST().visit(self.estimate_ast) 
        scala_reduce = ast_tools.ConvertPyAST_ScalaAST().visit(self.reduce_ast)
        scala_average =  ast_tools.ConvertPyAST_ScalaAST().visit(self.average_ast)

        TYPE_DECS = (['compute_estimate', ['BootstrapData'], 'double'],             
             ['reduce_bootstraps', [('list', 'double')], 'double'],
             ['average', [('array', 'double')], 'double'])
                  
        scala_gen = SourceGenerator(TYPE_DECS)    
        rendered_scala = scala_gen.to_source(scala_estimate)+'\n' + scala_gen.to_source(scala_reduce) \
                            +'\n'+ scala_gen.to_source(scala_average)
                
        rendered_scala = combine(rendered_scala)
        rendered = avro_backend.generate_scala_object("run","",rendered_scala)             
        #NOTE: must append outer to function name above to get the classname 
        # because of how scala_object created by avro_backend           
        mod.add_function("run_outer", rendered, backend = "scala")   
        
        email_filename = data[0]
        model_filename = data[1]
        return mod.run_outer(email_filename, model_filename, self.dim, self.num_subsamples, self.num_bootstraps, self.subsample_len_exp)  
 
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
 
        
        import asp.codegen.templating.template as template
        blb_template = template.Template(filename="templates/%s" % template_name, disable_unicode=True)
        impl_template = template.Template(filename="templates/blb_impl.mako", disable_unicode=True)

            
        if key[0] % self.dim != 0:
            raise ValueError( 'Data must be of dimension %d' % self.dim )
        n_vecs = key[0] / self.dim
        vec_n = int( pow( n_vecs, self.subsample_len_exp ) )
        impl_args = {'dim': self.dim}
        impl_args['n_data'] = key[0]
        impl_args['sub_n'] = int( pow( key[0], self.subsample_len_exp ) )
        impl_args['vec_n'] = vec_n
        impl_args['n_vecs'] = n_vecs 
        impl_funcs = []
        impl_args['scalar_type'] = 'double'
    
        estimate_converter = BLBConverter( key, input_size=vec_n, weighted=True )
        estimate_cpp = estimate_converter.render( self.estimate_ast )
        impl_args['classifier'] = estimate_cpp
        impl_funcs.extend( estimate_converter.desired_funcs )
        impl_args['bootstrap_dim'] = estimate_converter.output_dim()
    
        reduce_converter = BLBConverter( impl_args['bootstrap_dim'], False )
        reduce_cpp = reduce_converter.render( self.reduce_ast )
        impl_args['bootstrap_reducer'] = reduce_cpp
        impl_funcs.extend( reduce_converter.desired_funcs )
        impl_args['subsample_dim'] = reduce_converter.output_dim()
    
        average_converter = BLBConverter( impl_args['subsample_dim'], False )
        average_cpp = average_converter.render( self.average_ast )
        impl_args['subsample_reducer'] = average_cpp
        impl_funcs.extend( average_converter.desired_funcs )
        impl_args['average_dim'] = average_converter.output_dim()
    
        impl_args['desired_funcs'] = set(impl_funcs)
        fwk_args = self.set_framework_args(key, impl_args.copy())
        rendered = blb_template.render( **fwk_args )
    
    
        """
        === This was for non sejitizing code ===
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
        """

        rendered_impl = impl_template.render( **impl_args )
        
        import asp.jit.asp_module as asp_module
        mod = asp_module.ASPModule(specializer='BLB', cache_dir=cache_dir )
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
        subsample_indicies = random.sample(range(len(data) / self.dim), subsample_len)
        subsample = []
        for index in subsample_indicies:
            subsample.extend(data[index*self.dim: (index+1)*self.dim])
        return subsample

    def __bootstrap(self, data):
        bootstrap_vectors = [data[i*self.dim:(i+1)*self.dim] for i in xrange(len(data) / self.dim)]
        bootstrap = [random.choice(bootstrap_vectors) for i in xrange(len(data) / self.dim)]
        flat = []
        for item in bootstrap:
            flat.extend(item)
        return flat
        
    def set_includes(self, mod):
        mod.add_header('stdlib.h')
        mod.add_header('math.h')
        mod.add_header('time.h')
        mod.add_header('numpy/ndarrayobject.h')
        mod.add_header('gsl_rng.h')
        mod.add_header('gsl_randist.h')
        mod.add_library( 'blas', [], libraries=['blas'] )
        mod.add_header('gsl_vector.h')
        mod.add_header('gsl_matrix.h')
        mod.add_header('gsl_blas.h')
        mod.add_header('gsl_linalg.h')
        mod.add_library( 'cblas', [], libraries=['cblas'] )
        mod.add_library( 'gsl', [gslroot, gslroot+'/randist', gslroot+'/rng', gslroot+'/matrix', gslroot+'/vector', gslroot+'/blas', gslroot+'/linalg'],
                         [gslroot+'/.libs'], ['gsl'] )
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

    def set_framework_args(self, key, vars):
        '''
        Return a dictionary containing the appropriate kwargs for redering
        the framework template.

        the key argument is a fingerprint key for this specialiser
        '''
        platform = asp.config.PlatformDetector()
        # estimate cache line size
        if key[1] is list:
            vars['seq_type'] = 'list'
        elif key[1] is numpy.ndarray:
            vars['seq_type'] = 'ndarray'
    
        vars['n_subsamples'] = self.num_subsamples
        vars['n_bootstraps'] = self.num_bootstraps
    
        if self.with_openMP:
            # specialise this somehow.
            vars['parallel_loop'] = 'outer'
            vars['omp_n_threads'] = getattr(self, 'omp_n_threads', 1 )
        elif self.with_cilk:
            vars['cilk_n_workers'] = getattr(self, 'cilk_n_workers', 1)
            print 'DEBUG: cilk_nworkers = %d' % vars['cilk_n_workers']
            vars['parallel_loop'] = 'outer'
            
        return vars
    
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