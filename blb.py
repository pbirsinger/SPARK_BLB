"""
Main class representing the BLB algorithm.

"""

import random
import numpy
import time
#import asp.config
import inspect, ast
from convert.blb_convert import BLBConverter, create_data_model
import convert.blb_convert_data_model as data_model
from blb_setup import gslroot, cache_dir

import asp.jit.asp_module as asp_module
import asp.avro_inter.avro_backend as avro_backend    
import asp.codegen.ast_tools as ast_tools
from asp.codegen.codegen_scala import *
import ast

def combine(blb_funcs):
    parent = (open('distr_support/blb_core_parallel.scala')).read()  
    return parent + blb_funcs 
	
class BLB:
    known_reducers= ['mean', 'stdev', 'mean_norm', 'noop']
    def __init__(self, num_subsamples=25, num_bootstraps=100, 
                 subsample_len_exp=0.5, with_cilk=False, with_openMP=False,
	   	 dimension=1, pure_python=False, with_scala=False):

        self.dim = dimension
        self.with_cilk=with_cilk
	self.with_openMP = with_openMP
	self.pure_python= pure_python
	self.with_scala = with_scala

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
    	elif self.with_scala:
    	    self.run_distributed(data)
        else:
            f = self.fingerprint(data)
            mod = None
            if f in self.cached_mods:
                mod = self.cached_mods[f]
            else:
                mod = self.build_mod(f)
                self.cached_mods[f] = mod
            return mod.compute_blb(*data)
    
    def run_python(self, data):
        subsample_estimates = []
        for i in range(self.num_subsamples):
            subsample = self.__subsample(data, self.subsample_len_exp)
            bootstrap_estimates = [] 
            for j in range(self.num_bootstraps):
                bootstrap = self.__bootstrap(subsample)
                estimate = self.compute_estimate(bootstrap)
                bootstrap_estimates.append(estimate)
#                    print "***PYTHON bootstrap estimate for bootstrap " + str(j) + " and subsample " + str(i) + " is " + str(estimate)
            subsample_est = self.reduce_bootstraps(bootstrap_estimates)
            subsample_estimates.append(subsample_est)
#                print "***PYTHON subsample estimate for subsample " + str(i) + " is " + str(subsample_est)
        return self.average(subsample_estimates)
    
    def run_distributed(self,data): 
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
        
	#NOTE: must add dependencies in make_dependency_jar so that slave nodes will receive proper files
	time_stamp = str(int(round(time.time() * 1000)))
	os.system('/root/BLB/distr_support/make_dependency_jar ' + '/root/BLB/distr_support/dependencies/' + time_stamp)
	os.environ['DEPEND_LOC'] = '/root/BLB/distr_support/dependencies/' + time_stamp +'/depend.jar'

        email_filename = data[0]
        model_filename = data[1]
        return mod.run_outer(email_filename, model_filename, self.dim, self.num_subsamples, self.num_bootstraps, self.subsample_len_exp)  
 
    def compile_for( self, key  ):
	mod = None
	if key in self.cached_mods:
	    mod = self.cached_mods[key]
	else:
	    mod = self.build_mod(key)
	    mod.backends["c++"].compile()
	    self.cached_mods[key] = mod

    def create_input_model( self, key, vec_n ):
	#Subsample and load cutoff size
	cutoff = 10*1024*1024
	models = []
	for k in key:
	    print k
	    size = reduce( int.__mul__, k[0][1:] ) if len(k[0]) > 1 else 1
	    size *= vec_n
	    model = None
	    if size <= cutoff:
		model = data_model.DataModel( k[1], list(k[0]) ) 
	    else:
		model = data_model.IndirectDataModel( k[1], list(k[0]) )
	    model.set_len( k[0][0] )
	    model.dimensions[0] = vec_n
	    models.append( model )
	return models
	    
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

        #leading dimension of first data object tuple 
	n_vecs = key[0][0][0]

	vec_n = int( pow( n_vecs, self.subsample_len_exp ) )
	impl_args = {}
        #impl_args = {'dim': self.dim}
        #impl_args['n_data'] = key[0]
        #impl_args['sub_n'] = int( pow( key[0], self.subsample_len_exp ) )
        impl_args['vec_n'] = vec_n
	impl_args['n_vecs'] = n_vecs 
	impl_funcs = []
	impl_args['scalar_type'] = 'double'


	estimate_converter = BLBConverter( self.create_input_model(key, vec_n), input_size=vec_n, weighted=True )
	estimate_cpp = estimate_converter.render( self.estimate_ast )
	impl_args['classifier'] = estimate_cpp
	impl_funcs.extend( estimate_converter.desired_funcs )
 	impl_args['bootstrap_model'] = estimate_converter.get_ret_model()
	arg_model = [ model.clone() for model in estimate_converter.arg_model ]
	for i in range(len(arg_model)):
	    arg_model[i].dimensions = key[i][0]
	impl_args['arg_model'] = arg_model

	reduce_input = estimate_converter.get_ret_model().clone()
	reduce_input.dimensions.insert( 0, self.num_bootstraps )
	reduce_input.set_len( self.num_bootstraps )
	reduce_converter = BLBConverter( [reduce_input], input_size=self.num_bootstraps )
	reduce_cpp = reduce_converter.render( self.reduce_ast )
	impl_args['bootstrap_reducer'] = reduce_cpp
	impl_funcs.extend( reduce_converter.desired_funcs )
	impl_args['subsample_model'] = reduce_converter.get_ret_model()

	average_input = reduce_converter.get_ret_model().clone()
	average_input.dimensions.insert( 0, self.num_subsamples )
	average_input.set_len( self.num_subsamples )
	average_converter = BLBConverter( [average_input], input_size=self.num_subsamples )
	average_cpp = average_converter.render( self.average_ast )
	impl_args['subsample_reducer'] = average_cpp
	impl_funcs.extend( average_converter.desired_funcs )
	impl_args['average_model'] = average_converter.get_ret_model()

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
#	    mod.add_library( 'cblas', [], libraries=['cblas'] )
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
	n_vecs = key[0][0][0]
	models = create_data_model( key, int( n_vecs**self.subsample_len_exp ) )
	vars['seq_type'] = 'ndarray'
#	vars['arg_model'] = models
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
