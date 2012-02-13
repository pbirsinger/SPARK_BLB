<%doc>
Framework template for OpenMP implementation of BLB algorithm.
see 'blb_template.mako' for documentation on templating variables.
</%doc>
typedef ${scalar_type} scalar_t;
#define NPY_SCALAR ${ 'NPY_FLOAT32' if scalar_type is 'float' else 'NPY_FLOAT64' }

void subsample_and_load( scalar_t* data, scalar_t* out, gsl_rng* rng ){
        for( int i = 0; i<${vec_n}; i++ ){
	    unsigned int index = gsl_rng_get(rng) % ${n_vecs};
	    memcpy(out + i*${dim}, data + index*${dim}, ${dim}*sizeof(scalar_t) );
        }
}

void bootstrap( unsigned int* out, gsl_rng * rng ){
  double norm = ${ vec_n*(1.0/vec_n) };
  double sum_p = 0.0;
  double p = ${ 1.0/vec_n };
  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (int k = 0; k < ${vec_n}; k++)
    {
     out[k] = gsl_ran_binomial (rng, p / (norm - sum_p), ${n_vecs} - sum_n);
     sum_p += p;
     sum_n += out[k];
    }
}
<%doc>
determine header
allocate temporaries
effect computation
free temporaries
write to output
</%doc>

## list is the default type.
%if seq_type is UNDEFINED or seq_type == 'list':

PyObject * compute_blb(PyObject * data) {
    Py_INCREF( data );
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_SCALAR, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( py_arr );
%elif seq_type == 'ndarray':

PyObject* compute_blb( PyObject* data ){
//  printf("Made it into C\n");
  Py_INCREF( data );
  scalar_t * c_arr = (scalar_t*) PyArray_DATA( data );
%endif    

    scalar_t * const subsample_estimates = (scalar_t*) calloc(${n_subsamples*subsample_dim}, sizeof(scalar_t));
    scalar_t * const bootstrap_estimates = (scalar_t*) calloc(${n_bootstraps*omp_n_threads*bootstrap_dim}, sizeof(scalar_t));
    scalar_t * const subsample_values = (scalar_t*) calloc(${vec_n*omp_n_threads*dim}, sizeof(scalar_t));
    unsigned int* const bootstrap_weights = (unsigned int*) calloc(${vec_n*omp_n_threads}, sizeof(scalar_t)); 

%if parallel_loop is UNDEFINED or parallel_loop == 'outer':
    srand( time(NULL) );
    gsl_rng** rngs = (gsl_rng**) malloc( ${omp_n_threads}*sizeof(gsl_rng*) );
    for(int i = 0; i<${omp_n_threads}; i++){
	rngs[i] = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngs[i], rand());
    }
    #pragma omp parallel for schedule(dynamic) num_threads(${omp_n_threads})
    for (int i = 0; i < ${n_subsamples}; i++) {
        int tid = omp_get_thread_num();
        unsigned int* local_weights = bootstrap_weights+(${vec_n}*tid);
        scalar_t* local_values = subsample_values+(${vec_n*dim}*tid);
        scalar_t* local_estimates = bootstrap_estimates+(${n_bootstraps*bootstrap_dim}*tid);
	gsl_rng* rng = rngs[tid];
        subsample_and_load(c_arr, local_values, rng);
	for (int j=0; j < ${n_bootstraps}; j++) {
            bootstrap(local_weights, rng );
            compute_estimate(local_values, local_weights, ${vec_n}, local_estimates+j*${bootstrap_dim} );
        }
        reduce_bootstraps(local_estimates, ${n_bootstraps}, subsample_estimates+(i*${subsample_dim}) );
    }
    for( int i=0; i<${omp_n_threads};i++ ){
	gsl_rng_free(rngs[i]);
    }
    free(rngs);
%elif parallel_loop == 'inner':

    gsl_rng** rngs = (gsl_rng**) calloc( ${omp_n_threads}, sizeof(gsl_rng*));

    #pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
    for( int i = 0; i< ${omp_n_threads}; i++ ){
    rngs[i] = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set( rngs[i], time(NULL)*i );
    }
    for( int i = 0; i< ${n_subsamples}; i++ ){
	printf("Starting subsample %d\n", i);
	subsample_and_load( c_arr, subsample_values, rngs[i % ${omp_n_threads}] );
	#pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
	for( int j = 0; j< ${omp_n_threads}; j++ ){
	    int tid = omp_get_thread_num();
	    loaded_bootstrap( bootstrap_weights+(${sub_n}*tid), rngs[tid] );
	}
 	#pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
	for( int j = 0; j< ${n_bootstraps}; j++ ){
	    int tid = omp_get_thread_num();
	    unsigned int* local_weights = bootstrap_weights+(${sub_n}*tid);
	    unsigned int seed = tid* time(NULL);
	    permutation_bootstrap(local_weights, &seed);
	    bootstrap_estimates[j] = compute_estimate(subsample_values, local_weights, ${sub_n});
	}
	subsample_estimates[i] = reduce_bootstraps( bootstrap_estimates, ${n_bootstraps} );
    }
    for( int i=0; i<${omp_n_threads};i++ ){
	gsl_rng_free(rngs[i]);
    }
    free(rngs);
%endif
    %if average_dim == 1:
    scalar_t theta = 0.0;
    average(subsample_estimates, ${n_subsamples}, &theta ); 
    %else:
    scalar_t* theta = (scalar_t*) calloc( ${average_dim}, sizeof(scalar_t) );
    average( subsample_estimates, ${n_subsamples}, theta );
    %endif

    free(subsample_values);
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_weights);

%if seq_type is UNDEFINED or seq_type == 'list':
    Py_DECREF( py_arr );
%endif
    Py_DECREF( data );
    %if average_dim == 1:
    return PyFloat_FromDouble( theta );
    %else:
    npy_intp dim[1] = { ${average_dim} };
    return PyArray_SimpleNewFromData( 1, dim, NPY_SCALAR, theta );
    %endif
}



