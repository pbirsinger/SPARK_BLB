<%doc>
USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
  Templating variables in use:

  vec_n: The size in terms of number of vectors to be subsampled (sub_n / dim)
  n_vec: The number of vectors in the data
  n_bootstraps: The number of bootstraps to compute per subsample
  seq_type: The python type of the data sequence, should be list or ndarray
</%doc>

typedef ${scalar_type} scalar_t;
#define NPY_SCALAR ${ 'NPY_FLOAT32' if scalar_type is 'float' else 'NPY_FLOAT64' }

<%def name="bigrand( fname )">
inline unsigned long ${fname}( unsigned int* seed ){
        unsigned long ret = rand_r(seed);
        return ((ret << 32) | rand_r(seed));
}
</%def>

<%def name="littlerand( fname )">
inline unsigned int ${fname}( unsigned int* seed ){
    return rand_r(seed);
}
</%def>

void bootstrap( unsigned int* out, gsl_rng* rng ){
  double norm = ${ n_vecs*(1.0/n_vecs) };
  double sum_p = 0.0;
  double p = ${ 1.0/n_vecs };
  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (int k = 0; k < ${n_vecs}; k++)
    {
     out[k] = gsl_ran_binomial (rng, p / (norm - sum_p), ${n_vecs} - sum_n);
     sum_p += p;
     sum_n += out[k];
    }


}


## list is the default type.
%if seq_type == 'list':

PyObject* compute_bootstrap( PyObject*  data ){
    Py_INCREF(data);
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_SCALAR, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( py_arr );

%elif seq_type == 'ndarray':

PyObject* compute_bootstrap( PyObject* data ){
    printf("Inside Bootstrap\n"); 
    Py_INCREF( data );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( data );

%endif

    //note that these are never cleared; We always fill them up
    //with the appropriate data before perform calculations on them.
    unsigned int * bootstrap_weights = (unsigned int*) calloc( ${omp_n_threads*n_vecs}, sizeof(unsigned int) );
    scalar_t * bootstrap_estimates = (scalar_t*) calloc( ${n_bootstraps*bootstrap_dim}, sizeof(scalar_t) );
    gsl_rng** rngs = (gsl_rng**) malloc( ${omp_n_threads}*sizeof(gsl_rng*) );
    srand( time( NULL ) );
    for( int i = 0; i< ${omp_n_threads}; i++ ){
    rngs[i] = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set( rngs[i], rand() );
    }
    #pragma omp parallel for schedule(dynamic) num_threads(${omp_n_threads})
    for( int j=0; j<${n_bootstraps}; j++ ){
//      printf("Computing bootstrap number %d\n", j);
	int tid = omp_get_thread_num();
	unsigned int* local_weights = bootstrap_weights+(${n_vecs}*tid);
	gsl_rng* rng = rngs[tid];
	scalar_t* local_estimates = bootstrap_estimates+(j*${bootstrap_dim});
	bootstrap( local_weights, rng );
        //printf("bootstrap weights: ");
        //printArray(bootstrap_weights, 0, ${n_vecs});
	printf("Iteration %d, about to compute\n", j);
        compute_estimate(c_arr, local_weights, ${n_vecs}, local_estimates );
        //printf("Bootstrap estimate for bootstrap dim %d: ", ${bootstrap_dim});
        //printArray(bootstrap_estimates, 0, ${n_bootstraps*bootstrap_dim});
    }
    %if average_dim == 1:
    scalar_t theta = 0.0;
    average( bootstrap_estimates, ${n_bootstraps}, &theta );
    %else:
    scalar_t* theta = (scalar_t*) calloc( ${average_dim}, sizeof(scalar_t) );
    average( bootstrap_estimates, ${n_bootstraps}, theta );
    %endif

    for( int i = 0; i<${omp_n_threads}; i++ ){
	gsl_rng_free(rngs[i]);
    }
    free( rngs ); 
    free( bootstrap_weights );
    free( bootstrap_estimates );


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

void printArray(scalar_t* arr, int start, int end) {
     for (int i=start; i<end; i++) {
     	 printf("%f, ", arr[i]);
     }
     printf("\n");
}

void printArray(unsigned int* arr, int start, int end) {
     for (int i=start; i<end; i++) {
         printf("%d, ", arr[i]);
     }
     printf("\n");
}