<%doc>
USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
  Templating variables in use:
  sub_n: The size b(n) of data to be subsampled
  vec_n: The size in terms of number of vectors to be subsampled (sub_n / dim)
  n_data: The initial data size
  n_vec: The number of vectors in the data
  n_subsamples: The number of subsamples to take
  n_bootstraps: The number of bootstraps to compute per subsample
  subsmaple_threshold: the probability parameter for the subsample rng
  seq_type: The python type of the data sequence, should be list or ndarray
</%doc>
typedef ${scalar_type} scalar_t;
##define scalar_t ${scalar_type}
#define NPY_SCALAR ${ 'NPY_FLOAT32' if scalar_type is 'float' else 'NPY_FLOAT64' }
void printArray( scalar_t*, int, int );
void printArray( unsigned int*, int, int );

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

void subsample_and_load( scalar_t* data, scalar_t* out, gsl_rng* rng ){
        for( int i = 0; i<${vec_n}; i++ ){
	    unsigned int index = gsl_rng_get(rng) % ${n_vecs};
            vvc( data+(${dim}*index), out+(i*${dim}), ${dim} );
        }
}



## list is the default type.
%if seq_type == 'list':

PyObject* compute_blb( PyObject*  data ){
    Py_INCREF(data);
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_SCALAR, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( py_arr );

%elif seq_type == 'ndarray':

PyObject* compute_blb( PyObject* data ){
    Py_INCREF( data );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( data );

%endif
    //note that these are never cleared; We always fill them up
    //with the appropriate data before perform calculations on them.
    scalar_t * subsample_estimates = (scalar_t*) calloc( ${n_subsamples*subsample_dim}, sizeof(scalar_t) );
    scalar_t * subsample_values = (scalar_t*) calloc( ${vec_n*dim}, sizeof(scalar_t) );
    unsigned int * bootstrap_weights = (unsigned int*) calloc( ${vec_n}, sizeof(unsigned int) );
    scalar_t * bootstrap_estimates = (scalar_t*) calloc( ${n_bootstraps*bootstrap_dim}, sizeof(scalar_t) );
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set( rng, time(NULL) );

    for( int i=0; i<${n_subsamples}; i++ ){
//    	printf("Loading subsample number %d\n", i);
        subsample_and_load( c_arr, subsample_values, rng );
        for( int j=0; j<${n_bootstraps}; j++ ){
//	    printf("Computing bootstrap number %d\n", j);
            bootstrap( bootstrap_weights, rng );
            compute_estimate( subsample_values, bootstrap_weights, ${vec_n}, bootstrap_estimates + j*${bootstrap_dim} );
	    if( j == 0 ){
//		printf("Bootstrap_estimate: ");
//		printArray( bootstrap_estimates + j*${bootstrap_dim}, 0, ${bootstrap_dim} );
	    }
        }
    reduce_bootstraps( bootstrap_estimates, ${n_bootstraps}, subsample_estimates + i*${subsample_dim} );
//    printf("Subsample estimates: ");
//    printArray(subsample_estimates + i*${subsample_dim}, 0, ${subsample_dim});
    }
  
  %if average_dim == 1:
  scalar_t theta = 0.0;
  average( subsample_estimates, ${n_subsamples}, &theta );
  %else:
  scalar_t* theta = (scalar_t*) calloc( ${average_dim}, sizeof(scalar_t) );
  average( subsample_estimates, ${n_subsamples}, theta );
  %endif

  gsl_rng_free(rng);
  free( subsample_estimates );
  free( bootstrap_weights );
  free( bootstrap_estimates );
  free( subsample_values );

%if seq_type is UNDEFINED or seq_type == 'list':
  Py_DECREF( py_arr );
%endif
  Py_DECREF( data );
  %if average_dim == 1:
  return PyFloat_FromDouble( theta );
  %else:
  npy_intp dim[1] = { ${average_dim} };
  return PyArray_SimpleNewFromData( 1, dim  , NPY_SCALAR, theta ); 
  %endif
}

void printArray(scalar_t* arr, int start, int end) {
     for (int i=start; i<end; i++) {
     	 printf("%.4f, ", arr[i]);
     }
     printf("\n");
}

void printArray(unsigned int* arr, int start, int end) {
     for (int i=start; i<end; i++) {
         printf("%d, ", arr[i]);
     }
     printf("\n");
}
