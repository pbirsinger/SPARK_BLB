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

/*
void subsample_and_load( float* data, float* out, gsl_rng* rng ){
        for( int i = 0; i<${vec_n}; i++ ){
	    unsigned int index = gsl_rng_get(rng) % ${n_vecs};
            vvc( data+(${dim}*index), out+(i*${dim}), ${dim} );
        }
}
*/


## list is the default type.
%if seq_type == 'list':

PyObject* compute_blb( PyObject*  data ){
    Py_INCREF(data);
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_FLOAT32, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    float * c_arr = (float*) PyArray_DATA( py_arr );

%elif seq_type == 'ndarray':

PyObject* compute_blb( PyObject* data ){
    Py_INCREF( data );
    float * c_arr = (float*) PyArray_DATA( data );

%endif

    //note that these are never cleared; We always fill them up
    //with the appropriate data before perform calculations on them.
    unsigned int * bootstrap_weights = (unsigned int*) calloc( ${n_vecs}, sizeof(unsigned int) );
    float * bootstrap_estimates = (float*) calloc( ${n_bootstraps*bootstrap_dim}, sizeof(float) );
    float * result = (float*) calloc(${bootstrap_dim}, sizeof(float));
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set( rng, time(NULL) );

    for( int j=0; j<${n_bootstraps}; j++ ){
      printf("Computing bootstrap number %d\n", j);
      bootstrap( bootstrap_weights, rng );
      //printf("bootstrap weights: ");
      //printArray(bootstrap_weights, 0, ${n_vecs});
      compute_estimate(c_arr, bootstrap_weights, ${n_vecs}, bootstrap_estimates + j*${bootstrap_dim} );
      printf("Bootstrap estimate for bootstrap dim %d: ", ${bootstrap_dim});
      printArray(bootstrap_estimates, 0, ${n_bootstraps*bootstrap_dim});
    }
    reduce_bootstraps( bootstrap_estimates, ${n_bootstraps}, result);
    printf("Result is %f", result[0]);
//    printArray(result, 0, ${bootstrap_dim});
//    printf("Subsample estimates: ");
//    printArray(subsample_estimates + i*${subsample_dim}, 0, ${subsample_dim});
 
  float theta = result[0];
//  average( subsample_estimates, ${n_subsamples}, &theta );
 
  gsl_rng_free(rng);
  free( bootstrap_weights );
  free( bootstrap_estimates );
  free( result );

%if seq_type is UNDEFINED or seq_type == 'list':
  Py_DECREF( py_arr );
%endif
  Py_DECREF( data );
  return PyFloat_FromDouble(theta);
}

void printArray(float* arr, int start, int end) {
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