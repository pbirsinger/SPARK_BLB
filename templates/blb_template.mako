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

<<<<<<< HEAD
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

%if n_data >= 1<<32:
    ${bigrand('sub_rand')}
%else:
    ${littlerand('sub_rand')}
%endif
void subsample_and_load( float* data, float* out, unsigned int* seed ){
        for( int i = 0; i<${vec_n}; i++ ){
            memcpy(out + i*${dim}, data + (sub_rand(seed) % ${n_vec}) * ${dim}, ${dim} * sizeof(float));
        }
}
%if sub_n >= 1<<32:
    ${bigrand('boot_rand')}
%else:
    ${littlerand('boot_rand')}
%endif
void loaded_bootstrap( unsigned int* out, unsigned int * seed ){
        for( int i = 0; i<${vec_n}; i++ ){
            out[i] = boot_rand(seed) % ${vec_n};
=======
#define DATA_SIZE ${n_data}
#define SUBSAMPLE_SIZE ${sub_n}

// Canabalized from gsl_ran_multinomial
void bootstrap( unsigned int* out, gsl_rng* rng ){
  double norm = ${ sub_n*(1.0/sub_n) };
  double sum_p = 0.0;
  double p = ${ 1.0/sub_n };
  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (int k = 0; k < ${sub_n}; k++)
    {
     out[k] = gsl_ran_binomial (rng, p / (norm - sum_p), ${n_data} - sum_n);
     sum_p += p;
     sum_n += out[k];
    }


}

void subsample_and_load( float* data, float* out, gsl_rng* rng ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = data[ gsl_rng_get(rng) % ${n_data} ];
>>>>>>> 895e8963812d9de0345657e9377a5a5f2445a8f8
        }
}


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
<<<<<<< HEAD
  //float * subsample_estimates = (float*) calloc( ${n_subsamples}, sizeof(float) );
  float subsample_estimates[${n_subsamples * subsample_dim}];
  //float * subsample_values = (float*) calloc( ${sub_n}, sizeof(float) );
  float subsample_values[${sub_n * dim}];
  //unsigned int * bootstrap_indicies = (unsigned int*) calloc( ${sub_n}, sizeof(unsigned int) );
  unsigned int bootstrap_indicies[${vec_n}];
  //float * bootstrap_estimates =  (float*) calloc( ${n_bootstraps}, sizeof(float) );
  float bootstrap_estimates[${n_bootstraps * bootstrap_dim}];
  unsigned int cell = time(NULL);
  for( int i=0; i<${n_subsamples}; i++ ){
    subsample_and_load( c_arr, subsample_values, &cell );
    for( int j=0; j<${n_bootstraps}; j++ ){
      loaded_bootstrap( bootstrap_indicies, &cell );
      compute_estimate( subsample_values, bootstrap_indicies, ${vec_n}, bootstrap_estimates + j*${bootstrap_dim} );
      unsigned int k = j*${bootstrap_dim};
      bool hasNan = false;
      %for l in xrange(bootstrap_dim):
      	   if (bootstrap_estimates[k + ${l}] != bootstrap_estimates[k + ${l}]) {
	      hasNan = true;
	   }
      %endfor
      if (hasNan) {
            printf("Estimate computed for bootstrap %d in subsample %d: ", j, i);
      }
      %for l in xrange(bootstrap_dim):
      	   if (hasNan) {
	       printf("%f ", bootstrap_estimates[k + ${l}]);
	   }
      %endfor
      if (hasNan) {
          printf("\n");
      }
    }
    reduce_bootstraps( bootstrap_estimates, ${n_bootstraps}, subsample_estimates + i*${subsample_dim} );
    unsigned int m = i * ${subsample_dim};
    bool hasNan = false;
    %for l in xrange(subsample_dim):
    	 if (subsample_estimates[m+${l}] != subsample_estimates[m+${l}]) {
	    hasNan = true;
	 }
    %endfor
    if (hasNan) {
         printf("Bootstrap reduced for subsample %d: \n", i);
=======
  float * subsample_estimates = (float*) calloc( ${n_subsamples}, sizeof(float) );
  float * subsample_values = (float*) calloc( ${sub_n}, sizeof(float) );
  unsigned int * bootstrap_weights = (unsigned int*) calloc( ${sub_n}, sizeof(unsigned int) );
  float * bootstrap_estimates =  (float*) calloc( ${sub_n}, sizeof(float) );
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, time(NULL) ); 
  for( int i=0; i<${n_subsamples}; i++ ){

    subsample_and_load( c_arr, subsample_values, rng );
    for( int j=0; j<${n_bootstraps}; j++ ){

      bootstrap( bootstrap_weights, rng );
      bootstrap_estimates[j] = compute_estimate( subsample_values, bootstrap_weights, ${sub_n} );

>>>>>>> 895e8963812d9de0345657e9377a5a5f2445a8f8
    }
    %for l in xrange(subsample_dim):
    	 if (hasNan) {
    	    printf("%f ", subsample_estimates[m + ${l}]);
	 }
    %endfor
    printf("\n");
  }

  float theta = 0;
  average( subsample_estimates, ${n_subsamples}, &theta );

//  free( subsample_estimates );
//  free( bootstrap_indicies );
//  free( bootstrap_estimates );
//  free( subsample_values );

<<<<<<< HEAD
=======
  gsl_rng_free(rng);
  free( subsample_estimates );
  free( bootstrap_weights );
  free( bootstrap_estimates );
  free( subsample_values );
>>>>>>> 895e8963812d9de0345657e9377a5a5f2445a8f8
%if seq_type is UNDEFINED or seq_type == 'list':
  Py_DECREF( py_arr );
%endif
  Py_DECREF( data );
  return PyFloat_FromDouble(theta);
}




