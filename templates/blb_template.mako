<%doc>
USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
  Templating variables in use:
  sub_n: The size b(n) of data to be subsampled
  n_data: The initial data size
  n_subsamples: The number of subsamples to take
  n_bootstraps: The number of bootstraps to compute per subsample
  subsmaple_threshold: the probability parameter for the subsample rng
  seq_type: The python type of the data sequence, should be list or ndarray
</%doc>

void bootstrap( const unsigned int* in, unsigned int* out, unsigned int* seed ){
    <%
        if bootstrap_unroll is UNDEFINED:
            b = 1
        else:
            b = bootstrap_unroll
    %>
  for( int i=0; i< ${sub_n/b}; i++ ){
    % for i in range(b):
    out[i*${b} + ${i}] = in[ rand_r(seed) % ${sub_n} ];
    % endfor
  }
  % for i in range(sub_n % b):
  out[${sub_n-1-i}] = in[ rand_r(seed) % ${sub_n} ];
  % endfor
}

void loaded_bootstrap( unsigned int* out, unsigned int * seed ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = rand_r(seed) % ${sub_n};
        }
}


void subsample_and_load( float* data, float* out, unsigned int* seed ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = data[ rand_r(seed) % ${n_data} ];
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
  float * subsample_estimates = (float*) calloc( ${n_subsamples}, sizeof(float) );
  float * subsample_values = (float*) calloc( ${sub_n}, sizeof(float) );
  unsigned int * bootstrap_indicies = (unsigned int*) calloc( ${sub_n}, sizeof(unsigned int) );
  float * bootstrap_estimates =  (float*) calloc( ${sub_n}, sizeof(float) );
  unsigned int cell = time(NULL);
  for( int i=0; i<${n_subsamples}; i++ ){

    subsample_and_load( c_arr, subsample_values, &cell );
    for( int j=0; j<${n_bootstraps}; j++ ){

      loaded_bootstrap( bootstrap_indicies, &cell );
      bootstrap_estimates[j] = compute_estimate( subsample_values, bootstrap_indicies, ${sub_n} );

    }
    subsample_estimates[i] = reduce_bootstraps( bootstrap_estimates, ${n_bootstraps} );
  }

  float theta = average( subsample_estimates, ${n_subsamples} );

  free( subsample_estimates );
  free( bootstrap_indicies );
  free( bootstrap_estimates );
  free( subsample_values );
%if seq_type is UNDEFINED or seq_type == 'list':
  Py_DECREF( py_arr );
%endif
  Py_DECREF( data );
  return PyFloat_FromDouble(theta);
}




