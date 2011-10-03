<%doc>
USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
 Templating variables in use:
 sub_n: The size b(n) of data to be subsampled
 n_data: The initial data size
 n_subsamples: The number of subsamples to take
 n_bootstraps: The number of bootstraps to compute per subsample
</%doc>

void bootstrap( unsigned int* in, unsigned int* out ){
//  printf("About to bootstrap");
  for( int i=0; i< ${sub_n}; i++ ){
    int index = rand() % ${sub_n};
    out[i] = in[index];
  }
}

char subsampled[ ${n_data} ];
void subsample( unsigned int* out ){
//  printf("About to subsample");
  int size_out = ${sub_n};
  while( size_out > 0 ){
    unsigned int index = rand() % ${n_data};
    if( subsampled[index] ){
      continue;
    } else {
      subsampled[index] = 1;
      out[ ${sub_n} - size_out ] = index;
      size_out--;
    }
  }
  for( int i=0; i<${sub_n}; i++ ){
       subsampled[ out[i] ] = 0;
  }
}

PyObject* compute_blb( PyObject*  data ){
  Py_INCREF(data);
  PyObject * py_arr = PyArray_FROM_OTF( data, NPY_FLOAT32, NPY_IN_ARRAY );
  Py_INCREF( py_arr );
  float * c_arr = (float*) PyArray_DATA( py_arr );
  srand(time(NULL));
  memset( subsampled, 0, ${n_data} );

  //note that these are never cleared; We always fill them up
  //with the appropriate data before perform calculations on them.
  float * subsample_estimates = (float*) calloc( ${n_subsamples}, sizeof(float) );
  float * bootstrap_estimates =  (float*) calloc( ${sub_n}, sizeof(float) );
  unsigned int * subsample_indicies = (unsigned int*) calloc( ${sub_n}, sizeof(unsigned int) );
  unsigned int * bootstrap_indicies = (unsigned int*) calloc( ${sub_n}, sizeof(unsigned int) );
  for( int i=0; i<${n_subsamples}; i++ ){
    subsample( subsample_indicies );
    for( int j=0; j<${n_bootstraps}; j++ ){
      bootstrap( subsample_indicies, bootstrap_indicies );
      bootstrap_estimates[j] = compute_estimate( c_arr, bootstrap_indicies, ${sub_n} );
    }
    subsample_estimates[i] = reduce_bootstraps( bootstrap_estimates, ${n_bootstraps} );
  }
  float theta = average( subsample_estimates, ${n_subsamples} );
  free( subsample_indicies );
  free( subsample_estimates );
  free( bootstrap_estimates );
  free( bootstrap_indicies );
  Py_DECREF( py_arr );
  Py_DECREF( data );
  return PyFloat_FromDouble(theta);
}
