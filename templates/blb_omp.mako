<%doc>
Framework template for OpenMP implementation of BLB algorithm.
see 'blb_template.mako' for documentation on templating variables.
</%doc>


//USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
void bootstrap(unsigned int* in, unsigned int* out, unsigned int* seed) {
        for (int i = 0; i < ${sub_n}; i++) {
                int index = rand_r(seed) % ${sub_n};
                out[i] = in[index];
        }
}


void subsample_and_load( float* data, float* out, unsigned int* seed ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = data[ rand_r(seed) % ${n_data} ];
        }
}

void loaded_bootstrap( unsigned int* out, unsigned int * seed ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = rand_r(seed) % ${sub_n};
        }
}


<%doc>
We're going to try a little experiment.
We're goiing to play it fast and loose, relying upon the period
of our RNG to not select the same index twice in short sequence.
</%doc>
void subsample(unsigned int* out, unsigned int* seed) {
        unsigned int index = 0;
        for (int i=0; i < ${sub_n}; i++) {
                index =  rand_r(seed) % ${n_data};
                out[i] = index;
        }
}

## list is the default type.
%if seq_type == 'list':

PyObject * compute_blb(PyObject * data) {
    Py_INCREF( data );
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_FLOAT32, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    float * c_arr = (float*) PyArray_DATA( py_arr );
%elif seq_type == 'ndarray':

PyObject* compute_blb( PyObject* data ){
  Py_INCREF( data );
  float * c_arr = (float*) PyArray_DATA( data );
  printf("Successfully extracted data\n");
%endif


    float * subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float * bootstrap_estimates = (float*) calloc(${n_bootstraps*omp_n_threads}, sizeof(float));
    float * subsample_values = (float*) calloc(${sub_n*omp_n_threads}, sizeof(float));
    unsigned int * bootstrap_indicies = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically


    #pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
    for (int i = 0; i < ${n_subsamples}; i++) {
        int tid = omp_get_thread_num();
        unsigned int* local_indicies = bootstrap_indicies+(${sub_n}*tid);
        float* local_values = subsample_values+(${sub_n}*tid);
        float* local_estimates = bootstrap_estimates+(${n_bootstraps}*tid);
        unsigned int seed = time(NULL);
        subsample_and_load(c_arr, local_values, &seed);
        for (int j=0; j < ${n_bootstraps}; j++) {
            loaded_bootstrap(local_indicies, &seed);
            local_estimates[j] = compute_estimate(local_values, local_indicies, ${sub_n});
        }

        subsample_estimates[i] = reduce_bootstraps(local_estimates, ${n_bootstraps});
    }

    float theta = average(subsample_estimates, ${n_subsamples});
    free(subsample_values);
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_indicies);

    //printf("About to DECREF\n");
%if seqtype == 'list':
    Py_DECREF( py_arr );
%endif
    Py_DECREF( data );
    return PyFloat_FromDouble(theta);
}



