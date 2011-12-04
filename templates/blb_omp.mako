<%doc>
Framework template for OpenMP implementation of BLB algorithm.
see 'blb_template.mako' for documentation on templating variables.
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

%if n_data >= 1<<32:
    ${bigrand('sub_rand')}
%else:
    ${littlerand('sub_rand')}
%endif
void subsample_and_load( float* data, float* out, unsigned int* seed ){
        for( int i = 0; i<${sub_n}; i++ ){
	    memcpy(out + i*${dim}, data + (sub_rand(seed) % ${n_data}) * ${dim}, ${dim});
        }
}
%if sub_n >= 1<<32:
    ${bigrand('boot_rand')}
%else:
    ${littlerand('boot_rand')}
%endif
void loaded_bootstrap( unsigned int* out, unsigned int * seed ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = boot_rand(seed) % ${sub_n};
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
%if seq_type is UNDEFINED or seq_type == 'list':

PyObject * compute_blb(PyObject * data) {
    printf("Made it into C\n");
    Py_INCREF( data );
    printf("About to extract array\n");
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_FLOAT32, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    printf("Getting c_arr\n");
    float * c_arr = (float*) PyArray_DATA( py_arr );
%elif seq_type == 'ndarray':

PyObject* compute_blb( PyObject* data ){
  Py_INCREF( data );
  float * c_arr = (float*) PyArray_DATA( data );
%endif    



    //float * const subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float subsample_estimates[${n_subsamples}];
    //printf("subsample_estimates: %x\n", (unsigned int) subsample_estimates);
    //float * const bootstrap_estimates = (float*) calloc(${n_bootstraps*omp_n_threads}, sizeof(float));
    float bootstrap_estimates[${n_bootstraps * omp_n_threads}];
    //printf("bootstrap_estimates: %x\n", (unsigned int) bootstrap_estimates);
    //float * const subsample_values = (float*) calloc(${sub_n*omp_n_threads}, sizeof(float));
    float subsample_values[${sub_n * omp_n_threads}];
    //printf("subsample_values: %x\n", (unsigned int) subsample_values);
    //unsigned int * const bootstrap_indicies = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    unsigned int bootstrap_indicies[${sub_n*omp_n_threads}];
    //printf("bootstrap_indices: %x\n", (unsigned int) bootstrap_indicies);
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically

    printf("About to begin computation\n");
    #pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
    for (int i = 0; i < ${n_subsamples}; i++) {
        int tid = omp_get_thread_num();
        unsigned int* local_indicies = bootstrap_indicies+(${sub_n}*tid);
        float* const local_values = subsample_values+(${sub_n}*tid);
        float* const local_estimates = bootstrap_estimates+(${n_bootstraps}*tid);
        unsigned int seed = time(NULL);
        subsample_and_load(c_arr, local_values, &seed);
        for (int j=0; j < ${n_bootstraps}; j++) {
            loaded_bootstrap(local_indicies, &seed);
            compute_estimate(local_values, local_indicies, ${sub_n}, local_estimates + j);
	    //printf("Completed compute_estimate for bootstrap number %d, subsample number %d \n", j, i);
        }

        reduce_bootstraps(local_estimates, ${n_bootstraps}, subsample_estimates + i);
    }
    float theta = 0;
    average(subsample_estimates, ${n_subsamples}, &theta);

    printf("Freeing subsample_estimates \n");

    printf("subsample_estimates: %x\n", (unsigned int) subsample_estimates);
    printf("bootstrap_estimates: %x\n", (unsigned int) bootstrap_estimates);
    printf("subsample_values: %x\n", (unsigned int) subsample_values);
    printf("bootstrap_indicies: %x\n", (unsigned int) bootstrap_indicies);

    printf("1");
//    free(subsample_estimates);
//    free(bootstrap_estimates);
//    free(subsample_values);
//    free(bootstrap_indicies);

%if seq_type is UNDEFINED or seq_type == 'list':
    Py_DECREF( py_arr );
%endif
    Py_DECREF( data );
    return PyFloat_FromDouble(theta);
}



