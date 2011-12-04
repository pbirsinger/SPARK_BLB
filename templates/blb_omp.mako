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

void subsample_and_load( float* data, float* out, gsl_rng* rng  ){
        for( int i = 0; i<${sub_n}; i++ ){
            out[i] = data[ gsl_rng_get(rng) % ${n_data} ];
        }
}

void loaded_bootstrap( unsigned int* out, gsl_rng * rng ){
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

void permutation_bootstrap( unsigned int* out, unsigned int* seed){
    unsigned int index = rand_r(seed) % ${sub_n};
    unsigned int start = out[ index ];
    unsigned int carry = start;
    unsigned int temp;
    for( int i = 1; i< ${sub_n}; i++){
	index = rand_r(seed) % ${sub_n};
	temp = out[ index ];
	out[ index ] = carry;
	carry = temp;
    }
    out[ index ] += (temp - start);

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
//  printf("Made it into C\n");
  Py_INCREF( data );
  float * c_arr = (float*) PyArray_DATA( data );
%endif    
    unsigned int start = time(NULL);
    printf("Start of timing.\n");
    srand(start);    

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
%if parallel_loop is UNDEFINED or parallel_loop == 'outer':
    float * subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float * bootstrap_estimates = (float*) calloc(${n_bootstraps*omp_n_threads}, sizeof(float));
    float * subsample_values = (float*) calloc(${sub_n*omp_n_threads}, sizeof(float));
    unsigned int * bootstrap_weights = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically
    gsl_rng** rngs = (gsl_rng**) malloc( ${omp_n_threads}*sizeof(gsl_rng*));

    for(int i = 0; i<${omp_n_threads}; i++){
	rngs[i] = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngs[i], rand());
    }
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
	if( !tid )
	    printf("Thread %d starting subsample %d, time ellapsed %d\n", tid, i, time(NULL) - start);
        unsigned int* local_weights = bootstrap_weights+(${sub_n}*tid);
        float* local_values = subsample_values+(${sub_n}*tid);
        float* local_estimates = bootstrap_estimates+(${n_bootstraps}*tid);
	gsl_rng* rng = rngs[tid];
	unsigned int bootstrap_seed = (unsigned int) gsl_rng_get(rng);
        subsample_and_load(c_arr, local_values, rng);
	loaded_bootstrap( local_weights, rng );
	if( !tid )
	    printf("Thread %d preparation for subsample %d done, time ellapsed %d\n", tid, i, time(NULL) - start);
        for (int j=0; j < ${n_bootstraps}; j++) {
            permutation_bootstrap(local_weights, &bootstrap_seed);
            local_estimates[j] = compute_estimate(local_values, local_weights, ${sub_n});
        }
	if( !tid )
	    printf("Thread %d bootstraps accomplished for subsample %d, time ellapsed %d\n", tid, i, time(NULL) - start);
        subsample_estimates[i] = reduce_bootstraps(local_estimates, ${n_bootstraps});
	if( !tid )
	     printf("Thread %d finished subsample %d, time ellapsed %d\n",tid, i, time(NULL) - start);
    }

<<<<<<< HEAD
        reduce_bootstraps(local_estimates, ${n_bootstraps}, subsample_estimates + i);
=======
    for( int i=0; i<${omp_n_threads};i++ ){
	gsl_rng_free(rngs[i]);
    }
    free(rngs);
%elif parallel_loop == 'inner':
    float * subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float * bootstrap_estimates = (float*) calloc(${n_bootstraps}, sizeof(float));
    float * subsample_values = (float*) calloc(${sub_n}, sizeof(float));
    unsigned int * bootstrap_weights = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    unsigned int * seeds = (unsigned int *) calloc( ${omp_n_threads}, sizeof(unsigned int));
    //seed the RNGS a bit
    for( int i = 0; i<${omp_n_threads}; i++ ){
	seeds[i] = rand() % ${sub_n};
    }

    for( int i = 0; i< ${n_subsamples}; i++ ){
	subsample_and_load( c_arr, subsample_values, seeds );
	# pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
	for( int j = 0; j< ${n_bootstraps}; j++ ){
	    int tid = omp_get_thread_num();
	    unsigned int* local_weights = bootstrap_weights+(${sub_n}*tid);
	    loaded_bootstrap(local_weights, seeds+tid);
	    bootstrap_estimates[j] = compute_estimate(subsample_values, local_weights, ${sub_n});
	}
	subsample_estimates[i] = reduce_bootstraps( bootstrap_estimates, ${n_bootstraps} );
    }
    free( seeds );
%elif parallel_loop == 'teams':
    float * subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float * bootstrap_estimates = (float*) calloc(${n_bootstraps*omp_n_teams}, sizeof(float));
    float * subsample_values = (float*) calloc(${sub_n*omp_n_teams}, sizeof(float));
    unsigned int * bootstrap_weights = (unsigned int*) calloc(${sub_n*omp_n_teams*omp_team_size}, sizeof(unsigned int));
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically
    gsl_rng** rngs = (gsl_rng**) malloc( ${omp_n_teams*omp_team_size}*sizeof(gsl_rng*));

    for(int i = 0; i<${omp_n_teams*omp_team_size}; i++){
	rngs[i] = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngs[i], rand());
    }
    #pragma omp parallel for schedule(static) num_threads(${omp_n_teams*omp_team_size})
    for( int i = 0; i< ${omp_n_teams*omp_team_size}; i++ ){
	unsigned int* local_weights = bootstrap_weights+(i*${sub_n});
	loaded_bootstrap( local_weights, rngs[i] );
    }
    #pragma omp parallel for schedule(static) num_threads(${omp_n_teams})
    for (int i = 0; i < ${n_subsamples}; i++) {
        int tid = omp_get_thread_num();
	if( !tid ) printf("Thread %d starting subsample %d, time ellapsed %d\n", tid, i, time(NULL) - start);
        float* local_values = subsample_values+(${sub_n}*tid);
        float* local_estimates = bootstrap_estimates+(${n_bootstraps}*tid);
	gsl_rng* rng = rngs[tid];
        subsample_and_load(c_arr, local_values, rng);
	if( !tid ) printf("Thread %d preparation for subsample %d done, time ellapsed %d\n", tid, i, time(NULL) - start);
	#pragma omp parallel for schedule(static) num_threads(${omp_team_size})
        for (int j=0; j < ${n_bootstraps}; j++) {
	    int stid = omp_get_thread_num();
	    unsigned int seed = time(NULL)*stid;
	    unsigned int * local_weights = bootstrap_weights+(${sub_n}*stid);
	    permutation_bootstrap( local_weights, &seed  );
	    local_estimates[j] = compute_estimate(local_values, local_weights, ${sub_n});
        }
	if( !tid ) printf("Thread %d bootstraps accomplished for subsample %d, time ellapsed %d\n", tid, i, time(NULL) - start);
        subsample_estimates[i] = reduce_bootstraps(local_estimates, ${n_bootstraps});
	if( !tid ) printf("Thread %d finished subsample %d, time ellapsed %d\n",tid, i, time(NULL) - start);
>>>>>>> 895e8963812d9de0345657e9377a5a5f2445a8f8
    }
    float theta = 0;
    average(subsample_estimates, ${n_subsamples}, &theta);

    printf("Freeing subsample_estimates \n");

    printf("subsample_estimates: %x\n", (unsigned int) subsample_estimates);
    printf("bootstrap_estimates: %x\n", (unsigned int) bootstrap_estimates);
    printf("subsample_values: %x\n", (unsigned int) subsample_values);
    printf("bootstrap_indicies: %x\n", (unsigned int) bootstrap_indicies);

    for( int i=0; i<${omp_n_teams};i++ ){
	gsl_rng_free(rngs[i]);
    }
    free(rngs);
%endif
    float theta = average(subsample_estimates, ${n_subsamples}); 
    free(subsample_values);
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_weights);

%if seq_type is UNDEFINED or seq_type == 'list':
    Py_DECREF( py_arr );
%endif
    Py_DECREF( data );
    return PyFloat_FromDouble(theta);
}



