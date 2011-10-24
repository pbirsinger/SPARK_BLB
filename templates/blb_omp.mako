<%doc>
Framework template for OpenMP implementation of BLB algorithm.
see 'blb_template.mako' for documentation on templating variables.
</%doc>


## Always returns a complete statement.
<%def name="rng_init(rng, varname)">
    %if rng == 'gsl_taus':
        gsl_rng * ${varname} = gsl_rng_alloc(gsl_rng_taus);
    %endif
</%def>

<%def name="rng_gen(rng, varname, size)">   
    %if rng is UNDEFINED or rng == 'rand':
        (rand() % ${size})
    %elif rng == 'gsl_taus':
        gsl_rng_uniform_unt(${varname}, ${size})
    %endif
</%def>

#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)
#define BIT_OFFSET(b)  ((b) % BITS_PER_WORD)

typedef uint32_t word_t;
enum {
	BITS_PER_WORD = sizeof(word_t) * CHAR_BIT
};

//USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
void bootstrap(unsigned int* in, unsigned int* out, gsl_rng* rng) {
	for (int i = 0; i < ${sub_n}; i++) {
		int index = gsl_rng_uniform_int(rng, ${sub_n});
		out[i] = in[index];
	}
}


void set_bit(word_t *words, int n) {
        printf("Setting a bit\n");
	words[WORD_OFFSET(n)] |= (1 << BIT_OFFSET(n));
}

void clear_bit(word_t *words, int n) {
	words[WORD_OFFSET(n)] &= ~(1 << BIT_OFFSET(n));
}

int get_bit(word_t *words, int n) {
        printf("Getting a bit\n");
	word_t bit = words[WORD_OFFSET(n)] & (1 << BIT_OFFSET(n));
	return bit != 0;
}

<%doc>
We're going to try a little experiment.
We're goiing to play it fast and loose, relying upon the period
of our RNG to not select the same index twice in short sequence.
</%doc>
void subsample(unsigned int* out, gsl_rng* rng) {
     	unsigned int index = 0;
	for (int i=0; i < ${sub_n}; i++) {
		index =  gsl_rng_uniform_int(rng, ${n_data}); //Generate random integer between 0 and size_in
		out[i] = index;
	}
}

PyObject * compute_blb(PyObject * data) {
    Py_INCREF( data );
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_FLOAT32, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    float * c_arr = (float*) PyArray_DATA( py_arr );
    
    float * subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float * bootstrap_estimates = (float*) calloc(${n_bootstraps*omp_n_threads}, sizeof(float));
    unsigned int * subsample_indicies = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    unsigned int * bootstrap_indicies = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically

    #pragma omp parallel num_threads(${omp_n_threads})
{
    int tid = omp_get_thread_num();
    unsigned int* local_bi = bootstrap_indicies+(${sub_n}*tid);
    unsigned int* local_si = subsample_indicies+(${sub_n}*tid);
    float* local_be = bootstrap_estimates+(${n_bootstraps}*tid);
    gsl_rng * m_rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(m_rng, time(NULL)*(1+tid)); 

    #pragma omp for schedule(static) 
    for (int i = 0; i < ${n_subsamples}; i++) {
    	subsample(local_si, m_rng);
	for (int j=0; j < ${n_bootstraps}; j++) {
	    bootstrap(local_si, local_bi, m_rng);
	    local_be[j] = compute_estimate(c_arr, local_bi, ${sub_n});
	}

	subsample_estimates[i] = reduce_bootstraps(local_be, ${n_bootstraps});
    }
    //printf("About to free m_rng for thread %d\n", tid);
    gsl_rng_free(m_rng);
    //printf("Freed m_rng for thread %d\n", tid);    
}//end parallel
    float theta = average(subsample_estimates, ${n_subsamples});
    //printf("About to free subsample_indicies\n");
    free(subsample_indicies);
    //printf("About to free subsample_estimates\n");
    free(subsample_estimates);
    //printf("About to free bootstrap_estimates\n");
    free(bootstrap_estimates);
    //printf("About to free bootstrap_indicies\n");
    free(bootstrap_indicies);

    //printf("About to DECREF\n");
    Py_DECREF( py_arr );
    Py_DECREF( data );
    return PyFloat_FromDouble(theta);
}
