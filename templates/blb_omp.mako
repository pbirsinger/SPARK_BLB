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
     	//gsl_rng * r = gsl_rng_alloc(gsl_rng_taus); //Initialize random number generator
	//gsl_rng_set(r, time(NULL)); //Initialize with standard seed
	unsigned int index = 0;
	//Create a bit array to keep track of whether we have seen an index
	//char* bitArray = (char*) malloc(  ${n_data} * sizeof(char));
	//memset(bitArray, 0, ${n_data} );
	for (int i=0; i < ${sub_n}; i++) {
		//do {
			index =  gsl_rng_uniform_int(rng, ${n_data}); //Generate random integer between 0 and size_in
		//} while( bitArray[index]);
		//We have not seen index before so we update our bitArray
		//bitArray[index] = 1;
		out[i] = index;
	}
	//Free up memory
	//free(bitArray);
	//gsl_rng_free(r);
}

PyObject * compute_blb(PyObject * data) {
    Py_INCREF( data );
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_FLOAT32, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    float * c_arr = (float*) PyArray_DATA( py_arr );
    srand(time(NULL));
    
    float * subsample_estimates = (float*) calloc(${n_subsamples}, sizeof(float));
    float * bootstrap_estimates = (float*) calloc(${sub_n}, sizeof(float));
    unsigned int * subsample_indicies = (unsigned int*) calloc(${sub_n}, sizeof(unsigned int));
    unsigned int * bootstrap_indicies = (unsigned int*) calloc(${sub_n*omp_n_threads}, sizeof(unsigned int));
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically
    gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    for (int i = 0; i < ${n_subsamples}; i++) {
    	subsample(subsample_indicies, rng);

	#pragma omp parallel num_threads(${omp_n_threads})
	{
	unsigned int* x = bootstrap_indicies+(${sub_n}*omp_get_thread_num());
	gsl_rng * m_rng = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(m_rng, time(NULL)*(1+omp_get_thread_num())); 
	#pragma omp for schedule(static)
	for (int j=0; j < ${n_bootstraps}; j++) {
	    bootstrap(subsample_indicies, x, rng);
	    bootstrap_estimates[j] = compute_estimate(c_arr, x, ${sub_n});
	}
	gsl_rng_free(m_rng);
	} // end parallel
	subsample_estimates[i] = reduce_bootstraps(bootstrap_estimates, ${n_bootstraps});
    }
    
    float theta = average(subsample_estimates, ${n_subsamples});
    free(subsample_indicies);
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_indicies);
    gsl_rng_free(rng);

    Py_DECREF( py_arr );
    Py_DECREF( data );
    return PyFloat_FromDouble(theta);
}
