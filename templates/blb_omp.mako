<%doc>
Framework template for OpenMP implementation of BLB algorithm.
see 'blb_template.mako' for documentation on templating variables.
</%doc>


#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)
#define BIT_OFFSET(b)  ((b) % BITS_PER_WORD)

typedef uint32_t word_t;
enum {
	BITS_PER_WORD = sizeof(word_t) * CHAR_BIT
};

//USING ARRAYS OF INDICIES INSTEAD OF COPPYING DATA
void bootstrap(unsigned int* in, unsigned int* out) {
	int i = 0;
	for (i; i < ${sub_n}; i++) {
		int index = rand() % ${sub_n};
		out[i] = in[index];
	}
}


void set_bit(word_t *words, int n) {
	words[WORD_OFFSET(n)] |= (1 << BIT_OFFSET(n));
}

void clear_bit(word_t *words, int n) {
	words[WORD_OFFSET(n)] &= ~(1 << BIT_OFFSET(n));
}

int get_bit(word_t *words, int n) {
	word_t bit = words[WORD_OFFSET(n)] & (1 << BIT_OFFSET(n));
	return bit != 0;
}

void subsample(unsigned int* out) {
     	//gsl_rng * r = gsl_rng_alloc(gsl_rng_taus); //Initialize random number generator
	//gsl_rng_set(r, 0); //Initialize with standard seed
	int i = 0;
	unsigned int index = 0;
	//Create a bit array to keep track of whether we have seen an index
	word_t * bitArray = (word_t *) malloc( (((int) ${sub_n} / BITS_PER_WORD) + 1) * sizeof(word_t));
	for (i; i < ${sub_n}; i++) {
		do {
			index =  rand() % ${n_data};//gsl_rng_uniform_unt(r, size_in); //Generate random integer between 0 and size_in
		} while (get_bit(bitArray, i));
		//We have not seen index before so we update our bitArray
		set_bit(bitArray, i);
		out[i] = index;
	}
	//Free up memory
	free(bitArray);
//	gsl_rng_free(r);
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
    unsigned int * bootstrap_indicies = (unsigned int*) calloc(${sub_n}, sizeof(unsigned int));
    //We use static scheduling to avoid the overhead of synchronization and assigning tasks to threads dynamically
    #pragma omp parallel num_threads(${omp_n_threads})
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < ${n_subsamples}; i++) {
    	subsample(subsample_indicies);
	int j = 0;
	for (j; j < ${n_bootstraps}; j++) {
	    bootstrap(subsample_indicies, bootstrap_indicies);
	    bootstrap_estimates[j] = compute_estimate(c_arr, bootstrap_indicies, ${sub_n});
	}
	subsample_estimates[i] = reduce_bootstraps(bootstrap_estimates, ${n_bootstraps});
    }
    
    float theta = average(subsample_estimates, ${n_subsamples});
    free(subsample_indicies);
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_indicies);
    Py_DECREF( py_arr );
    Py_DECREF( data );
    return PyFloat_FromDouble(theta);
}
