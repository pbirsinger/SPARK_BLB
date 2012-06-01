<%doc>
Framework template for OpenMP implementation of BLB algorithm.
see 'blb_template.mako' for documentation on templating variables.
</%doc>
typedef ${scalar_type} scalar_t;
#define NPY_SCALAR ${ 'NPY_FLOAT32' if scalar_type is 'float' else 'NPY_FLOAT64' }

##', '.join(['c_arr'+str(i) for i in len(range(arg_model))] + [ 'local_values'+str(i) for i in len(range(arg_model)) ])


//New-style subsample
<% 
    s_model = filter(lambda x: x.should_subsample, arg_model)
    data_args = [ '%s* %s_data' % ( arg.scalar_t.ctype(), arg.ref_name() ) for arg in s_model ]
    out_args = [ '%s* %s_out' % ( arg.ctype(), arg.ref_name() ) for arg in s_model ]
%>
void subsample_and_load( ${ ', '.join( data_args + out_args ) }, gsl_rng* rng ){
    for( int i = 0; i<${vec_n}; i++ ){
	unsigned int index = gsl_rng_get(rng) % ${n_vecs};
%for arg in s_model:
    <% name, dim = arg.ref_name(), arg.element_size() %>
	${arg.subsample_method()}
##    %if arg.is_scalar():
##	${name}_out[i] = ${name}_data[index];
##    %else:
##	memcpy(${name}_out + (i*${dim}), ${name}_data + index*${dim}, ${dim}*sizeof(${arg.scalar_t.ctype()}) );
 ##   %endif
%endfor
    }
}

void bootstrap( unsigned int* out, gsl_rng * rng ){
  double norm = ${ vec_n*(1.0/vec_n) };
  double sum_p = 0.0;
  double p = ${ 1.0/vec_n };
  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (int k = 0; k < ${vec_n}; k++)
    {
     out[k] = gsl_ran_binomial (rng, p / (norm - sum_p), ${n_vecs} - sum_n);
     sum_p += p;
     sum_n += out[k];
    }
}

<%doc>
determine header
allocate temporaries
effect computation
free temporaries
write to output
</%doc>

## list is the default type.
%if seq_type is UNDEFINED or seq_type == 'list':

PyObject * compute_blb(PyObject * data) {
    Py_INCREF( data );
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_SCALAR, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( py_arr );
%elif seq_type == 'ndarray':

PyObject* compute_blb( ${ ', '.join( [ 'PyObject* %s' % arg.ref_name() for arg in arg_model ] ) } ){
//  printf("Made it into C\n");
%for arg in arg_model:
Py_INCREF( ${arg.ref_name()} );
${arg.scalar_t.ctype()}* c_${arg.ref_name()} = (${arg.scalar_t.ctype()}*) PyArray_DATA( ${arg.ref_name()} );
%endfor

%endif    


%for arg in s_model:
##<% scalar_t = arg.scalar_t.ctype() %>
  <% scalar_t = arg.ctype() %>
      ${scalar_t}* const ${arg.ref_name()}_subsamples = (${scalar_t}*) calloc(${vec_n*omp_n_threads*arg.element_size()}, sizeof(${scalar_t}));
%endfor

    <%  sub_type, boot_type = subsample_model.scalar_t.ctype(), bootstrap_model.scalar_t.ctype() %> 
    ${sub_type} * const subsample_estimates = (${sub_type}*) calloc(${n_subsamples*subsample_model.dimension()}, sizeof(${sub_type}));
    ${boot_type} * const bootstrap_estimates = (${boot_type}*) calloc(${n_bootstraps*omp_n_threads*bootstrap_model.dimension()}, sizeof(${boot_type}));
    unsigned int* const bootstrap_weights = (unsigned int*) calloc(${vec_n*omp_n_threads}, sizeof(unsigned int)); 
##<% assert False, 'this happened here' %>

%if parallel_loop is UNDEFINED or parallel_loop == 'outer':
    srand( time(NULL) );
    gsl_rng** rngs = (gsl_rng**) malloc( ${omp_n_threads}*sizeof(gsl_rng*) );
    for(int i = 0; i<${omp_n_threads}; i++){
	rngs[i] = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngs[i], rand());
    }
    #pragma omp parallel for schedule(static) num_threads(${omp_n_threads})
    for (int i = 0; i < ${n_subsamples}; i++) {
	printf("Starting subsample %d\n", i );
	int tid = omp_get_thread_num();
        unsigned int* local_weights = bootstrap_weights+(${vec_n}*tid);

%for arg in s_model:
	${arg.ctype()}* local_${arg.ref_name()} = ${arg.ref_name()}_subsamples+(${vec_n*arg.element_size()}*tid);
%endfor
        ${bootstrap_model.scalar_t.ctype()}* local_estimates = bootstrap_estimates+(${n_bootstraps*bootstrap_model.dimension()}*tid);
	gsl_rng* rng = rngs[tid];
	subsample_and_load( ${ ', '.join(['c_'+arg.ref_name() for arg in s_model] + ['local_'+arg.ref_name() for arg in s_model])}, rng);
	for (int j=0; j < ${n_bootstraps}; j++) {
	    printf("Starting bootstrap %d\n", j );
            bootstrap(local_weights, rng );
            compute_estimate( ${ ', '.join( [ ('local_' if arg.should_subsample else 'c_')+arg.ref_name() for arg in arg_model ] )}, local_weights, local_estimates+(j*${bootstrap_model.dimension()}) );
        }
        reduce_bootstraps(local_estimates, subsample_estimates+(i*${subsample_model.dimension()}) );
	printf("Done with subsample %d\n", i);
    }
    for( int i=0; i<${omp_n_threads};i++ ){
	gsl_rng_free(rngs[i]);
    }
    free(rngs);
%endif
    %if average_model.dimension() == 1:
    ${average_model.scalar_t.ctype()} theta = (${average_model.scalar_t.ctype()})0;
    average(subsample_estimates, &theta ); 
    %else:
    <% scalar_t = average_model.scalar_t.ctype() %>
    ${scalar_t}* theta = (${scalar_t}*) calloc( ${average_model.dimension()}, sizeof(${scalar_t}) );
    average( subsample_estimates, theta );
    %endif

%for arg in s_model:
free(${arg.ref_name()}_subsamples);
%endfor
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_weights);

%if seq_type is UNDEFINED or seq_type == 'list':
    Py_DECREF( py_arr );
%endif
%for arg in arg_model:
    Py_DECREF( ${arg.ref_name()} );
%endfor
    %if average_model.dimension() == 1:
    return PyFloat_FromDouble( theta );
    %else:
    npy_intp dim[1] = { ${average_model.dimension()} };
    return PyArray_SimpleNewFromData( 1, dim, NPY_SCALAR, theta );
    %endif
}



