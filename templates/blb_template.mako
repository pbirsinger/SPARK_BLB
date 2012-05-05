
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

#define NPY_SCALAR ${ 'NPY_FLOAT32' if average_model.scalar_t.ctype() is 'float' else 'NPY_FLOAT64' }
void printArray( unsigned int*, int, int );

void bootstrap( unsigned int* out, gsl_rng* rng ){
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
void subsample_and_load( scalar_t* data, scalar_t* out, gsl_rng* rng ){
        for( int i = 0; i<${vec_n}; i++ ){
	    unsigned int index = gsl_rng_get(rng) % ${n_vecs};
            vvc( data+(${dim}*index), out+(i*${dim}), ${dim} );
        }
}
</%doc>

<%
    data_args = [ '%s* data%d' % ( arg_model[i].scalar_t.ctype(), i ) for i in range(len(arg_model) ) ]
    out_args = [ '%s* out%d' % ( arg_model[i].scalar_t.ctype(), i ) for i in range(len(arg_model)) ]
%> 
void subsample_and_load( ${ ', '.join( data_args + out_args  )}, gsl_rng* rng ){
        for( int i = 0; i<${vec_n}; i++ ){
	    unsigned int index = gsl_rng_get(rng) % ${n_vecs};
%for i in range(len(arg_model)):
    %if len(arg_model[i].dimensions) == 1:
	    out${i}[i] = data${i}[index];
    %else:
        <% dim = arg_model[i].dimensions[1]  %>
	memcpy(out${i} + (i*${dim}), data${i} + index*${dim}, ${dim}*sizeof(${arg_model[i].scalar_t.ctype()}) );
    %endif
%endfor
        }
}

int check_nans( double* data, int n ){
    for( int i = 0; i<n; i++ ){
	if( data[i] != data[i] ) return 1;
    }
    return 0;
}
## list is the default type.
%if seq_type == 'list':

PyObject* compute_blb( PyObject*  data ){
    Py_INCREF(data);
    PyObject * py_arr = PyArray_FROM_OTF( data, NPY_SCALAR, NPY_IN_ARRAY );
    Py_INCREF( py_arr );
    scalar_t * c_arr = (scalar_t*) PyArray_DATA( py_arr );

%elif seq_type == 'ndarray':
PyObject* compute_blb( ${ ', '.join( [ 'PyObject* arg%d' % i for i in range(len(arg_model)) ] ) } ){
//  printf("Made it into C\n");
%for i in range(len(arg_model)):
  <% argname = 'arg' + str(i) %>	
  Py_INCREF( ${argname} );
  ${arg_model[i].scalar_t.ctype()} * c_arr${i} = (${arg_model[i].scalar_t.ctype()}*) PyArray_DATA( ${argname} );

##  printf("c_arr${i} nominal length: ${reduce(int.__mul__, arg_model[i].dimensions )}\n");
##  printf("actual length: %d\n", PyArray_Size( ${argname} ) );
%endfor
%endif
    //note that these are never cleared; We always fill them up
    //with the appropriate data before perform calculations on them.
%for i in range(len(arg_model)):
    <% 
	model = arg_model[i]
        scalar_t = model.scalar_t.ctype() 
        dim = model.element_size()
    %>
    ${scalar_t} * subsample_values${i} = (${scalar_t}*) calloc(${vec_n*dim}, sizeof(${scalar_t}));
%endfor

    <%  
	sub_type = subsample_model.scalar_t.ctype() 
	boot_type = bootstrap_model.scalar_t.ctype()
    %>	

    ${sub_type}* subsample_estimates = (${sub_type}*) malloc(${n_subsamples*subsample_model.dimension()}*sizeof(${sub_type}));
    ${boot_type}* bootstrap_estimates = (${boot_type}*) malloc(${n_bootstraps*bootstrap_model.dimension()}*sizeof(${boot_type}));
    unsigned int* bootstrap_weights = (unsigned int*) malloc(${vec_n}*sizeof(unsigned int)); 

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set( rng, time(NULL) );

    for( int i=0; i<${n_subsamples}; i++ ){
	subsample_and_load( ${ ', '.join(['c_arr'+str(i) for i in range(len(arg_model))] + [ 'subsample_values'+str(i) for i in range(len(arg_model)) ]) }, rng);
        for( int j=0; j<${n_bootstraps}; j++ ){
            bootstrap( bootstrap_weights, rng );
            compute_estimate( ${ ', '.join( [ 'c_arr'+str(i) for i in range(len(arg_model)) ] )}, bootstrap_weights, bootstrap_estimates+j*${bootstrap_model.dimension()} );
	    if( check_nans( bootstrap_estimates+j*${bootstrap_model.dimension()}, ${bootstrap_model.dimension()} ) ){
		printf("nan detected in subsample %d, bootstrap %d!\n", i, j );
	    }
}
    reduce_bootstraps( bootstrap_estimates, subsample_estimates + i*${subsample_model.dimension()} );
}
  
    gsl_rng_free(rng);

    %if average_model.dimension() == 1:
    ${average_model.scalar_t.ctype()} theta = (${average_model.scalar_t.ctype()})0;
    average(subsample_estimates, &theta ); 
    %else:
        <% scalar_t = average_model.scalar_t.ctype() %>
    ${scalar_t}* theta = (${scalar_t}*) malloc( ${average_model.dimension()}*sizeof(${scalar_t}) );
    average( subsample_estimates, theta );
    %endif

%for i in range(len(arg_model)):
    free(subsample_values${i});
%endfor
    free(subsample_estimates);
    free(bootstrap_estimates);
    free(bootstrap_weights);
##<% assert False, 'made it through free' %>
%if seq_type is UNDEFINED or seq_type == 'list':
    Py_DECREF( py_arr );
%endif
%for i in range(len(arg_model)):
    Py_DECREF( arg${i} );
%endfor
##<% assert False, 'made it through decrefs' %>
    %if average_model.dimension() == 1:
    return PyFloat_FromDouble( theta );
    %else:
    npy_intp dim[1] = { ${average_model.dimension()} };
    return PyArray_SimpleNewFromData( 1, dim, NPY_SCALAR, theta );
    %endif
}

void printArray(unsigned int* arr, int start, int end) {
     for (int i=start; i<end; i++) {
         printf("%d, ", arr[i]);
     }
     printf("\n");
}
