<%doc>
  This is a dummy impl file which allows the BLB
  C framework code to compile when user-defined
  reducer functions aren't supplied. A substitute
  for this file should generally be provided in real
  JIT compilation.
 
  @author David Howard
 
Templating varriables in use
 sub_n: The size of subsampled data
 n_bootstraps: the number of bootstraps to preform per subsample
 n_subsamples: the numer of subsamples to perform
 bootstrap_sreducer: the name of the function with which to reduce bootstraps
 subsample_reducer: the name of the function with which to reduce subsamples
 classifier: the name of the function to calculate on the bootstraps
 
 attributes: a dictionary of optional flags that are passed into the classifier constructors
  with_cilk: whether or not the rendered functions should make use of cilk
 
Reducers defined by this file follow the interfaces

void compute_estimate( scalar_t* data, unsigned int* indicies, unsigned int size, scalar_t* const result );

scalar_t* reduce_bootstraps( scalar_t* data, unsigned int size );

scalar_t* average( scalar_t * data, unsigned int size );

</%doc>

<%def name="noop( weighted, input_dim, output_dim )" >
</%def>

<%def name="linreg( weighted, input_dim, output_dim )" >
// Form X'X in temporary memory
//printf("Checking data...\n");
//for(int i = 0; i<size; i++ ){
//printf("Checking vector %d\n", i );	
// double * vec = data+i*${input_dim};
// vprint( vec, ${input_dim} );
// printf("\n");	
//}
double _XX[ ${(input_dim-1)*(input_dim-1)} ];
memset( _XX, 0, ${(input_dim-1)*(input_dim-1)}*sizeof( double ) );
_gsl_matrix_view XX = gsl_matrix_view_array( _XX, ${input_dim-1}, ${input_dim-1} );
    for( unsigned int i = 0; i<size; i++ ){
	double* vec = data+i*${input_dim};
        _gsl_vector_const_view v = gsl_vector_const_view_array( vec+1, ${input_dim-1} );
        %if weighted:
        gsl_blas_dger( (double)weights[i]*weights[i], & v.vector, & v.vector, & XX.matrix );
        %else:
        gsl_blas_dger( 1.0, & v.vector, & v.vector, & XX.matrix );
        %endif
    }
// Form X'*y in temporary memory
//printf("About to allocate vector\n");
double _Xy[${input_dim-1}];
memset( _Xy, 0.0, ${input_dim-1}*sizeof(double) );
_gsl_vector_view Xy = gsl_vector_view_array( _Xy, ${input_dim-1} );
//printf("About to calculate vector\n");
    for( unsigned int i = 0; i<size; i++ ){
	double* vec = data + i*${input_dim};
        _gsl_vector_const_view v = gsl_vector_const_view_array( vec+1, ${input_dim-1} );
        %if weighted:
        gsl_blas_daxpy( weights[i]*weights[i]*vec[0], & v.vector, & Xy.vector );
        %else:
        gsl_blas_daxpy( vec[0], & v.vector, & Xy.vector );
        %endif
    }
// Solve X'X*x = X'*y , store in result
//printf("About to calculate result \n");
gsl_linalg_HH_svx( & XX.matrix, & Xy.vector );
//printf("about to copy into result\n");
memcpy( result , Xy.vector.data , ${input_dim-1}*sizeof(double) );
</%def>

<%def name="noop( weighted, input_dim, output_dim )" >
</%def>

##Spits out the body of a mean norm calculation, sans declaration or return statement.
<%def name="mean_norm(weighted, input_dim, output_dim )">
   <%
        access = 'data + i*%s ' % input_dim
	cardinality = 'DATA_SIZE' if weighted else 'size'
   %>
   %if input_dim == 1:
      double mean = 0.0;
      for (unsigned int i=0; i<size; i++) {
           mean += *(${access});
      }
      mean /= ${cardinality};
      result[0] = abs( mean ); //norm of scalar = absolute value
   %else:
    double mean_vec[${input_dim}];
    %if weighted:
    vsp( data, weights[0], mean_vec, ${input_dim} );
    %else:
    vvc( data, mean_vec, ${input_dim} );
    %endif
    for (unsigned int i=1; i<size; i++) {
        %if weighted:
	vspa( ${access}, weights[i], ${input_dim}, mean_vec );
	%else:
	vva( ${access}, mean_vec, ${input_dim} );
	%endif	
    }
    vsid( mean_vec, ${cardinality}, ${input_dim} );
    //Take the norm of the mean vector
    *result =  norm( mean_vec, ${input_dim} );
    
    %endif
</%def>


##Spits out the body of a mean calculation, sans declaration or return statement.
<%def name="mean(weighted, input_dim, output_dim )">
    <%
	access = 'data + i*%s ' % input_dim
	cardinality = 'DATA_SIZE' if weighted else 'size'
    %>
    %if input_dim == 1:
    double mean = 0.0;
    for (unsigned int i=0; i<size; i++) {
    	printf("%f ", ${access});
        mean += *(${access});
    }
    mean /= ${cardinality};
    result[0] = mean;
    %else:
    	%if weighted:
    vsp( data, weights[0], result, ${input_dim} );
    	%else:
    vvc( data, result, ${input_dim} );
    	%endif
    for (unsigned int i=1; i<size; i++) {
        double* vec = ${access};
    	%if weighted:
	vspa( vec, weights[i], ${input_dim}, result );
	%else:
	vva( vec, result, ${input_dim} );
	%endif
    }
    vsid( result, ${cardinality}, ${input_dim} );
    %endif
</%def>


#define LINE_SIZE 64
#define MIN(a,b) ((a>b)?(b):(a))
#define DATA_SIZE ${n_vecs}
#define SUBSAMPLE_SIZE ${sub_n}

//vector print
void vprint( double* const vec, const unsigned int dim ){
   for( unsigned int i = 0; i<dim; i++ ){
	printf("%f ", vec[i]);
   } 
}
// vector-scalar product
inline double* vsp( double* const vec, const int a, double* out, const unsigned int dim ){
    for( unsigned int i = 0; i<dim; i++ ){
        out[i] = vec[i]*a;
    } 
    return out;
}
// vector-vector copy
inline double* vvc( double* const vec, double* out, const unsigned int dim ){
   memcpy( out, vec, dim*sizeof(double));
   return out;
}
// vector-scalar product & add
inline double* vspa( double* vec, int a, unsigned int dim, double* out ){
    for( unsigned i = 0; i<dim; i++ ){
	out[i] += vec[i]*a;
    }
    return out;
}
// vector-vector add
inline double* vva( double* vec, double* out, unsigned int dim ){
    for( unsigned int i = 0; i<dim; i++ ){
	out[i] = out[i] + vec[i];
    }
    return out;
}
// vector-scalar in place divide
inline double* vsid( double* vec, const int a, const unsigned int dim ){
    for( unsigned int i = 0; i<dim; i++ ){
	vec[i] /= a;
    }
    return vec;
}
// l-2 norm
inline double norm( double* vec, unsigned int dim ){
    double  norm = 0.0;
    for( unsigned int i = 0; i<dim; i++ ){
	norm += vec[i]*vec[i];
    }
    return sqrt( norm );	 
}
// weighted vector variance and add
inline double* wvvara( double* vec, unsigned int a, double* mean, double* out, unsigned int dim ){
	for( unsigned int i = 0; i<dim; i++ ){
	    double residual = vec[i] - mean[i];
	    out[i] += a*residual*residual;
	}
	return out;
}
// vector variance and add
inline double* vvara( double* vec, double* mean, double* out, unsigned int dim ){
	for(unsigned int i = 0; i<dim; i++ ){
	    double residual = vec[i] - mean[i];
	    out[i] += residual*residual;
	}
	return out;
}
// vector vector in place divide
inline double* vvid( double* vec, double* quot, unsigned int dim ){
	for( unsigned int i = 0; i<dim; i++){
	    vec[i] /= quot[i];
	}
	return vec;
}
// vector sqrt in place
inline double* vsqrti( double* vec, unsigned int dim ){
	for( unsigned int i = 0; i<dim; i++ ){
	    vec[i] = sqrt( vec[i] );
	}
	return vec;
}
inline double update_mean( const double mu1, const double mu2, const unsigned int n1, const unsigned int n2 ){
     double delta = mu2 - mu1;
     return mu1 + (n2*delta)/(n1 + n2);
}
inline double update_var( const double mu1, const double mu2, const double var1, const double var2, const unsigned int n1, const unsigned n2 ){
    int n = n1 + n2;
    double delta = mu2 - mu1;
    return (n1*var1 + n2*var2 + ((n1*n2*delta*delta)/n) )/n;
}
<%def name="std(weighted, input_dim, output_dim)">
    <%
	n = 'sum_weights' if weighted else 'LINE_SIZE'
	n_left = 'sum_weights' if weighted else 'size - k'
	n_cum = 'sum_weights_cum' if weighted else 'k'
    %>
    double mu_carry[${input_dim}]; //Mean of all the data elements seen before current cache block
    double mu_curr[${input_dim}]; //Accumulator for current cache block
    double var_carry[${input_dim}];
    double var_curr[${input_dim}];
    memset( mu_carry, 0, ${input_dim}*sizeof(double) );
    memset( mu_curr, 0, ${input_dim}*sizeof(double) );
    memset( var_carry, 0, ${input_dim}*sizeof(double) );
    memset( var_curr, 0, ${input_dim}*sizeof(double) ); 
    unsigned int k = 0;
    %if weighted:
    unsigned int sum_weights = 0; //For current cache block
    unsigned int sum_weights_cum = 0;
    %endif
    for( unsigned int i = 0; i<size/LINE_SIZE; i++ ){
	for( unsigned int j = k; j<k+LINE_SIZE; j++ ){
	    %if weighted:
	    vspa( data+j*${input_dim}, weights[j], ${input_dim}, mu_curr );
	    sum_weights += weights[j];
	    %else:
	    vva( data+j*${input_dim}, mu_curr, ${input_dim} );
	    %endif
	}
	vsid( mu_curr, ${n}, ${input_dim} );
	for( unsigned int j= k; j<k+LINE_SIZE; j++ ){
	    %if weighted:
	    wvvara( data+j*${input_dim}, weights[j], mu_curr, var_curr, ${input_dim} );
	    %else:
	    vvara( data+j*${input_dim}, mu_curr, var_curr, ${input_dim} );
	    %endif
	}
	vsid( var_curr, ${n}, ${input_dim} );

	for ( unsigned int l=0; l<${input_dim}; l++) {
	    var_carry[l] = update_var( mu_carry[l], mu_curr[l], var_carry[l], var_curr[l], ${n_cum}, ${n} );
	    mu_carry[l] = update_mean( mu_carry[l], mu_curr[l], ${n_cum}, ${n} );
	}
	k += LINE_SIZE; 
	memset( var_curr, 0, ${input_dim}*sizeof(double));
	memset( mu_curr, 0, ${input_dim}*sizeof(double));
	%if weighted:
	sum_weights_cum += sum_weights;
	sum_weights = 0;
	%endif

    }
    // The leftovers
    for( unsigned int j=k; j<size; j++ ){
	%if weighted:
	vspa( data+j*${input_dim}, weights[j], ${input_dim}, mu_curr );
	sum_weights += weights[j];
	%else:
	vva( data+j*${input_dim}, mu_curr, ${input_dim} );
	%endif
    }
    vsid( mu_curr, ${n_left}, ${input_dim} );
    for( unsigned int j = k; j<size; j++ ){
	%if weighted:
	wvvara( data+j*${input_dim}, weights[j], mu_curr, var_curr, ${input_dim} );
	%else:
	vvara( data+j*${input_dim}, mu_curr, var_curr, ${input_dim} );
	%endif
    }
    vsid( var_curr, ${n_left}, ${input_dim} );

    for ( unsigned int l=0; l<${input_dim}; l++) {
    	var_carry[l] = update_var( mu_carry[l], mu_curr[l], var_carry[l], var_curr[l], ${n_cum}, ${n_left});
    }

    memcpy( result, vsqrti( var_carry, ${input_dim} ), ${input_dim} * sizeof(double) );
</%def>


##produce the classifier from the requested function
<%def name="make_classifier(func_name, input_dim, output_dim)">
    <%
        body = self.template.get_def(func_name).render(True, input_dim, output_dim)
    %>
void compute_estimate( double * const data, unsigned int * const weights, unsigned int size, double* const result ){
      	 ${body}
}
</%def>

##produce the bootstrap reducer from the requested function
<%def name="make_reduce_bootstraps( func_name, input_dim, output_dim )">
    <%
        body = self.template.get_def(func_name).render(False, input_dim, output_dim)
    %>
void reduce_bootstraps( double * const data, unsigned int size, double* const result ){
     	 ${body}
}
</%def>

##produce the subsample reducer from the requested function
<%def name="make_average( func_name, input_dim, output_dim )">
    <%
        body = self.template.get_def(func_name).render(False, input_dim, output_dim)
    %>
void average( double * const data, unsigned int size, double* const result ){
   	 ${ body }
}
</%def>


%for func_desc in desired_funcs:
void ${func_desc[1]}( double* const data, ${ "unsigned int* const weights," if func_desc[2] else "" } const unsigned int size, double* result ){
${self.template.get_def(func_desc[0]).render(func_desc[2], func_desc[3], func_desc[4])}
}
%endfor
%if classifier is not UNDEFINED:
void compute_estimate( double* data, unsigned int* weights,  unsigned int size, double* const result ){
    ${classifier}
}
%elif use_classifier is not UNDEFINED:
    ${make_classifier(use_classifier, dim, bootstrap_dim )}
%endif

%if bootstrap_reducer is not UNDEFINED:
void reduce_bootstraps( double* data, unsigned int size, double* result ){
    ${bootstrap_reducer}
}
%elif use_bootstrap_reducer is not UNDEFINED:
    ${make_reduce_bootstraps(use_bootstrap_reducer, bootstrap_dim, subsample_dim)}
%endif

%if subsample_reducer is not UNDEFINED:
void average( double* data, unsigned int size, double* result ){
    ${subsample_reducer}
}
%elif use_subsample_reducer is not UNDEFINED:
    ${make_average(use_subsample_reducer, subsample_dim, average_dim)}
%endif