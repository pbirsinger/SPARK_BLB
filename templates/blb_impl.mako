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

void compute_estimate( float* data, unsigned int* indicies, unsigned int size, float* const result );

float* reduce_bootstraps( float* data, unsigned int size );

float* average( float * data, unsigned int size );

</%doc>

void printArray(float* arr, int start, int end);
void printArray(unsigned int* arr, int start, int end);

<%def name="noop( weighted, input_dim, output_dim )" >
</%def>

<%def name="linreg(wieghted, input_dim, output_dim)" >
// form X'X

// take X*y 
// solve X'X = X*y with cholskey decomposition
// copy to result
</%def>
##Spits out the body of a mean norm calculation, sans declaration or return statement.
<%def name="mean_norm(weighted, input_dim, output_dim )">
   printf("Mean called with input_dim %d\n", ${input_dim});
   <%
        access = 'data + i*%s ' % input_dim
	cardinality = 'DATA_SIZE' if weighted else 'size'
   %>
   %if input_dim == 1:
      float mean = 0.0;
      for (unsigned int i=0; i<size; i++) {
           mean += *(${access});
      }
      mean /= ${cardinality};
      result[0] = abs( mean ); //norm of scalar = absolute value
      printf("Mean is: %f\n", result[0]);
   %else:
    printf("mean with dimension %d\n", ${input_dim});
    float mean_vec[${input_dim}];
    %if weighted:
    vsp( data, weights[0], mean_vec, ${input_dim} );
    printf("Mean vec is: ");
    printArray(mean_vec, 0, ${input_dim});
    %else:
    vvc( data, mean_vec, ${input_dim} );
    printf("Mean vec is: ");
    printArray(mean_vec, 0, ${input_dim});
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
    printf("Mean vector is: ");
    printArray(mean_vec, 0, ${input_dim});
    *result =  norm( mean_vec, ${input_dim} );
    
    %endif
</%def>


##Spits out the body of a mean calculation, sans declaration or return statement.
<%def name="mean(weighted, input_dim, output_dim )">
    printf("Mean called with input_dim %d\n", ${input_dim});
    <%
	access = 'data + i*%s ' % input_dim
	cardinality = 'DATA_SIZE' if weighted else 'size'
    %>
    %if input_dim == 1:
    float mean = 0.0;
    for (unsigned int i=0; i<size; i++) {
    	printf("%f ", ${access});
        mean += *(${access});
    }
    mean /= ${cardinality};
    result[0] = mean;
    %else:
    	%if weighted:
    vsp( data, weights[0], result, ${input_dim} );
    printf("result of vsp is: ");
    printArray(result, 0, ${input_dim});
    	%else:
    vvc( data, result, ${input_dim} );
    printf("result of vvc is: ");
    printArray(result, 0, ${input_dim});
    	%endif
    for (unsigned int i=1; i<size; i++) {
        float* vec = ${access};
    	%if weighted:
	vspa( vec, weights[i], ${input_dim}, result );
	printf("Result of vspa is: ");
	printArray(result, 0, ${input_dim});
	%else:
	vva( vec, result, ${input_dim} );
	printf("result of vva is: ");
	printArray(result, 0, ${input_dim});
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
void vprint( float* const vec, const unsigned int dim ){
   for( unsigned int i = 0; i<dim; i++ ){
	printf("%f ", vec[i]);
   } 
}
// vector-scalar product
inline float* vsp( float* const vec, const int a, float* out, const unsigned int dim ){
    for( unsigned int i = 0; i<dim; i++ ){
        out[i] = vec[i]*a;
    } 
    return out;
}
// vector-vector copy
inline float* vvc( float* const vec, float* out, const unsigned int dim ){
   memcpy( out, vec, dim*sizeof(float));
   return out;
}
// vector-scalar product & add
inline float* vspa( float* vec, int a, unsigned int dim, float* out ){
    for( unsigned i = 0; i<dim; i++ ){
	out[i] += vec[i]*a;
    }
    return out;
}
// vector-vector add
inline float* vva( float* vec, float* out, unsigned int dim ){
    for( unsigned int i = 0; i<dim; i++ ){
	out[i] = out[i] + vec[i];
    }
    return out;
}
// vector-scalar in place divide
inline float* vsid( float* vec, const int a, const unsigned int dim ){
    for( unsigned int i = 0; i<dim; i++ ){
	vec[i] /= a;
    }
    return vec;
}
// l-2 norm
inline float norm( float* vec, unsigned int dim ){
    float  norm = 0.0;
    for( unsigned int i = 0; i<dim; i++ ){
	norm += vec[i]*vec[i];
    }
    return sqrt( norm );	 
}
// weighted vector variance and add
inline float* wvvara( float* vec, unsigned int a, float* mean, float* out, unsigned int dim ){
	for( unsigned int i = 0; i<dim; i++ ){
	    float residual = vec[i] - mean[i];
	    out[i] += a*residual*residual;
	}
	return out;
}
// vector variance and add
inline float* vvara( float* vec, float* mean, float* out, unsigned int dim ){
	for(unsigned int i = 0; i<dim; i++ ){
	    float residual = vec[i] - mean[i];
	    out[i] += residual*residual;
	}
	return out;
}
// vector vector in place divide
inline float* vvid( float* vec, float* quot, unsigned int dim ){
	for( unsigned int i = 0; i<dim; i++){
	    vec[i] /= quot[i];
	}
	return vec;
}
// vector sqrt in place
inline float* vsqrti( float* vec, unsigned int dim ){
	for( unsigned int i = 0; i<dim; i++ ){
	    vec[i] = sqrt( vec[i] );
	}
	return vec;
}
inline float update_mean( const float mu1, const float mu2, const unsigned int n1, const unsigned int n2 ){
     float delta = mu2 - mu1;
     return mu1 + (n2*delta)/(n1 + n2);
}
inline float update_var( const float mu1, const float mu2, const float var1, const float var2, const unsigned int n1, const unsigned n2 ){
    int n = n1 + n2;
    float delta = mu2 - mu1;
    return (n1*var1 + n2*var2 + ((n1*n2*delta*delta)/n) )/n;
}
<%def name="stdev(weighted, input_dim, output_dim)">
    <%
	n = 'sum_weights' if weighted else 'LINE_SIZE'
	n_left = 'sum_weights' if weighted else 'size - k'
	n_cum = 'sum_weights_cum' if weighted else 'k'
    %>
    float mu_carry[${input_dim}]; //Mean of all the data elements seen before current cache block
    float mu_curr[${input_dim}]; //Accumulator for current cache block
    float var_carry[${input_dim}];
    float var_curr[${input_dim}];
    memset( mu_carry, 0, ${input_dim}*sizeof(float) );
    memset( mu_curr, 0, ${input_dim}*sizeof(float) );
    memset( var_carry, 0, ${input_dim}*sizeof(float) );
    memset( var_curr, 0, ${input_dim}*sizeof(float) ); 
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
	memset( var_curr, 0, ${input_dim}*sizeof(float));
	memset( mu_curr, 0, ${input_dim}*sizeof(float));
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

    memcpy( result, vsqrti( var_carry, ${input_dim} ), ${input_dim} * sizeof(float) );
</%def>


##produce the classifier from the requested function
<%def name="make_classifier(func_name, input_dim, output_dim)">
    <%
        body = self.template.get_def(func_name).render(True, input_dim, output_dim)
    %>
void compute_estimate( float * const data, unsigned int * const weights, unsigned int size, float* const result ){
      	 ${body}
}
</%def>

##produce the bootstrap reducer from the requested function
<%def name="make_reduce_bootstraps( func_name, input_dim, output_dim )">
    <%
        body = self.template.get_def(func_name).render(False, input_dim, output_dim)
    %>
void reduce_bootstraps( float * const data, unsigned int size, float* const result ){
     	 ${body}
}
</%def>

##produce the subsample reducer from the requested function
<%def name="make_average( func_name, input_dim, output_dim )">
    <%
        body = self.template.get_def(func_name).render(False, input_dim, output_dim)
    %>
void average( float * const data, unsigned int size, float* const result ){
   	 ${ body }
}
</%def>

%if classifier is not UNDEFINED:
void compute_estimate( float* data, unsigned int* indicies,  unsigned int size, float* const result ){
    ${classifier}
}
%elif use_classifier is not UNDEFINED:
    ${make_classifier(use_classifier, dim, bootstrap_dim )}
%endif

%if bootstrap_reducer is not UNDEFINED:
float* reduce_bootstraps( float* data, unsigned int size ){
    ${bootstrap_reducer}
%elif use_bootstrap_reducer is not UNDEFINED:
    ${make_reduce_bootstraps(use_bootstrap_reducer, bootstrap_dim, subsample_dim)}
%endif

%if subsample_reducer is not UNDEFINED:
float average( float* data, unsigned int size ){
    ${subsample_reducer}
}
%elif use_subsample_reducer is not UNDEFINED:
    ${make_average(use_subsample_reducer, subsample_dim, average_dim)}
%endif