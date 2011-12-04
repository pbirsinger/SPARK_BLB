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

##Spits out the body of a mean norm calculation, sans declaration or return statement.
<%def name="mean_norm(iter_direct, input_dim, output_dim )">
    <%
        i = 'i' if iter_direct else 'indicies[i]'
        access = 'data + %s*%s ' % (i, input_dim)
   %>
   %if input_dim == 1:
      float mean = 0.0;
      for (unsigned int i=0; i<size; i++) {
          mean += *(${access});
      }
      mean /= size;
      result[0] = sqrt(mean * mean); //norm
   %else:
   float mean_vec[${input_dim}];
    %if iter_direct:
    float* initial=data;
    %else:
    float* initial= data+(indicies[0]*${input_dim});
    %endif
    %for k in xrange(input_dim):
    mean_vec[${k}] = initial[${k}];
    %endfor

    for (unsigned int i=1; i<size; i++) {
         float* vec = ${access};
         %for k in xrange(input_dim):
             mean_vec[${k}] +=  vec[${k}];
         %endfor
    }
    %for j in xrange(input_dim):
         mean_vec[${j}] /= size;
    %endfor
    //Take the norm of the mean vector
    float norm = 0.0;
    printf("Average vector is ");
    %for j in xrange(input_dim):
    	 printf("%f ", mean_vec[${j}]);
    	 norm += mean_vec[${j}] * mean_vec[${j}];
    %endfor
    printf("\n");
    if (norm < 0) {
       printf("negative norm");
    }
    result[0] = sqrt(norm);
    if (result[0] != result[0]) {
       printf("nan detected");
    }
    %endif
</%def>


##Spits out the body of a mean calculation, sans declaration or return statement.
<%def name="mean(iter_direct, input_dim, output_dim )">
    <%
        i = 'i' if iter_direct else 'indicies[i]'
	access = 'data + %s*%s ' % (i, input_dim)
   %>
   %if input_dim == 1:
      float mean = 0.0;
      for (unsigned int i=0; i<size; i++) {
      	  mean += *(${access});
      }
      mean /= size;
      result[0] = mean;
   %else:
    %if iter_direct:
    float* initial=data;
    %else:
    float* initial= data+(indicies[0]*${input_dim});
    %endif
    %for k in xrange(input_dim):
    result[${k}] = initial[${k}];
    %endfor

    for (unsigned int i=1; i<size; i++) {
         float* vec = ${access};
    	 %for k in xrange(input_dim):
             result[${k}] +=  vec[${k}];
	 %endfor
    }
    %for j in xrange(input_dim):
    	 result[${j}] /= size;
    //printf("Mean computed: %f\n", result[0]); 
    %endfor
    %endif

#define LINE_SIZE 64
#define MIN(a,b) ((a>b)?(b):(a))
#define DATA_SIZE ${n_data}
#define SUBSAMPLE_SIZE ${sub_n}

##Spits out the body of a mean calculation, sans declaration or return statement.
<%def name="mean()">
    float mean = 0.0;
    for( unsigned int i=0; i<size; i++ ){
       mean += data[i];
    }			 
    mean /= size;
</%def>
<%def name="weighted_mean()" >
    float mean = 0.0;
    for(unsigned int i=0; i<size; i++){
	mean += data[i]*weights[i];
    }
    mean /= DATA_SIZE;
</%def>
inline float update_mean( const float mu1, const float mu2, const unsigned int n1, const unsigned int n2 ){
     float delta = mu2 - mu1;
     return mu1 + (n2*delta)/(n1 + n2);
}
inline float update_var( const float mu1, const float mu2, const float var1, const float var2, const unsigned int n1, const unsigned n2 ){
    int size = n1 + n2;
    float delta = mu2 - mu1;
    return (n1*var1 + n2*var2 + (n1*delta*n2*delta)/size)/size;
}
<%def name="weighted_stdev()">
    float mu_carry = 0.0;
    float mu_curr = 0.0;
    float var_carry = 0.0;
    float var_curr = 0.0;
    int k = 0;
    unsigned int sum_weights = 0;
    unsigned int sum_weights_cum = 0;
    for( int i= 0; i<size/LINE_SIZE; i++ ){
	for( int j = k; j<k+LINE_SIZE; j++ ){
	    mu_curr += data[j]*weights[j];
	    sum_weights += weights[j];
	}
	mu_curr /= sum_weights;
	for( int j= k; j<k+LINE_SIZE; j++ ){
	    float residual = data[j] - mu_curr;
	    var_curr += weights[j]*residual*residual;
	}
	var_curr /= sum_weights;
	var_carry = update_var( mu_carry, mu_curr, var_carry, var_curr, sum_weights_cum, sum_weights );
	mu_carry = update_mean( mu_carry, mu_curr, sum_weights_cum, sum_weights );
	sum_weights_cum += sum_weights;
	mu_curr =  var_curr = 0.0;
	sum_weights = 0;
	k += LINE_SIZE;
    }
    // The leftovers
    for( int j=k; j<size; j++ ){
	mu_curr += weights[j]*data[j];
	sum_weights += weights[j];
    }
    mu_curr /= sum_weights;
    for( int j = k; j<size; j++ ){
	float residual = data[j] - mu_curr;
	var_curr += weights[j]*residual*residual;
    }
    var_curr /= sum_weights;
    var_carry = update_var( mu_carry, mu_curr, var_carry, var_curr, sum_weights_cum, sum_weights);
    float stdev = sqrt( var_carry );
</%def>
##Spits out the body of a standard deviation calculation, like unto mean defined above.
<%def name="stdev(iter_direct, input_dim, output_dim )">
      <%
	access = 'data+i*%d' %input_dim if iter_direct else 'data + indicies[i]*%d ' %input_dim
      %>
    float mean[${input_dim}];
    %if iter_direct:
    	float* initial=data;
    %else:
        float* initial=data+(indicies[0]*${input_dim});
    %endif
    %for k in xrange(input_dim):
    mean[${k}] = initial[${k}];
    %endfor
    for (unsigned int i=1; i<size; i++) {
         float* vec = ${access};
         %for k in xrange(input_dim):
             mean[${k}] +=  vec[${k}];
         %endfor
    }
    %for j in xrange(input_dim):
         mean[${j}] /= size;
    %endfor
    
    %for l in xrange(input_dim):
    	 result[${l}] = 0.0;
    %endfor

    for( unsigned int i=0; i<size; i++ ){
	 float* vec = ${access};
	 for (unsigned int j=0; j<${input_dim}; j++) {
	     float x = vec[j] - mean[j];
	     x *= x;
	     if (x < 0) {
	     	printf("x*=x makes x negative\n");
	     }
	     result[j] += x;
	 }
    }
    for ( unsigned int j=0; j<${input_dim}; j++) {
    	if (result[j] < 0) {
	   printf("negative sum of square of deviations\n");
	}
	if (size <= 0) {
	   printf("non positive size\n");
	}
    	result[j] /= size;
	result[j] = sqrt(result[j]);
    }
<%def name="stdev()">
    float mu_carry = 0.0;
    float mu_curr = 0.0;
    float var_carry = 0.0;
    float var_curr = 0.0;
    int k = 0;
    for( int i = 0; i<size/LINE_SIZE; i++ ){
	for( int j = k; j<k+LINE_SIZE; j++ ){
	    mu_curr += data[j];
	}
	mu_curr /= LINE_SIZE;
	for( int j= k; j<k+LINE_SIZE; j++ ){
	    float residual = data[j] - mu_curr;
	    var_curr += residual*residual;
	}
	var_carry = update_var( mu_carry, mu_curr, var_carry, var_curr, k, LINE_SIZE );
	mu_carry = update_mean( mu_carry, mu_curr, k, LINE_SIZE );
	mu_curr =  var_curr = 0.0;
	k += LINE_SIZE;
    }
    // The leftovers
    for( int j=k; j<size; j++ ){
	mu_curr += ${datum};
    }
    mu_curr /= size - k;
    for( int j = k; j<size; j++ ){
	float residual = ${datum};
	var_curr += residual*residual;
    }
    mu_curr /= size - k;
    var_carry = update_var( mu_carry, mu_curr, var_carry, var_curr, k, size-k);
    float stdev = sqrt( var_carry );
</%def>

##produce the classifier from the requested function
<%def name="make_classifier(func_name, input_dim, output_dim)">
    <%
        body = self.template.get_def(func_name).render(False, input_dim, output_dim)
    %>
void compute_estimate( float * const data, unsigned int * const indicies, unsigned int size, float* const result ){
      	 ${body}
}
</%def>

##produce the bootstrap reducer from the requested function
<%def name="make_reduce_bootstraps( func_name, input_dim, output_dim )">
    <%
        body = self.template.get_def(func_name).render(True, input_dim, output_dim)
    %>
void reduce_bootstraps( float * const data, unsigned int  size, float* const result ){
     	 ${body}
}
</%def>

##produce the subsample reducer from the requested function
<%def name="make_average( func_name, input_dim, output_dim )">
    <%
        body = self.template.get_def(func_name).render(True, input_dim, output_dim)
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