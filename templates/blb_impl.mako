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
 bootstrap_reducer: the name of the function with which to reduce bootstraps
 subsample_reducer: the name of the function with which to reduce subsamples
 classifier: the name of the function to calculate on the bootstraps
 
 attributes: a dictionary of optional flags that are passed into the classifier constructors
  with_cilk: whether or not the rendered functions should make use of cilk
 
Reducers defined by this file follow the interfaces

float compute_estimate( float* data, unsigned int* indicies, unsigned int size );

float reduce_bootstraps( float* data, unsigned int size );

float average( float * data, unsigned int size );
</%doc>

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
<%def name="make_classifier(func_name)">
    <%
        body = self.template.get_def("weighted_" + func_name).render()
    %>
float compute_estimate( float * data, unsigned int * weights, unsigned int size ){
      	 ${body}
	 return ${func_name};
}
</%def>

##produce the bootstrap reducer from the requested function
<%def name="make_reduce_bootstraps( func_name )">
    <%
        body = self.template.get_def(func_name).render()
    %>
float reduce_bootstraps( float * data, unsigned int size ){
     	 ${body}
	 return ${func_name};
}
</%def>

##produce the subsample reducer from the requested function
<%def name="make_average( func_name )">
    <%
        body = self.template.get_def(func_name).render()
    %>
float average( float * data, unsigned int size ){
   	 ${ body }
	 return ${func_name};
}
</%def>


%if classifier is not UNDEFINED:
float compute_estimate( float* data, unsigned int* weights,  unsigned int size ){
    ${classifier}
}
%elif use_classifier is not UNDEFINED:
    ${make_classifier(use_classifier)}
%endif

%if bootstrap_reducer is not UNDEFINED:
float reduce_bootstraps( float* data, unsigned int size ){
    ${bootstrap_reducer}
%elif use_bootstrap_reducer is not UNDEFINED:
    ${make_reduce_bootstraps(use_bootstrap_reducer)}
%endif

%if subsample_reducer is not UNDEFINED:
float average( float* data, unsigned int size ){
    ${subsample_reducer}
}
%elif use_subsample_reducer is not UNDEFINED:
    ${make_average(use_subsample_reducer)}
%endif