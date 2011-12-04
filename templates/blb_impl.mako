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