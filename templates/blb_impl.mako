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

##Spits out the body of a mean calculation, sans declaration or return statement.
<%def name="mean(iter_direct, dim)">
    <%
        j = 'j' if iter_direct else 'indicies[j]'
	access = 'data[ %s ]' % j
    %>
    float* mean = malloc(sizeof(float) * ${dim});
    for (unsigned int i=0; i<${dim}; i++) {
    	float sub_mean = 0.0;
	for (unsigned int j=0; j<size; j+=${dim}) {
	    sub_mean += ${access};
	}
	sub_mean /= (size / ${dim});
	mean[i] = sub_mean;
    }
</%def>

##Spits out the body of a standard deviation calculation, like unto mean defined above.
<%def name="stdev(iter_direct, dim)">
    ${mean(iter_direct)}
      <%
	access = 'data[i+j]' if iter_direct else 'data[ indicies[i+j] ]'
      %>
    float stdev = 0.0;
    for( unsigned int i=0; i<size; i+=${dim} ){
    	 float squared_norm = 0.0
	 for (unsigned int j=0; j<${dim}; j++) {
	     squared_norm += ${access} - mean[j];
	     squared_norm *= squared_norm;
         stdev += squared_norm;
    }
    stdev = sqrt( stdev / (size / ${dim}) );
</%def>

##produce the classifier from the requested function
<%def name="make_classifier(func_name)">
    <%
        body = self.template.get_def(func_name).render(False)
    %>
float compute_estimate( float * data, unsigned int * indicies, unsigned int size ){
      	 ${body}
	 return ${func_name};
}
</%def>

##produce the bootstrap reducer from the requested function
<%def name="make_reduce_bootstraps( func_name )">
    <%
        body = self.template.get_def(func_name).render(True)
    %>
float reduce_bootstraps( float * data, unsigned int size ){
     	 ${body}
	 return ${func_name};
}
</%def>

##produce the subsample reducer from the requested function
<%def name="make_average( func_name )">
    <%
        body = self.template.get_def(func_name).render(True)
    %>
float average( float * data, unsigned int size ){
   	 ${ body }
	 return ${func_name};
}
</%def>


%if classifier is not UNDEFINED:
float compute_estimate( float* data, unsigned int* indicies,  unsigned int size ){
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