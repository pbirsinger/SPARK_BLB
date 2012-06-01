import blb
import numpy

class SVMVerifierBLB( blb.BLB ):
    def compute_estimate( emails, tags, models = ('models', 'nosubsample') ):
	errors = 0.0
	for email, tag in emails, tags:
	    choice = 0
	    max_match = -1
	    for model in models:
		match = dot( model, email )
		if match > max_match:
		    choice = index() + 1
		    max_match = match
	    if choice != tag:
		errors += 1
	return errors / len( emails )

    def reduce_bootstraps( bootstraps ):
	mean = 0.0
	for bootstrap in bootstraps:
	    mean += bootstrap
	return mean / len(bootstraps)

    def average( subsamples ):
	mean = 0.0
	for subsample in subsamples:
	    mean += subsample
	return mean / len( subsamples )


if __name__ == '__main__':
    X = numpy.zeros(( 10000, 1000))
    Y = numpy.zeros(1000)
    models = numpy.zeros((5000, 1000))
    svm = SVMVerifierBLB( with_openMP=True )
    svm.compile_for(svm.fingerprint( [ X, Y, models ] ) )
