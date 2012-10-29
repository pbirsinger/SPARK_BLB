
import java.util.ArrayList;

def bootstrap(data: ArrayList[Double]): ArrayList[Double] = {
	var bootstrap_vectors = new ArrayList[ArrayList[Double]]()
	for (i <- Range(0, data.size()/DIM)){
		var store_arr = new ArrayList[Double]()
		for (j <- Range(i*DIM, (i+1)*DIM)){
			store_arr.add(data.get(j))
		}			
		bootstrap_vectors.add(store_arr)
	}
	var bootstrap = new ArrayList[Double]()
	var len = bootstrap_vectors.size()
	var gen = new java.util.Random()
	for (i <- Range(0, data.size()/DIM)){
		bootstrap.addAll(bootstrap_vectors.get(gen.nextInt(len)))
	}
	return bootstrap
}


def subsample(data: ArrayList[Double], subsample_len_exp:Double): ArrayList[Double] = {
	var subsample_len = math.pow(data.size(), subsample_len_exp).asInstanceOf[Int]
	var subsample_indicies = new ArrayList[Int]()
	var gen = new java.util.Random()
	for (i <- Range(0, subsample_len)){
		subsample_indicies.add(gen.nextInt(data.size() / DIM))
	}
	var subsample = new ArrayList[Double]()
	//subsample_indicies = subsample_indicies.toArray(new Array[Int])
	for (s <- Range(0, subsample_indicies.size())){
		var index = subsample_indicies.get(s)
		var store_arr = new ArrayList[Double]()
		for (j <- Range((index*DIM).asInstanceOf[Int], ((index+1)*DIM).asInstanceOf[Int])){
			store_arr.add(data.get(j))
		}
		subsample.addAll(store_arr)
	}
	return subsample
}

def run(data: ArrayList[Double], num_subsamples:Int, num_bootstraps:Int, subsample_len_exp:Double):Double={
	var subsample_estimates = new ArrayList[Double]()
	for (i<- Range(0, num_subsamples)){
		var subsamp = subsample(data, subsample_len_exp)
		var bootstrap_estimates = new ArrayList[Double]()
		for (j <- Range(0, num_bootstraps)){
			var btstrap = bootstrap(subsamp)
			var estimate = compute_estimate(btstrap)
			bootstrap_estimates.add(estimate)
		}
		var subsample_est = reduce_bootstraps(bootstrap_estimates)
		subsample_estimates.add(subsample_est)
	}
	return average(subsample_estimates)
}





