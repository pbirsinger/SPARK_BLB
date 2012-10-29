
import java.util.ArrayList;
import spark._
import SparkContext._
import javro.scala_arr
import javro.JAvroInter
import org.apache.hadoop.io._
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;

/**
*  converts email from string format into actual Email class
**/
def formatEmail(input: String): Email={
	var vector = input.split(" ")
	var em = new Email()
	em.vec_indices = new Array[Int](vector.length-1)
	em.vec_weights = new Array[Int](vector.length-1)

	var first = true
	var num = 0
	var weight = 0
	var count = 0
	for (elem <- vector){
		if (first){
			em.tag = Integer.parseInt(elem)
			first = false
		}
		else {
			num = Integer.parseInt(elem.substring(0, elem.indexOf(':')))
			weight = Integer.parseInt(elem.substring(elem.indexOf(':')+1, elem.length))
			em.vec_indices(count) = num
			em.vec_weights(count) = weight
			count += 1
		}
	}
	return em
}

/**
*  formatting of the input data can be done here
*  note that the input types may have to be adjusted 
**/
def formatInputItem(input: String): Email={
	return formatEmail(input)
}

/**
 * performs a dot product between a compressed model vector and a
 * compressed email vector stored in the Email class
**/
def custom_dot(model: ArrayList[Float], email: Email): Double ={
    var email_indices = email.get_vec_indices()
    var email_weights = email.get_vec_weights()
    var total =0.0
    var email_index = 0.0
    var email_weight = 0.0
    var model_index = -1.0
    var model_weight =0.0
    var model_index_counter = -2

    for (i <- Range(0, email_indices.length)){
            email_index = email_indices(i)
            while (model_index < email_index && model_index_counter+2 < model.size){
                    model_index_counter += 2
                    model_index = model.get(model_index_counter)
            }
            if (model_index == email_index){
                    email_weight = email_weights(i)
                    model_weight = model.get(model_index_counter+1)
                    total += email_weight * model_weight
            }
    }
    return total
}

def run(email_filename: String, model_filename:String, DIM: Int, 
			num_subsamples:Int, num_bootstraps:Int, subsample_len_exp:Double):Double={
	
	System.setProperty("spark.default.parallelism", "32")// probably want to set to num_nodes * num_cores/node

	// FILE_LOC is set to the file_path of the jar containing the necessary files in /asp/jift/scala_module.py 
	val sc = new SparkContext(System.getenv("MASTER"), "Blb", "/root/spark", List(System.getenv("FILE_LOC")))
	val bnum_bootstraps = sc.broadcast(num_bootstraps)
	val bsubsample_len_exp = sc.broadcast(subsample_len_exp)
	val bnum_subsamples = sc.broadcast(num_subsamples)

	// read data from file to be operated on
	// data needn't necessarily be from a sequence file
    var distData = sc.sequenceFile[Int, String](email_filename)
    var reader =(new JAvroInter("res.avro", "args.avro")).readModel(model_filename)
    var models_arr = List[java.util.ArrayList[Float]]()
    while (reader.hasNext()){
        models_arr = models_arr :+ new ArrayList(reader.next().get(1).asInstanceOf[org.apache.avro.generic.GenericData.Array[Float]].asInstanceOf[java.util.List[Float]])
   }
    var models = sc.broadcast(models_arr)
    
    
    val data_count = distData.count().asInstanceOf[Int]
    val broadcast_data_count = sc.broadcast(data_count)
    val rand_prob = sc.broadcast(math.pow(data_count, subsample_len_exp)/data_count)


    var subsamps = distData.flatMap(item =>{
    	val gen = new java.util.Random()
    	var subsamp_count = 0
    	var outputs = List[(Int, Email)]()
    	var prob =0.0
    	val funcs = new run_outer_data()

    	for (i <- Range(0, bnum_subsamples.value)){
    		prob = gen.nextDouble()
    		if (prob < rand_prob.value){
    			// item = (key, value) so value = item._2
    			outputs ::= (subsamp_count, funcs.formatInputItem(item._2))
    		}
    		subsamp_count += 1
    	}
    	outputs  	
    }).groupByKey().cache()

    // --> (Subsamp Id, subsample) --> (Int, Seq[Email])
    //assuming DIM=1 in map

    var subsamps_out = subsamps.flatMap(subsamp => {
    	//List (subsamp id, subsample) --> (subsample#, Seq(emails))
    	var outputs = List[(Int, Seq[Email])]()
    	for (i <- Range(0, bnum_bootstraps.value)){
    		outputs ::= subsamp
    	}
    	outputs
    }).map(subsamp => {
    	val funcs = new run_outer_data()
    	val subsamp_id = subsamp._1
    	//convert into indexed seq, or even array .toIndexedSeq()
    	var subsamp_vec = subsamp._2.toIndexedSeq

    	val gen = new java.util.Random()
    	var email =""
    	var estimate = 0.0
    	var len = subsamp_vec.size
    	var subsamp_weights = new Array[Int](subsamp_vec.size)

    	for (i <- Range(0, broadcast_data_count.value)){
    		subsamp_weights(gen.nextInt(len)) += 1
    	}
    	
    	//for (temp <- subsamp_vec zip subsamp_weights)
    	for (i <- Range(0, subsamp_weights.length)){
    		subsamp_vec(i).weight = subsamp_weights(i)
    	}
    	
    	var btstrap_data = new BootstrapData()
    	btstrap_data.emails = subsamp_vec.toList
    	btstrap_data.models = models.value
    	(subsamp_id, funcs.compute_estimate(btstrap_data))

    }).groupByKey().map(bootstrap_estimates =>{
    	val funcs = new run_outer_data()
    	funcs.reduce_bootstraps(bootstrap_estimates._2.toList)
    }).collect()
    
    var res = average(subsamps_out)
    System.err.println("res is:" + res)
    return res
}

