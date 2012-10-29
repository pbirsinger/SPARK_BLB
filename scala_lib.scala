
import org.apache.avro.generic.GenericData;
import org.apache.avro.Schema
import util.Random
import java.util.ArrayList;


object scala_lib{

	def getType(item: Object): String ={
		if (item == null){
			return "\"null\"";
		}
		var c = item.getClass();
		var name = c.getName();

		if (name == "java.lang.Double" || name == "double"){
			return "\"double\"";
		}
		else if (name == "java.lang.Integer" || name == "int"){
			return "\"int\"";
		}
		else if (name == "java.lang.String" || name == "org.apache.avro.util.Utf8"){
			return "\"string\"";
		}
		else if (name == "java.lang.Float" || name == "float"){
			return "\"float\"";
		}
		else if (name == "java.lang.Long" || name == "long"){
			return "\"long\"";
		}
		else if (name == "java.lang.Boolean" || name == "boolean"){
			return "\"boolean\"";
		}
		else if (name == "org.apache.avro.generic.GenericData$Array"){
			//return "{\"type\": \"array\", \"items\": {\"type\": \"array\", \"items\": %s }}".format(getType((item)))
			//FIX THIS SO IT DOESN'T NEED TO BE HARDCODED IN 
			return "{\"type\": \"array\", \"items\": \"double\" }"
	
		}
		else {
			System.err.println(name);
			throw new RuntimeException("Unknown Argument Type to Write to Avro File");
		}
	}	

	def rand_choice[T](arr: org.apache.avro.generic.GenericData.Array[T]): T = {
		var len = arr.size()
		var gen = Random
		return arr.get(gen.nextInt(len))
	}

	def rand_sample[T](data: org.apache.avro.generic.GenericData.Array[T], k :Int): 
	    org.apache.avro.generic.GenericData.Array[T] = {	    
	    var len = data.size()                                                
	    var gen = Random
	    var index = 0
	    
	    var item: Object = data.get(0).asInstanceOf[Object]
	    var typ = getType(item)                                           
	    var str = "{\"type\": \"array\", \"items\": %s}".format(typ)
	    var schema: Schema = (new Schema.Parser()).parse(str)
	    var arr = new org.apache.avro.generic.GenericData.Array[T](k, schema);

	    for (s <- Range(0,k)){
	    	index = gen.nextInt(len)
	    	arr.add(data.get(index))
	    }
	    return arr
	  }	
	
	def extend[T](arr1: org.apache.avro.generic.GenericData.Array[T], arr2: org.apache.avro.generic.GenericData.Array[T]):
	    org.apache.avro.generic.GenericData.Array[T] = {	 		
		def len = arr2.size()
		for (s <- Range(0, len)){
			arr1.add(arr2.get(s).asInstanceOf[T])
		}
		return arr1
	}

	def slice[T](arr: org.apache.avro.generic.GenericData.Array[T], start: Int, end:Int):
		org.apache.avro.generic.GenericData.Array[T] = {
		var item: Object = arr.get(0).asInstanceOf[Object]
		var typ = getType(item)
	    var str = "{\"type\": \"array\", \"items\": %s}".format(typ)
	    var schema: Schema = (new Schema.Parser()).parse(str)
	    var arr1 = new org.apache.avro.generic.GenericData.Array[T](end-start, schema)
	    var count = start
	    while (count < end){
	    	arr1.add(arr.get(count))
	    	count +=1
	    }
	    return arr1
	}
	
	def main(args: Array[String]) {
		
	    //var str = makeSchema(12)
	    var str = "{\"type\": \"array\", \"items\": \"double\"}"	    
	    var schema: Schema = (new Schema.Parser()).parse(str)
	    var arr =new org.apache.avro.generic.GenericData.Array[Double](1, schema);
		for (s <- Range(0,10)){
			arr.add((s*3).asInstanceOf[Double])
		}
		println("ARR1 IS:" + arr.toString())
		/**
		var arr2 = rand_sample(arr, 5)
		println("ARR2 IS:"+arr2.toString())
		var arr3 = extend(arr, arr2)
		println("ARR3 IS:" +arr3.toString())
		var arr4 = arr3.iterator()
		println("arr4 is:" + arr4)
		println("next is:" + arr4.next())
		**/
		println(slice(arr, 2, arr.size()))
		println("starting")
		println(rand_choice(arr))
		println(rand_choice(arr))
		println(rand_choice(arr))
		println(rand_choice(arr))
		println("done with random")
		//var a = new org.apache.avro.generic.GenericData.Array[Double](12,schema).add(3.0).add(3.0)
	}
}