import java.util.ArrayList

class Email extends Serializable{
	var weight = 0
	var tag = 0
	var vec_indices = Array[Int]()
	var vec_weights = Array[Int]()
	
	def get_weight():Int={ return this.weight}
	def get_tag():Int={ return this.tag}
	def get_vec_indices() :Array[Int]={return this.vec_indices}
	def get_vec_weights() :Array[Int]={return this.vec_weights}
}


// This class should be customized to allow for application-specific
// input to compute_estimate
class BootstrapData extends Serializable{
	var emails = List[Email]()
	var models = List[ArrayList[Float]]()
	
	def get_emails(): List[Email]={
		return this.emails
	}

	def get_models():List[ArrayList[Float]]={
		return this.models
	}		
}