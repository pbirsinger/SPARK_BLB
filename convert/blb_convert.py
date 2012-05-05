
import ast
import asp.codegen.cpp_ast as cpp_ast
import asp.codegen.ast_tools as ast_tools 
from blb_convert_data_model import DataModel, DataModelView, ReturnModel, ConstDataModel
import blb_convert_functions as functions

class BLBConverter( ast_tools.ConvertAST ):
    def __init__( self, data_model, input_size, weighted=False ):
	self.weighted = weighted
	self.desired_funcs = []
	# mapping from names to model elements
	self.data_model = {}
	# raw in-order list representation
	self._data_model = data_model
	self.loopvar = [ '_blb_i' ]
	self.registers = {}
	self.retModel = None
	self.input_size = input_size
	self.aggregates = []
	self.arg_model = []

    def aggregate_on( self, ag_string ):
	self.aggregates.append( ag_string )
	for model in self.data_model.itervalues():
	    if model.should_declare():
   	        model.weight_with( ag_string )

    def aggregate_off( self, ag_string ):
	if ag_string in self.aggregates:
	    self.aggregates.remove( ag_string )
	    for model in self.data_model.itervalues():
		if model.weight_index == ag_string:
		    model.weight_with( None )

    def get_or_create_model( self, node ):
	""" Returns the data model appropriate to the given node.
	    Node is assumed to have been already transformed.
	"""
	if node in self.data_model:
	    return self.data_model[ node ]
	elif type( node ) == cpp_ast.CName:
	    if node.name in self.data_model:
		return self.data_model[ node.name ]
	    elif node.name in self.loopvar:
		return DataModel( int, [1], None, node.name )
	    else:
		raise ValueError( 'Unknown data object: %s' % node.name )
	elif type( node ) == cpp_ast.Subscript:
	    return self.get_or_create_model( node.value ).branch()
	elif isinstance( node, DataModel ):
	    return node
	else:
	    raise TypeError( '%s does not represent a data object' % str( node ) )

    def get_register( self, dtype, dim ):
	""" Returns a pair of name, index representing a register. """
	idx = -1
	if not (dtype, dim) in self.registers:
	    self.registers[(dtype,dim)] = [] 
	reglist = self.registers[(dtype,dim)]
	for i in range( len( reglist ) ):
	    if reglist[i]:
		reglist[i] = False
		idx = i
		break
	if idx < 0:
	    reglist.append( False )
	    idx = len( reglist ) - 1
	regname = register_from_spec( dtype, dim, idx )
	try:
	    return self.get_or_create_model( cpp_ast.CName( regname ) ), idx
	except ValueError:
	    register = DataModel( dtype, [dim], name = regname )
	    self.data_model[ regname ] = register
	    register._declare = True
	    return register, idx

    def get_ret_model( self ):
	return self.retModel

    def release_register( self, dtype, dim, i ):
	self.registers[(dtype, dim)][i] = True
		
    def inject_registers( self, body ):
	for model in self.data_model.itervalues():
	    if model.should_declare():
		#TODO: This is an awkward way to do this.
		body.contents.insert(0, cpp_ast.Call( 'memset', [ model.ref_name(), 0, model.size()*model.scalar_t.csize() ] ) )
		body.contents.insert(0, cpp_ast.Value( model.scalar_t.ctype(), cpp_ast.CName("%s [%d]" % ( model.ref_name(), model.size() ) ) ) )

    def visit_arguments( self, node ):
	args = map( lambda arg: self.visit( arg ), node.args )
	defaults = map( lambda tup: self.visit( tup ), node.defaults )
	annotations = dict()
	for tup in defaults:
	    annotations[tup[0].text] = tuple([elt.text for elt in tup])
	print 'annotations', annotations
	ret = []
	if len( self._data_model ) != len( args ):
	    raise TypeError( 'Expected %d arguments, received %d' % ( len( args ), len( self._data_model ) ) )
	foo = zip( args, self._data_model )
	for arg, model in foo:
	    self.data_model[ arg.name ] = model
	    model.name = arg.name
	    ret.append(cpp_ast.Pointer(cpp_ast.Value(model.scalar_t.ctype(), arg.name)))
	    print "Argument '%s': %s" % ( arg.name, repr( model ) )
	    self.arg_model.append( model )
	    model.should_subsample = not ( model.name in annotations and 'nosubsample' in annotations[model.name] )
	if self.weighted:
	    ret.append(cpp_ast.Pointer(cpp_ast.Value("const unsigned int", cpp_ast.CName( '_blb_weights' ) ) ) )
	return ret

    def visit_Assign( self, node ):
	lhs = self.visit( node.targets[0] )
	rhs = self.visit( node.value )
	#find or create target model
	#pass pointer into opencall, or do what is necessary
	target = None
	try:
	    target = self.get_or_create_model( lhs )
	except ValueError:
	    pass
	value = None
	try:
	    value = self.get_or_create_model( rhs )
	except TypeError:
	    pass
	if value is not None: #RHS is a data object
	    if target is not None:
		if type( target ) != DataModelView:
		    raise ValueError( "'%s' refers to a real buffer and cannot be reassigned to '%s'" % ( str(lhs), str(rhs) ) )
		#TODO we allways assume here that the name must be assigned, which is not always the case.
		name = target.name
		self.data_model[ name ] = DataModelView( value, name )
		return cpp_ast.Assign( cpp_ast.CName( target.ref_name() ), cpp_ast.CName( value.ref_name() ) )
	    elif type( lhs ) == cpp_ast.CName:
		self.data_model[ lhs.name ] = DataModelView( value, lhs.name )
		data_type = value.scalar_t.ctype() if value.is_scalar() else value.scalar_t.ctype() + "*"
		return cpp_ast.Assign( cpp_ast.Value( data_type, lhs.name ), cpp_ast.CName(value.ref_name()) )
	    else:
		raise ValueError( "could not assign to '%s': not a data object or name" %  str( lhs ) )
	elif type( rhs ) == OpenCall:
	    if target is None and type( lhs ) == cpp_ast.CName:
		params = rhs.get_output_params()
		self.data_model[ lhs.name ] = target = DataModel( params['type'], params['shape'], None, lhs.name )
		pre = target.declare( self )
		post =  rhs.write_to( target, self )
		return cpp_ast.UnbracedBlock( [pre, post] )
	    elif target:
		return rhs.write_to( target, self )
	    else:
		raise ValueError( "could not assign to '%s': not a data object or name" % str(lhs) ) 
	elif type( rhs ) == cpp_ast.CNumber:
	    if target is None and type( lhs ) == cpp_ast.CName:
		self.data_model[ lhs.name ] = target = DataModel(  type( rhs.num ) , [1], None, lhs.name )
		return cpp_ast.Initializer( cpp_ast.Value( target.scalar_t.ctype(), lhs.name ), rhs.num )
	    elif target:
		assert target.is_scalar(), "'%s' is not a scalar data object" % target.name 
		assert target.scalar_t.matches( type( rhs.num ) ), \
		     "Type mismatch: '%s' is type '%s', but '%s' is type '%s'" % ( target.name, target.scalar_t, rhs.num, type(rhs.num ) ) 
		return cpp_ast.Initializer( cpp_ast.Value( taget.scalar_t.ctype(), target.name ), rhs.num )
	    else:
		raise ValueError( "could not assign to '%s': not a data object or name" % str(lhs) ) 
	else:
	    raise ValueError( "could not assign from '%s'" % str( rhs ) )
	
    
    def visit_AugAssign( self, node ):
	left = self.visit( node.target )
	right = self.visit( node.value )
	op = self.visit( node.op )
	target = self.get_or_create_model( left )
	return OpenCall( self, left, op, right ).write_to( target, self )

    def visit_BinOp( self, node ):
	left = self.visit( node.left )
	right = self.visit( node.right )
	op = self.visit( node.op )
	if isinstance( left, cpp_ast.CNumber ) and isinstance( right, cpp_ast.CNumber ):
	    num = eval( '%s %s %s' ) % ( left.num, op, right.num )
	    return ConstDataModel( num, type(num) )
	else:
	    return OpenCall( self, left, op, right )
	
    def visit_Pow( self, node ):
	return '**'

    def visit_Call(self, node):
	func = node.func.id
	args = [ self.visit( arg ) for arg in node.args ]
	if functions.is_utility( func ):
	    kargs = dict([ (keyword.arg, self.visit(keyword.value) ) for keyword in node.keywords ])
	    return functions.do_utility_func( func, args, self, kargs )
	elif functions.is_operation( func ):
	    return OpenCall( self,  args[0], func, args[1] if len(args) > 1 else None )
	elif functions.is_productivity( func ):
	    return functions.do_productivity_func( func, args, self )
	else:
	    raise ValueError( "Unsupported function: '%s'" % func )

	
    def visit_DataModel( self, node ):
	return cpp_ast.CName( node.ref_name() )

    def visit_For( self, node ):
	_iters = self.visit( node.iter )
	targets = self.visit( node.target )

	#are we iterating over a data object?
	if type( _iters ) is cpp_ast.CName and _iters.name in self.data_model:
	    _iters = [ _iters ]
	    targets = [ targets ]	
	#co-iterating over a set of datasets!
	if type( _iters ) is list:
	    #set up the loop variable and weight
	    self.loopvar.append( self.loopvar[-1] + 'i' )
	    loopvar = cpp_ast.CName( self.loopvar[-1] )
	    body = []

	    if self.weighted:
        	self.aggregate_on( loopvar )	
	
	    # add the temporary children to the data model
	    for target, _iter in zip( targets, _iters ):
		if not (type(_iter) is cpp_ast.CName and _iter.name in self.data_model):
		    raise ValueError( "Cannot iterate over unknown data object: '%s'" % _iter.name )  
		parent_model = self.data_model[ _iter.name ]
		if type( target ) is not cpp_ast.CName:
		    raise ValueError( "Not a valid iteration target: '%s'" % str(target) )
		if target.name in self.data_model:
		    raise ValueError( "Cannot reuse iteration target: '%s'" % target.name )
	        target_model = self.data_model[ target.name ] = parent_model.branch( name = target.name )
 	        ctype = target_model.scalar_t.ctype()
	        decl = cpp_ast.Value( ctype, target.name ) if target_model.is_scalar() else cpp_ast.Pointer( cpp_ast.Value( ctype, target.name ) ) 
	        init = cpp_ast.Subscript( _iter, loopvar) if target_model.is_scalar() \
		   else cpp_ast.BinOp( _iter, '+', cpp_ast.BinOp( loopvar , '*', parent_model.element_size() ) )
	        body.append( cpp_ast.Assign( decl, init ) )
	    

	    # visit the body of the for
	    body += [ self.visit( x ) for x in node.body ]
  

	    # generate the C for intializer
	    ret = cpp_ast.RawFor( cpp_ast.Assign( cpp_ast.Value( 'int', loopvar ), cpp_ast.CNumber(0) ), \
		   cpp_ast.Compare( loopvar, '<', cpp_ast.CNumber( self.data_model[ _iter.name ].dimensions[0] ) ), \
		   cpp_ast.UnaryOp( '++', loopvar ), cpp_ast.Block( body ) )

	    # remove temporary children from data model and otherwise clean up
	    self.loopvar = self.loopvar[:-1]
	    for target in targets:
	        del self.data_model[ target.name ]
	    if self.weighted:	
	        self.aggregate_off( loopvar )
	    # return
	    return ret

	#or perhaps, a known range.
	elif isinstance( _iters, cpp_ast.CNumber ):
	    return cpp_ast.RawFor( cpp_ast.Assign( cpp_ast.Value( 'int', targets ), cpp_ast.CNumber(0) ), \
		    cpp_ast.Compare( targets, '<', _iter ), cpp_ast.UnaryOp( '++', target ), \
		    cpp_ast.Block([ self.visit( x ) for x in node.body ] ) )

	else:
	    raise ValueError('Loop iterand "%s" is not a numeric expression or known data object' % str( _iters ) )

    def visit_FunctionDef(self, node):
	declarator = cpp_ast.Value( 'void', node.name )
	args = self.visit( node.args )
        body = cpp_ast.Block([self.visit(x) for x in node.body])
	return declarator, args, body

    def visit_Num( self, node ):
	return ConstDataModel( node.n, type(node.n) )

    def visit_Return( self, node ):
	""" A return statement ought to be checked for safety, then
	    translated into a memcpy to output memory.
	"""
	stmt = self.visit( node.value )
	
	if isinstance( stmt, OpenCall ):
	    output = stmt.get_output_params()
	    self.retModel = ReturnModel( output['type'], output['shape'] )
	    pre = stmt.write_to( self.retModel )
	    pre.contents.append( cpp_ast.ReturnStatement( "" ) ) 
	    return pre
	elif isinstance( stmt, cpp_ast.CNumber ):
	    self.retModel = ReturnModel( type( stmt.num ), [1] )
	    return cpp_ast.UnbracedBlock( [ cpp_ast.Statement('*_blb_result = %s;' % stmt.num ), cpp_ast.ReturnStatement( "" ) ] )
	else:
	    try:
		source = self.get_or_create_model( stmt )
		self.retModel = ReturnModel( source.scalar_t, source.dimension() )
		return cpp_ast.UnbracedBlock( [cpp_ast.FunctionCall( 'memcpy', [ '_blb_result', source.ref_name(), source.dimension()*source.scalar_t.csize() ] ), cpp_ast.ReturnStatement( "" ) ] )
	    except ValueError, TypeError:
		raise ValueError( "Invalid return object: %s" % str( stmt ) )

    def visit_Subscript( self, node ):
	index = self.visit( node.slice )
	value = self.visit( node.value )
	if type( node.ctx ) == ast.Load:
	    if type( value ) == cpp_ast.CName:
		target_model = self.get_or_create_model( value )
		if target_model.is_scalar():
		    raise ValueError( 'Invalid target for indexing: %s is scalar' % value.name )
		child_model = target_model.branch( idx = index.name )
                if not child_model.weight_index:
		    child_model.weight_with( index.name )
		ret = cpp_ast.Subscript( value, index ) 
		if child_model.is_scalar() and self.weighted:

		    ret = WeightOp( child_model.weight_index ,  ret )
		    self.data_model[ ret ] = child_model
		    return ret
		else:
		    self.data_model[ ret ] = child_model
		    return ret
	    elif type( value ) == cpp_ast.Subscript:
		target_model = self.get_or_create_model( value )
		if target_model.is_scalar():
		    raise ValueError( 'Invalid target for indexing: %s' % value.name )
		child_model = target_model.branch()
		new_index = cpp_ast.CName( '(((%d)*(%s))+(%s))' % ( target_model.dimensions[0], value.index.generate(), index.generate() ) )		
		ret = cpp_ast.Subscript( value, new_index )
		if child_model.is_scalar() and self.weighted:
		    # This isn't fully general, and in fact only works for two level indexing. 
		    ret = WeightOp( child_model.weight_index if child_model.weight_index else index,  value )
		    self.data_model[ ret ] = child_model
		    return ret
		else:
		    ret = cpp_ast.Subscript( value.value, new_index )
		    self.data_model[ ret ] = child_model
		    return ret
	    else:
		raise ValueError( 'Invalid target for indexing: %s' % str( value ) )
	else:
	    #Nothing fancy here, just make sure everything gets written out right
	    if type( value ) == cpp_ast.CName:
		return cpp_ast.Subscipt( value, index )
	    elif type( value ) == cpp_ast.Subscript:

		if not target_model or target_model.is_scalar():
		    raise ValueError( 'Invalid target for indexing: %s' % value.name )
		child_model = target_model.branch()
		new_index = cpp_ast.CName( '(((%d)*(%s))+(%s))' % ( target_model.dimensions[0], value.index.generate(), index.generate() ) )
		return cpp_ast.Subscript( value.value, new_index )
	    else:
		raise ValueError( 'Invalid target for indexing: %s' % str( value ) )

	

    def visit_Tuple( self, node ):
	return [ self.visit(elt) for elt in node.elts ]

    def render( self, node ):
        declarator, args, body = self.visit(node)
	self.inject_registers( body )
	args.append(cpp_ast.Pointer(cpp_ast.Value( self.retModel.scalar_t.ctype(), '_blb_result' ) ) )  
	model = cpp_ast.FunctionBody( cpp_ast.FunctionDeclaration( declarator, args ), body )
	self.desired_funcs.extend( functions.flush_requests() )
	return str( model )
	

    def mangle( self, func_name ):
	return "%s%s_%d" % ( "weighted_" if self.weighted else "", func_name, self.dim )

    def output_dim( self ):
	return self.retModel.dimension()

def create_data_model( args, sub_len ):
    """ Creates an abstract representation of the scalar type and dimensionality of argument
	data needed by the Nodetransformer to do its job.
    """
    models = []
    for arg in args:
	model = DataModel( arg[1], list(arg[0]) )
	model.set_len( arg[0][0] )
	model.dimensions[0] = sub_len
	assert len( model ) == arg[0][0], "blb_convert: 328: set_len didn't work right..."
	models.append( model )
    return models
  			
	
class WeightOp( cpp_ast.BinOp ):
    def __init__(self, weight_index, target ):
	self.left = cpp_ast.Subscript( cpp_ast.CName( '_blb_weights' ), weight_index )
	self.right = target
	self.op = '*'
	super(cpp_ast.BinOp, self).__init__()

def register_from_spec( dtype, dim, num ):
    return '_blb_%s_%d_reg%d' % ( dtype.ctype(), dim, num )

class OpenCall( ast.AST ):
    literal_types = set([ cpp_ast.CNumber, cpp_ast.CName, cpp_ast.Subscript, WeightOp, ConstDataModel])
    def __init__( self, converter, left, op, right = None ):
	self.left = left
	self.op = op
	self.right = right
	self.conv = converter

    @classmethod
    def is_literal( cls, thing ):
	return type(thing) in cls.literal_types

    def literal_value( self, literal ):
	if not self.is_literal( literal ):
	    raise TypeError( 'Internal Error: cannot get literal value of %s' % str( literal ) )
	if type( literal ) == cpp_ast.CNumber:
	    return literal
	elif type( literal ) == WeightOp:
	    return self.conv.get_or_create_model( literal )
	else:
	    return self.conv.get_or_create_model( literal )  

    def get_output_params( self ):
	""" Returns a dictionary containing the shape and type of the output, using 'shape' and 'type' as keys """
	left = right = None
	if self.is_literal( self.left ):
	    left = self.literal_value( self.left )
	else:
	    left = self.left.get_output_params()
	if self.is_literal( self.right ):
	    right = self.literal_value( self.right )
	elif self.right is None:
	    pass
	else:
	    right = self.right.get_output_params()
	return functions.get_output_params( self.op, [ left, right ] )

    def write_to( self, pointer, op=ast.Store() ):
	res = self.conv.get_or_create_model( pointer )
	if self.right is None:
	    if self.is_literal( self.left ):
		return cpp_ast.UnbracedBlock( [ functions.function_from_op( self.op,  [ self.literal_value(self.left) ], res  ) ]  )
	    elif type( self.left ) == OpenCall:
		pre = self.left.write_to( pointer, self.conv )
		pre.contents.append( functions.function_from_op( self.op, [ res ], res, in_place = True ) )
		return pre
	    else:
		raise TypeError( 'Invalid argument to unary operation: %s' % str( self.left ) )
	elif self.is_literal( self.right ):
	    if self.is_literal( self.left ):
		return cpp_ast.UnbracedBlock( [ functions.function_from_op( self.op, [ self.literal_value( self.left ), self.literal_value( self.right ) ], res ) ] )
	    elif type( self.left ) == OpenCall:
		pre = self.left.write_to( pointer, self.conv )
		pre.contents.append( functions.function_from_op( self.op, [ res, self.literal_value( self.right )], res, in_place = True ) )
		return pre
	    else:
		raise TypeError( 'Invalid left argument to binary operation: %s' % str( self.left ) ) 
	elif  type( self.right ) == OpenCall:
	    if self.is_literal( self.left ):
		register, regnum = self.conv.get_register( res.scalar_t, res.dimension() )
		post = self.right.write_to( register, self.conv )
		post.contents.append( functions.function_from_op( self.op, [ self.literal_value( self.left ), register ], res, in_place = True ) )
		self.conv.release_register( res.scalar_t, res.dimension(), regnum )
		return post
	    elif type( self.left ) == OpenCall:
		register, regnum = self.conv.get_register( res.scalar_t, res.dimension() )
		
		pre = self.left.write_to( pointer, self.conv )
		post = self.right.write_to( register, self.conv )
		pre.contents.extend( post.content )
		pre.contents.append( functions.function_from_op( self.op, [ res, register ], res, in_place = True ) )
		self.conv.release_register( res.scalar_t, res.dimension(), regnum )
		return pre
	    else:
		raise TypeError( 'Invalid left argument to binary operation: %s' % str( self.left ) )
	else:
	    raise TypeError( 'Invalid right argument to binary operation: %s' % str( self.right ) ) 



"""CODE DUNGEON"""
''' 
This is where old code is left to rot until it is determined to be
truly unnecessary
''' 

"""

old visit_Call

	if node.func.id == 'linreg':
	    self.out_dim = self.dim - 1
	elif node.func.id == 'mean_norm':
	    self.out_dim = 1
	params=[ 'data' ]
	if self.weighted:
	    params.append('weights')
	params.extend( [ 'size', 'result' ] )  
	func_name = self.mangle( node.func.id )
	self.desired_funcs.append( ( node.func.id, func_name, self.weighted, self.dim, self.dim ) )
	return cpp_ast.FunctionCall(cpp_ast.CName(func_name), params )
"""