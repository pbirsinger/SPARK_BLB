
import ast
import asp.codegen.cpp_ast as cpp_ast
import asp.codegen.ast_tools as ast_tools 
from blb_convert_data_model import DataModel, RobustType, ConstDataModel


## TODO: add more type/bounds checking!
  
""" Module level data structures """
OP_HANDLERS = {}
REQUESTED_FUNCS = []


""" Module level functions """
def function_from_op( op, op_args, op_res, in_place = False ):
    """ Dispatcher for op to c func conversions.

	INPUTS:
	op: a string representing a vectorizable operation or function
	op_args: a list of arguments to the c func, in invocation order.
	op_res: a DataModel containing the data type and dimensions of the computation,
		as well as the name of the pointer to write to
	in_place: whether or not this should be an inout function in _pointer_
	
	OUTPUT:
	an ast Call node with the appropriate configuration 

	SIDE EFFECTS:
	adds a function request to the appropriate data structure for templatization
    """
    return OP_HANDLERS[ op ]( op_args, op_res, in_place )


def register_op_handler( op_name, op_transformer ):
    OP_HANDLERS[ op_name ] = op_transformer

def flush_requests():
    global REQUESTED_FUNCS
    tmp = set(REQUESTED_FUNCS)	
    REQUESTED_FUNCS = []
    return tmp


def get_output_params( op, args ):
    handler = OP_HANDLERS[ op ]
    return handler.get_output_params( args )

def is_operation( func ):
    return func in OP_HANDLERS

""" Op Handler implementations """

#### vecinit : initializes a vector to a scalar value.

def is_matrix( thing ):
    return isinstance( thing, DataModel ) and len(thing.dimensions) == 2

class BLBCall( cpp_ast.FunctionCall ):
    @classmethod
    def get_output_params( cls, args ):
	if isinstance( args[0], DataModel ):
	    return { 'shape':args[0].dimensions, 'type':args[0].scalar_t }
	elif isinstance( args[1], DataModel ):
	    return { 'shape':args[1].dimensions, 'type':args[1].scalar_t }
	elif type( args[0] ) == dict:
	    return args[0]
	elif type( args[1] ) == dict:
	    return args[1]
	elif type( args[0] ) == cpp_ast.CNumber:
	    return { 'shape':[1], 'type':RobustType( type(args[0].num) ) }

class DotCall( BLBCall ):
    def __init__( self, op_args, op_res, op ):
	self.left = op_args[0]
	self.right = op_args[1]
	self.res = op_res
	if not isinstance( self.left, DataModel ) or not isinstance( self.right, DataModel ):
	    raise TypeError( "Unsupported operands for function 'dot': %s, %s" % ( type(self.left), type(self.right) ) )
	elif self.left.dimension() != self.right.dimension() or not self.res.is_scalar():
	    raise ValueError( "Incompatible dimensions for function 'dot': %d, %d, %d" % \
		 ( self.left.dimension(), self.right.dimension(), self.res.dimension() ) )
	elif op == ast.Add():
	    self.prepare_dota()
	else:
	    self.prepare_dot()

    def prepare_dot( self ):	
	self.fname = '_blb_dot'
	self.params = [ self.left.ref_name(), self.right.ref_name(), '&%s' % self.res.ref_name(), self.left.dimension() ]
	REQUESTED_FUNCS.append( ( '_blb_dot', ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype() ) ) )

    def preapre_dota( self ):
	self.fname = '_blb_dota'
	if self.res.is_weighted():
	    self.params = [ self.left.ref_name(), self.right.ref_name(), '&%s' % self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(), self.left.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), True ) ) )
	else:
    	    self.params = [ self.left.ref_name(), self.right.ref_name(), '&%s' % self.res.ref_name(), self.left.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), False ) ) )

    @classmethod
    def get_output_params( cls, args ):
	return { 'shape':[1], 'type':args[0].scalar_t } 

register_op_handler( 'dot', DotCall )
 
class PlusCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.left = op_args[0]
	self.right = op_args[1]
	self.res = op_res
	# cases to consider
	# is one or both arguments numeric?
	# is one or both arguments wieghted?
	# is the operation performed in-place?
	assert isinstance( self.left, DataModel ) and isinstance( self.right, DataModel ), \
	    "Invalid operands to operation '+': %s, %s" % ( type(self.left), type(self.right) )
	assert self.left.dimension() == self.right.dimension(), \
	    "incompatible data lengths for operation '+': %s,%d %s,%d" % ( self.left.ref_name(), self.left.dimension(), self.right.ref_name(), self.right.dimension()) 
	if self.left.is_scalar() and self.right.is_scalar(): 
	    self.prepare_simple()    
	else:
	    self.prepare_vaxbyo()

    def prepare_simple( self ):
	self.fname = '_blb_PLUS_MACRO'
	a = 1.0
	if self.res.is_weighted():
	    a = '_blb_weights[%s]' % self.res.get_weight()
	self.params = [ self.left.ref_name(), self.right.ref_name(), a, self.res.ref_name() ]
	REQUESTED_FUNCS.append( (self.fname, (self.res.scalar_t.ctype(),) ) )

    def prepare_vaxbyo( self ):
	self.fname = '_blb_vaxbyo'
	if self.res.is_weighted():
	    self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(), self.res.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), True ) ) )
	else:
    	    self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), self.res.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), False ) ) )

register_op_handler( '+', PlusCall )
	    
class MinusCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.left = op_args[0]
	self.right = op_args[1]
	self.res = op_res
	if not (isinstance( self.left, DataModel ) and isinstance( self.right, DataModel )):
	    raise ValueError( 'Unsupported operands for operator "-": %s , %s' % ( type( self.left ), type( self.right ) ) )
	elif not self.left.dimension() == self.right.dimension() and self.left.dimension() == self.res.dimension():
	    raise ValueError( "Operands to operator '-' have incompatible dimension: %s, %s, %s" \
		% (self.left.ref_name(), self.right.ref_name(), self.res.ref_name() ) ) 
	if self.left.is_scalar() and self.right.is_scalar():
	    self.prepare_simple()
	else: 
	    self.prepare_vaxbyo()
	    

    def prepare_simple( self ):
	self.fname = '_blb_MINUS_MACRO'
	a = ('_blb_weights[%s]' % self.res.get_weight()) if self.res.is_weighted() else 1.0 
	self.params = [ self.left.ref_name(), self.right.ref_name(), a, self.res.ref_name() ]
	REQUESTED_FUNCS.append( ( self.fname, ( self.res.scalar_t.ctype(), ) ) )

    def prepare_vaxbyo( self ):
	self.fname = '_blb_vaxmyo'
	if self.res.is_weighted():
	    self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(),  self.res.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), True ) ) )
	else:
	    self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(),  self.res.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), False ) ) )

register_op_handler( '-', MinusCall )

class MultCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.left = op_args[0]
	self.right = op_args[1]
	self.res = op_res
	if not (isinstance( self.left, DataModel ) and isinstance( self.right, DataModel )):
	    raise TypeError( "invalid operands to opeartion '*': %s, %s" % ( type(self.left), type(self.right) ) )
	elif self.left.is_scalar() and self.right.is_scalar():
	    self.prepare_simple()
	elif self.left.is_scalar():
	    self.prepare_vscaleo( self.right, self.left )
	elif self.right.is_scalar():
	    self.prepare_vscaleo( self.left, self.right )
	else:
	    self.prepare_vmulte()

    def prepare_simple( self ):
	self.fname = '_BLB_MULT_MACRO'
	a = ('_blb_weights[%s]' % self.res.get_weight()) if self.res.is_weighted() else 1.0
	self.params = [ self.left.ref_name(), self.right.ref_name(), a, self.res.ref_name() ]
	REQUESTED_FUNCS.append( ( self.fname, (self.res.scalar_t.ctype(),) ) )

    def prepare_vscaleo( self, vector, scalar ):
	self.fname = '_blb_vscaleo'
	a = scalar.ref_name()
	if self.res.is_weighted():
	    a = '%s*_blb_weights[%s]' % ( a, self.res.get_weight() )
	self.params = [ vector.ref_name(), a, self.res.ref_name(), self.res.dimension() ]
	REQUESTED_FUNCS.append( ( self.fname, ( vector.scalar_t.ctype(), None, self.res.scalar_t.ctype() )))

    def prepare_vmulte( self ):
	self.fname = '_blb_vmulte'
	if self.res.is_weighted():
            self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(), self.res.dimension() ]
	    REQUESTED_FUNCS.append ( (self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), True ) ) )
	else:
            self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), self.res.dimension() ]
	    REQUESTED_FUNCS.append ( (self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), False ) ) )

register_op_handler( '*', MultCall )

class OuterProductCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	if len(op_args) != 2:
	    raise TypeError( "'outer_product' takes exactly 2 arguments, %d supplied" % len(op_args) )
	self.u = op_args[0]
	self.v = op_args[1]
	if not ( isinstance( self.u, DataModel ) and isinstance( self.v, DataModel ) ):
	    raise TypeError( "unsupported operands to function 'outer_product': %s, %s" % ( type( self.u ), type( self.v ) ) )	
	if not is_matrix( op_res ):
	    raise TypeError( "invalid output buffer for 'outer_product': %s" % str( op_res ) )
	self.fname = '_blb_outer_product'
	self.params = [ self.u.ref_name() ]
	if self.u.is_weighted():
	    self.params.append( '_blb_weights[%s]' % self.u.weight_index )
	self.params.append( self.v.ref_name() )
	if self.v.is_weighted():
	    self.params.append( '_blb_weights[%s]' % self.v.weight_index )
	self.params.extend( [ op_res.ref_name(), self.v.dimension(), self.u.dimension() ] )
	REQUESTED_FUNCS.append( ('_blb_outer_product', ( self.u.scalar_t.ctype(), self.u.is_weighted(), self.v.scalar_t.ctype(), self.v.is_weighted(), op_res.scalar_t.ctype() ) ) )

    @classmethod
    def get_output_params( cls, args ):
	return { 'shape':[ args[0].dimension(), args[1].dimension() ], 'type':args[0].scalar_t }

register_op_handler( 'outer_product', OuterProductCall )

class DivCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.left = op_args[0]
	self.right = op_args[1]
	self.res = op_res
	if not( isinstance( self.left, DataModel ) and isinstance( self.right, DataModel ) ):
	    raise ValueError( "Unsupported operand types for operator '/': %s, %s" % ( type( self.left ), type( self.right ) ) )
	elif self.left.is_scalar() and self.right.is_scalar():
	    self.prepare_simple()
	elif self.right.is_scalar():
	    self.prepare_vscaleo()
	else:
	    self.preapre_vdive()

    def prepare_simple( self ):
	self.fname = '_blb_DIV_MACRO'
	a = ('_blb_weights[%s]' % self.res.get_weight() )  if self.res.is_weighted() else 1.0
	self.params = [ self.left.ref_name(), self.right.ref_name(), a, self.res.ref_name() ]
	REQUESTED_FUNCS.append( ( self.fname, (self.res.scalar_t.ctype(),) ) )

    def prepare_vscaleo( self ): 
	self.fname = '_blb_vscaleo'
	a = '1.0/%s' % self.right.ref_name()
	if self.res.is_weighted():
	    a = '%s*_blb_weights[%s]' % ( a, self.res.get_weight() )
	self.params = [ self.left.ref_name(), a, self.res.ref_name(), self.res.dimension() ]
	REQUESTED_FUNCS.append( ( self.fname, ( self.left.scalar_t.ctype(), None, self.res.scalar_t.ctype() )))

    def prepare_vdive( self ):
	if self.left.dimension() != self.right.dimension():
	    raise ValueError( "incompatible dimensions for operation '/': %s,%d ,  %s,%d" \
		% ( self.left.ref_name(), self.left.dimension(), self.right.ref_name(), self.right.dimension() ) )

	self.fname = '_blb_vdive'
	if self.res.is_weighted():
            self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(), self.res.dimension() ]
	    REQUESTED_FUNCS.append ( (self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), True ) ) )
	else:
            self.params = [ self.left.ref_name(), self.right.ref_name(), self.res.ref_name(), self.res.dimension() ]
	    REQUESTED_FUNCS.append ( (self.fname, ( self.left.scalar_t.ctype(), self.right.scalar_t.ctype(), self.res.scalar_t.ctype(), False ) ) )

register_op_handler( '/', DivCall )

""" Supported elementwise mathematical functions """

class PowCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.base = op_args[0]
	self.exponent = op_args[1]
	self.res = op_res
	if not( isinstance( self.base, DataModel ) and isinstance( self.exponent, DataModel ) ):
	    raise ValueError( "Invalid operand types to opreation '**': %s, %s" % ( self.base, self.exponent ) )
	elif  self.left.dimension() != self.res.dimension():
	    raise ValueError( "Incompatible argumnet/result dimensions for operation '**': %d, %d" % (self.left.dimension(), self.res.dimension() ) ) 
	if self.left.is_scalar() and self.right.is_scalar():
	    self.prepare_pow() 
	elif self.right.is_scalar():
	    self.prepare_vexpa()
	elif self.left.dimension() == self.right.dimension():
	    self.prepare_vexpe()
	else:
	    raise ValueError( "Incompatible argument dimensions for operation '**': %d, %d" % (self.left.dimension(), self.right.dimension() ) )

    def prepare_pow( self ):
	self.fname = '_blb_pow'
	if self.res.is_weighted():
	    self.params = [ self.base.ref_name(), self.exponent.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight()  ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.base.scalar_t.ctype(), self.exponent.scalar_t.ctype(), self.res.scalar_t.ctype(), True) ) )
	else:
	    self.params = [ self.base.ref_name(), self.exponent.ref_name(), self.res.ref_name() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.base.scalar_t.ctype(), self.exponent.scalar_t.ctype(), self.res.scalar_t.ctype(), False) ) )



    def prepare_vexpa( self ):
	self.fname = '_blb_vexpa'
	if self.res.is_weighted():
	    self.params = [ self.base.ref_name(), self.exponent.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(), self.res.dimension()  ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.base.scalar_t.ctype(), self.exponent.scalar_t.ctype(), self.res.scalar_t.ctype(), True) ) )
	else:
	    self.params = [ self.base.ref_name(), self.exponent.ref_name(), self.res.ref_name(), self.res.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.base.scalar_t.ctype(), self.exponent.scalar_t.ctype(), self.res.scalar_t.ctype(), False) ) )

    def prepare_vexpe( self ):
	self.fname = '_blb_vexpe'
	if self.res.is_weighted():
	    self.params = [ self.base.ref_name(), self.exponent.ref_name(), self.res.ref_name(), '_blb_weights[%s]' % self.res.get_weight(), self.res.dimension()  ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.base.scalar_t.ctype(), self.exponent.scalar_t.ctype(), self.res.scalar_t.ctype(), True) ) )
	else:
	    self.params = [ self.base.ref_name(), self.exponent.ref_name(), self.res.ref_name(), self.res.dimension() ]
	    REQUESTED_FUNCS.append( ( self.fname, ( self.base.scalar_t.ctype(), self.exponent.scalar_t.ctype(), self.res.scalar_t.ctype(), False) ) )


register_op_handler( '**', PowCall )
register_op_handler( 'pow', PowCall )

class SqrtCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.base = op_args[0]
	self.res = op_res
	if not isinstance( self.base, DataModel ):
	    raise ValueError( "Unsupported operand to function 'sqrt': '%s'" % type( self.base ) )
	elif self.base.is_scalar():
	    self.prepare_sqrt()
	else:
	    self.prepare_vsqrte()

    def prepare_sqrt( self ):
	self.fname = '_blb_sqrt'
	self.params = [ self.base.ref_name(), '&%s' % self.res.ref_name()]
	REQUESTED_FUNCS.append( (self.fname, (self.base.scalar_t.ctype(),) ) )

    def prepare_vsqrte( self ):
	self.fname = '_blb_vsqrte'
	self.params = [ self.base.ref_name(), self.res.ref_name(), self.res.dimension() ]
	REQUESTED_FUNCS.append( (self.fname, (self.base.scalar_t.ctype(), self.res.scalar_t.ctype()) )  )
	
register_op_handler( 'sqrt', SqrtCall )



class MVSolveCall( BLBCall ):
    def __init__( self, op_args, op_res, in_place ):
	self.A = op_args[0]
	self.b = op_args[1]
	self.res = op_res
	if not (( isinstance( self.A, DataModel ) and isinstance( self.b, DataModel ) )):
	    raise TypeError( "Unsupported operands to 'mv_solve': %s, %s" % ( type(self.A), type( self.b ) ) )
	elif not isinstance( self.res, DataModel ):
	    raise TypeError( "Unsupported output for 'mv_solve': %s" % type(self.res) ) 
	else:
	    self.prepare_mvsolve( in_place or op_res.ref_name() == self.b.ref_name() )
	
    def prepare_mvsolve( self, in_place ):
	assert len(self.A.dimensions) == 2, "mv_solve: Data object does not represent a matrix: %s" % self.A.ref_name() 
	m, n = self.A.dimensions[0], self.A.dimensions[1]
	assert self.b.dimension() == m, "mv_solve: Invalid dimension for argument vector: should be %d but is %d" % ( m, self.b.dimension() ) 
	assert self.res.dimension() == n, "mv_solve: Invalid dimension for result vector: should be %d but is %d" % ( n, self.res.dimension() ) 
	self.fname = '_blb_mvsolve'
	self.params = [ self.A.ref_name(), self.b.ref_name() ]
	if self.b.is_weighted():
	    raise ValueError( "mv_solve: cannot use weighted vector as argument: %s" % self.b.ref_name() )
	if not in_place:
	    self.params.append( self.res.ref_name() )
        self.params.extend( [ m, n] )
	REQUESTED_FUNCS.append( ( self.fname, ( in_place, ) ) )

    @classmethod
    def get_output_params( cls, args ):
	return { 'shape':[args[0].dimensions[1]], 'type':args[0].scalar_t }
	
register_op_handler( 'mv_solve', MVSolveCall )

""" Supported utility functions """

UTILITY_HANDLERS = {}
def register_utility_handler( func, handler ):
    UTILITY_HANDLERS[ func ] = handler

def is_utility( func ):
    return ( func in UTILITY_HANDLERS )

def do_utility_func( func, args, converter, kargs ):
     return UTILITY_HANDLERS[func]( args, converter, kargs )

def handle_copy( args, converter, kargs ):
    if len(args) < 2:
	raise TypeError( "copy() requires two arguments, %d given" % len(args) )
    target = converter.get_or_create_model( converter.visit(args[0]) )
    source = converter.get_or_create_model( converter.visit(args[1]) )
    assert target.scalar_t == source.scalar_t, "copy(): Type mismatch between %s and %s" % ( target.ref_name(), source.ref_name() ) 
    assert target.dimension() == source.dimentsion(), "copy(): Length mismatch: %s and %s " % (target.ref_name(), source.ref_name() ) 
    return cpp_ast.FunctionCall( 'memcpy', [ target.ref_name(), source.ref_name(), target.dimension()*target.scalar_t.csize() ] )

register_utility_handler( 'copy', handle_copy )

def handle_dim( args, converter, kargs ):
    if len( args ) < 1:
	raise TypeError( "dim() requires at least 1 argument, %d given." % len( args ) )
    target = converter.get_or_create_model( converter.visit(args[0]) )
    n = 0
    try:
	idx = converter.visit(args[1])
	n = idx.num % len(target.dimensions)
    except IndexError:
	pass
    return ConstDataModel( target.dimensions[n], long )

register_utility_handler( 'dim', handle_dim )

def handle_dtype( args, converter, kargs ):
    if len(args) < 1:
	raise TypeError("dtype takes one argument, %d given" % len(args) )
    return cpp_ast.CName( converter.get_or_create_model( args[0] ).scalar_t.ctype() )
    
register_utility_handler( 'dtype', handle_dtype )

def handle_index( args, converter, kargs ):
    if len( args ) == 0:
	return cpp_ast.CName(converter.loopvar[-1])
    elif len( args ) == 1:
	return cpp_ast.CName(converter.loopvar[-1 - args[0]])
    else:
	raise TypeError("Too many arguments to function 'index': %d" % len(args) )
register_utility_handler( 'index', handle_index)
	
def handle_len( args, converter, kargs ):
    if len(args) < 0:
	raise TypeError('len takes one argument, %d given' % len( args ) )
    model = converter.get_or_create_model( args[0] )
    return ConstDataModel( len( model), long )

register_utility_handler( 'len', handle_len )

def handle_range( args, converter, kargs ):
    if len(args) != 1:
	raise TypeError('range takes exactly one argument, %d given' % len(args) )
    return ConstDataModel( args[0], int )
register_utility_handler( 'range', handle_range )
	    
anon_vecs = 0
def handle_vector( args, converter, kargs ):
    global anon_vecs
    
    length = 0
    if isinstance( args[0], cpp_ast.CNumber ):
	length = args[0].num
    else:
	raise TypeError( "vector: invalid length argument: %s" % ars[0] )
    tp = None
    if type( args[1] ) == cpp_ast.CName:
	try:
	    tp = RobustType( args[1].name ) 
	except KeyError:
	    raise TypeError( "vector: invalid data type argument: %s" % args[1] )
    else:
	raise TypeError( "vector: invalid data type argument: %s" % args[1] )
    name = "_blb_anon_vec%d" % anon_vecs
    anon_vecs += 1
    model = DataModel( tp, [length], None, name )
    model._declare = True
    converter.data_model[ name ] = model
    return cpp_ast.CName( name )
    

register_utility_handler( 'vector', handle_vector )



def handle_initialize( args, converter, kargs ):
    if len( args ) < 2:
	raise TypeError( "initialize() requires 2 arguments, %d given." % len(args) )
    target = converter.get_or_create_model( args[0] )
    assert type( args[1] )== cpp_ast.CNumber, "initialize(): Invalid initial value: %s" % str( args[1] ) 
    REQUESTED_FUNCS.append( ( '_blb_vecinit', ( target.scalar_t.ctype(), ) ) )
    return cpp_ast.FunctionCall( '_blb_vecinit', [ target.ref_name(), args[1].num, target.dimension() ] ) 

register_utility_handler( 'initialize', handle_initialize )



def handle_matrix( args, converter, kargs ):
    global anon_vecs
    assert len(args) == 3, "matrix takes 3 arguments, %d supplied" % len(args)
    assert isinstance( args[0], cpp_ast.CNumber )
    assert isinstance( args[1], cpp_ast.CNumber )
    assert isinstance( args[2], cpp_ast.CName )
    name = "_blb_anon_vec%d" % anon_vecs
    anon_vecs += 1
    model = DataModel( args[2].name, [ args[0].num, args[1].num], None, name )
    model.declare = True
    converter.data_model[ name ] = model
    return cpp_ast.CName( name )

register_utility_handler( 'matrix', handle_matrix )

""" Supported productivity functions """
# TODO: implement these
PRODUCTIVITY_HANDLERS =  {}


def is_productivity( func ):
    return ( func in PRODUCTIVITY_HANDLERS )

def do_utiltiy_func( func, args, converter ):
   raise NotImplementedError()

  
"""

Original code for productivity functions.


	 node.func.id == 'linreg':
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