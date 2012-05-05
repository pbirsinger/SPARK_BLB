import ast, numpy
import asp.codegen.cpp_ast as cpp_ast

class DataModel( ast.AST ):
    def __init__( self, scalar_type, dimensions, parent = None, name = "" ):
	self.scalar_t = scalar_type if isinstance( scalar_type, RobustType ) else RobustType( scalar_type )
	self.dimensions = dimensions
	self.parent = parent
	self.weight_index = None
	self.name = name
	self.index_names = []
	self.index_multipliers = []
	self._declare = False
	self.length = dimensions[0]
	self.should_subsample = False
	self._fields = [ 'scalar_t', 'dimensions', 'parent', 'weight_index', 'name', 'index_names', 'index_multipliers', '_declare', 'length', 'should_subsample' ]

    def __len__( self ):
	return self.length

    def __repr__( self ):
	return "<DataModel %s, %s, %s, %s >" % ( self.scalar_t.npy_type(), self.dimensions, self.parent, self.name )
 
    def __str__( self ):
	return self.ref_name()
 
    def branch( self, name = "", idx = None ):
	if self.is_scalar():
	    # TODO: figure out how to make this more informative
	    raise TypeError( 'Scalar datum cannot be decomposed' )
	dims = self.dimensions[1:] if len( self.dimensions ) > 1 else [1]
	_name = name if name else self.name
	child = DataModel( self.scalar_t, dims, self, _name )
	child.weight_with( self.weight_index )
	child.index_names = [ idxn for idxn in self.index_names ]
	if idx is not None:
            child.index_names.append( idx )
	    child.index_multipliers.append( self.element_size() ) 
	return child

    def is_scalar( self ):
	return (len( self.dimensions ) == 1) and (self.dimensions[0] == 1) 
	     
    def print_name( self ):
	""" For debugging purposes """
	print self.name
	print self.index_names
	print self.index_multipliers

    def weight_with( self, index ):
	self.weight_index = index  

    def declare( self, converter ):
	""" Returns a c++ statement declaring this object. """
	if self.is_scalar():
	    return cpp_ast.Assign( cpp_ast.Value( self.scalar_t.ctype(), self.ref_name() ), cpp_ast.CNumber(0) )
	else:
	    self._declare = True
	    return cpp_ast.Expression()

    def dimension( self ):
	return self.dimensions[0]

    def element_size( self ):
	if len( self.dimensions ) == 1:
	    return 1
	else:
	    return reduce( int.__mul__, self.dimensions[1:]  )

    def is_weighted( self ):
	return self.weight_index is not None

    def get_weight( self ):
	return self.weight_index[0] if hasattr(self.weight_index, '__iter__') else self.weight_index

    def get_aggregates( self ):
	return self.weight_index

    def ref_name( self ):
	if len( self.index_names ) == 0:
	    return self.name
	if self.is_scalar():
	    parts = map( lambda idx, idn: "(%d*%s)" % ( idx, idn ), self.index_multipliers, self.index_names )
	    return self.name + "[" + '+'.join( parts ) + "]"
	else:
	    parts = [ self.name ]
	    parts.extend( map( lambda idx, idn: "(%d*%s)" % ( idx, idn ), self.index_multipliers, self.index_names ) )
	    return "+".join( parts )

    def set_len( self, n ):
	self.length = n

    def should_declare( self ):
	return self._declare

    def size( self ):
	return reduce( int.__mul__, self.dimensions )
    def clone( self ):
	_clone = DataModel( self.scalar_t, [1], self.parent, self.name )
	for field in self._fields:
	    setattr( _clone, field, getattr( self, field ) )
	_clone.dimensions = self.dimensions[:]
	return _clone

class DataModelView( DataModel ):
    fields = ['base', 'name']
    def __init__( self, base, name ):
	object.__setattr__( self, 'base', base )
	object.__setattr__( self, 'name', name )
	
    def __setattr__( self, name, value ):
	setattr( self.base, name, value )

    def __getattribute__( self, name ):
	try:
	    return object.__getattribute__( self, name )
	except AttributeError:
	    return getattr( self.base, name )
	
    def branch( self, name = "", idx = None ):
	self.base.branch( name, idx )

    def is_scalar( self ):
	return self.base.is_scalar()

    def weight_with( self, index ):
	self.base.weight_with( index )

    def dimension( self ):
	return self.base.dimension()

    def is_weighted( self ):
	return self.base.is_weighted()

    def ref_name( self ):
	return self.name

    def should_declare( self ):
	return ( self.base.name == "" ) and self.base.should_declare()

class ReturnModel( DataModel ):
    def __init__( self, scalar_type, dimension ):
	DataModel.__init__( self, scalar_type, dimension if type(dimension) is list else [dimension], None, '_blb_result' )

    def should_declare( self ):
	return False

    def ref_name( self ):
	return '*_blb_result' if self.is_scalar() else '_blb_result'

class ConstDataModel( DataModel, cpp_ast.CNumber ):
    def __init__( self, num, dtype ):
	DataModel.__init__( self, dtype, [1] )
	cpp_ast.CNumber.__init__( self,  num )

    def ref_name( self ):
	return str( self.num )

# ( ctype, numpy type, python type, c size in bytes )
__ALIAS_CATALOGUE = {}

def register_robust_type( tp_tuple ):
    for tp in tp_tuple:
	__ALIAS_CATALOGUE[ tp ] = tp_tuple

def get_type_aliases( tp ):
    return __ALIAS_CATALOGUE[ tp ]

register_robust_type( ( 'int', numpy.dtype('int32'), int, 4 ) )
register_robust_type( ( 'float', numpy.dtype('float32'), float, 4 ) )
register_robust_type( ( 'double', numpy.dtype('float64'), float, 8 ) )
register_robust_type( ( 'long', numpy.dtype('int64'), long, 8 ) )

class RobustType( object ):
    def __init__( self, tp ):
	self.aliases = get_type_aliases( tp )
	
    def ctype( self ):
	return self.aliases[0]

    def npy_type( self ):
	return self.aliases[1]

    def py_type( self ):
	return self.aliases[2]

    def matches( self, tp ):
	return ( tp in self.aliases )

    def csize( self ):
	return self.aliases[3]

    def __eq__( self, other ):
        return self.aliases == other.aliases

    def __hash__( self ):
	return hash( self.aliases[0] ) 

    def zero( self ):
	return self.tp_aliases[2](0)