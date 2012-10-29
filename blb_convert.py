
import ast
import asp.codegen.cpp_ast as cpp_ast
import asp.codegen.ast_tools as ast_tools 
from blb_convert_data_model import DataModel, DataModelView, ReturnModel
import blb_convert_functions as functions

class BLBConverter( ast_tools.ConvertAST ):
    def __init__( self, data_model, input_size, weighted=False ):
        self.weighted = weighted
        self.desired_funcs = []
        # mapping from names to model elements
        self.data_model = {}
        # raw in-order list representation
        self._data_model = data_model
        self.loopvar = '_blb_i'
        self.registers = {}
        self.retModel = None
        self.input_size = input_size
    
    def get_or_create_model( self, node ):
        """ Returns the data model appropriate to the given node.
            Node is assumed to have been already transformed.
        """
        if node in self.data_model:
            return self.data_model[ node ]
        elif type( node ) == cpp_ast.CName:
            try:
                return self.data_model[ node.name ]
            except KeyError:
                raise ValueError( 'Unknown data object: %s' % node.name )
        elif type( node ) == cpp_ast.Subscript:
            return self.get_or_create_model( node.value ).branch()
        
        elif isinstance( node, DataModel ):
            return node
        else:
            raise TypeError( '%s does not represent a data object' % str( node ) )

    def get_register( self, dtype, dim ):
        idx = -1
        reglist = self.registers.getdefault( (dtype, dim), [] )
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
        except KeyError:
            register = DataModel( dtype, dim, name = regname )
            self.data_model[ regname ] = register
            return register, idx

    def get_ret_model( self ):
        return self.retModel.clone()

    def release_register( self, dtype, dim, i ):
        self.registers[(dtype, dim)][i] = True
        
    def visit_arguments( self, node ):
        args = map( lambda arg: self.visit( arg ), node.args )
        ret = []
        if len( self._data_model ) != len( args ):
            raise TypeError( 'Expected %d arguments, received %d' % ( len( args ), len( self._data_model ) ))
        for arg, model in zip( args, self._data_model ):
            self.data_model[ arg.value ] = model
            model.name = arg.value
            ret.append(cpp_ast.Pointer(cpp_ast.Value(model.ctype(), arg)))
            mode.dimensions[0] = min( model.dimension(), self.input_size ) 
        if self.weighted:
            ret.append(cpp_ast.Pointer(Value("const unsigned int", cpp_as.CName( '_blb_weights' ) ) ) )
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
            value = get_or_create_model( rhs )
        except ValueError:
            pass
        except TypeError:
            pass
        if value is not None: #RHS is a data object
            if target is not None:
                if type( target ) != DataModelView:
                    raise ValueError( "'%s' refers to a real buffer and cannot be reassigned to '%s'" % ( str(lhs), str(rhs) ) )
                else:
                    #TODO we allways assume here that the name must be assigned, which is not always the case.
                    name = target.name
                    self.data_model[ name ] = DataModelView( value, name )
                    return cpp_ast.Initializer( cpp_ast.Value( '%s*' % value.scalar_t.ctype(), name ), value.name )
            elif type( lhs ) == cpp_ast.CName:
                self.data_model[ lhs.value ] = DataModelView( value, lhs.value )
                return cpp_ast.Initializer( cpp_ast.Value( '%s*' % value.scalar_t.ctype(), lhs.value ), value.name )
            else:
                raise ValueError( "could not assign to '%s': not a data object or name" %  str( lhs ) )
        elif type( rhs ) == OpenCall:
            if target is None and type( lhs ) == cpp_ast.CName:
                params = rhs.get_output_params()
                self.data_model[ lhs.value ] = target = DataModel( params[1], params[0], None, lhs.value )
                return rhs.write_to( target, self )
            elif target:
                return rhs.write_to( target, self )
            else:
                raise ValueError( "could not assign to '%s': not a data object or name" % str(lhs) ) 
        elif type( rhs ) == cpp_ast.CNumer:
            if target is None and type( lhs ) == cpp_ast.CName:
                self.data_model[ lhs.value ] = target = DataModel(  type( rhs.num ) , [1], None, lhs.value )
                return cpp_ast.Initializer( cpp_ast.Value( target.scalar_t.ctype(), lhs.value ), rhs.num )
            elif target:
                assert( target.is_scalar(), "'%s' is not a scalar data object" % target.name )
                assert( target.scalar_t.matches( type( rhs.num ) ), "Type mismatch: '%s' is type '%s', but '%s' is type '%s'" % ( target.name, target.scalar_t, rhs.num, type(rhs.num ) ) )
                return cpp_ast.Initializer( cpp_ast.Value( taget.scalar_t.ctype(), target.name ), rhs.num )
            else:
                raise ValueError( "could not assign to '%s': not a data object or name" % str(lhs) ) 
        else:
            raise ValueError( "could not assign from '%s'" % str( rhs ) )
    

    def visit_AugAssign( self, node ):
        left = self.visit( node.left )
        right = self.visit( node.right )
        op = self.visit( node.op )
        target = self.get_or_create_model( left )
        return OpenCall( left, op, right ).write_to( target, self )

    def visit_BinOp( self, node ):
        left = self.visit( node.left )
        right = self.visit( node.right )
        op = self.visit( node.op )
        if type( left ) == cpp_ast.CNumber and type( right ) == cpp_ast.CNumber:
            num = eval( '%s %s %s' ) % ( left.num, op, right.num )
            return cpp_ast.CNumber( num )
        else:
            return OpenCall( left, op, right )
    
    def visit_Call(self, node):
        func = node.func.id
        args = map( self.visit, node.args )
        if functions.is_utility( func ):
            return functions.do_utility_func( func, args, self )
        elif functions.is_operation( func ):
            return OpenCall( args[0], func, args[1] if len(args) > 1 else None )
        elif functions.is_productivity( func ):
            return functions.do_productivity_func( func, args, self )
        else:
            raise ValueError( "Unsupported function: '%s'" % func )
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

    def visit_For( self, node ):
        iter = self.visit( node.iter )
        target = self.visit( node.target )
        if type( target ) != cpp_ast.CName:
            raise ValueError( '"%s" is not a valid for loop target' % str(target) )
    
        #are we iterating over a data object?
        if type( iter ) is  cpp_ast.CName and iter.name in self.data_model:
            self.loopvar += 'i'
            loopvar = cpp_ast.CName( self.loopvar )
            # add the temporary child to the data model
            target_model = data_model[ target.name ] = data_model[ iter.name ].branch( name == target.name )
            target_model.weight_with( loopvar.name )
            decl = cpp_ast.Value('scalar_t', target.name ) if target_model.is_scalar()    else cpp_ast.Pointer( cpp_ast.Value('scalar_t', target.name ) ) 
            init = cpp_ast.Subscript( iter, loopvar) if target_model.is_scalar() \
               else cpp_ast.BinOp( iter, '+', cpp_ast.BinOp( loopvar , '*', target_model.dimensions[0] ) )
            body = [ cpp_ast.Assign( decl, init ) ]
            
            # visit the body of the for
            body += [ self.visit( x ) for x in node.body ]  
            # generate the C for intializer
            ret = cpp_ast.RawFor( cpp_ast.Assign( Value( 'int', loopvar ), CNumber(0) ), \
               cpp_ast.Compare( loopvar, '<', cpp_ast.CNumber( data_model[ iter.name ].dimensions[0] ) ), \
               cpp_ast.UnaryOp( '++', loopvar ), body )
    
            # remove temporary child from data model
            self.loopvar = self.loopvar[:-1]
            del data_model[ target.name ]
            # return
            return ret
        #or perhaps, a known range.
        elif type( iter ) is cpp_ast.CNumber:
            return cpp_ast.RawFor( cpp_ast.Assign( Value( 'int', target ), CNumber(0) ), \
                cpp_ast.Compare( target, '<', iter ), cpp_ast.UnaryOp( '++', target ), \
                [ self.visit( x ) for x in node.body ] )
        else:
            raise ValueError('Loop iterand "%s" is not a numeric expression or known data object' % str( node.iter ) )

    def visit_FunctionDef(self, node):
        declarator = cpp_ast.Value( 'void', node.name )
        args = self.visit( node.args )
        body = cpp_ast.Block([self.visit(x) for x in node.body])
        return declarator, args, body

    def visit_Return( self, node ):
        """ A return statement ought to be checked for safety, then
            translated into a memcpy to output memory.
        """
        stmt = self.visit( node.value )
    
        if insinstance( stmt, OpenCall ):
            tp, dim = stmt.get_output_params()
            self.retModel = ReturnModel( tp, dim )
            pre = stmt.write_to( retModel )
            pre.body.append( cpp_ast.ReturnStatement( "" ) ) 
            return pre
        elif isinstance( stmt, cpp_ast.CNumber ):
            self.retModel = ReturnModel( type( stmt.num ), [1] )
            return cpp_ast.UnbracedBlock( [ cpp_ast.Statement('*_blb_result = %s;' % stmt.num ), cpp_ast.ReturnStatement( "" ) ] )
        else:
            try:
                source = self.get_or_create_model( stmt )
                self.retModel = ReturnModel( source.scalar_t, source.dimension() )
                return cpp_ast.UnbracedBlock( [cpp_ast.FunctionCall( 'memcpy', [ '_blb_result', source.ref_name(), source.dimension()*source.scalar_t.csize() ] ), cpp_ast.ReturnStatement( "" ) ] )
            except  ValueError, TypeError:
                raise ValueError( "Invalid return object: %s" % str( stmt ) )

    def visit_Subscript( self, node ):
        index = self.visit( node.slice )
        value = self.visit( node.value )
        if type( node.ctx ) == ast.Load:
            if type( value ) == cpp_ast.CName:
                target_model = self.get_or_create_model( value )
                if target_model.is_scalar():
                    raise ValueError( 'Invalid target for indexing: %s is scalar' % value.name )
                child_model = target_model.branch()
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

    

    def render( self, node ):
        declarator, args, body = self.visit(node)
        self.inject_registers( body )
        args.append(cpp_ast.Pointer(cpp_ast.Value( self.retModel.ctype(), '_blb_result' ) ) )  
        model = cpp_ast.FunctionBody( cpp_ast.FunctionDeclaration( declarator, args ), body )
        return str( model )
    

    def mangle( self, func_name ):
        return "%s%s_%d" % ( "weighted_" if self.weighted else "", func_name, self.dim )

    def output_dim( self ):
        return self.retModel.dimension()

def create_data_model( *args ):
    """ Creates an abstract representation of the scalar type and dimensionality of argument
    data needed by the Nodetransformer to do its job."""
    return [ DataModel( arg[1], arg[0] ) for arg in args ]
              
    
class WeightOp( cpp_ast.BinOp ):
    def __init__(self, weight_index, target ):
        super(cpp_ast.BinOp, self).__init__( cpp_ast.Subscript( '_my_weights', weight_index ), '*', target )

class MultiSubscript( cpp_ast.Expression ):
    def __init__(self, value, indicies, offsets):
        self.value = value
        self.indicies = indicies
        self.offsets = offsets
        self._fields = ['value', 'indicies', 'offsets']

    def generate(self, with_semicolon=False):
        components =  self.indicies[0].name  + [ '(%d*%s)' % (offset, index.name) for (offset, index) in zip( self.offsets, self.indicies[1:] ) ]
        yield "%s[%s]" % (self.value, '+'.join( components ) )

    def from_Subscript( node ):
        return MultiSubscript( node.value, [ node.index ], [] )

def register_from_spec( dtype, dim, num ):
    return '_blb_%s_%d_reg%d' % ( dtype, dim, num )

class OpenCall( ast.Node ):
    literal_types = set([ cpp_ast.CNumber, cpp_ast.CName, cpp_ast.Subscript, WeightOp])
    def __init__( self, left, op, right = None ):
        self.left = left
        self.op = op
        self.right = right

    @classmethod
    def is_literal( cls, thing ):
        return type(thing) in cls.literal_types
    @classmethod
    def literal_value( cls, literal, converter ):
        if not is_literal( literal ):
            raise TypeError( 'Internal Error: cannot get literal value of %s' % str( literal ) )
        if type( literal ) == cpp_ast.CNumber:
            return thing
        elif type( literal ) == WeightOp:
            # Think I'm going to get rid of WeightOp
            pass
        else:
            return converter.get_or_create_data_model( literal )  

    def output_params( self ):
        pass
        """ Returns a tuple of the length and scalar type of this call's output"
        return functions.get_output_dim( self.op, [ self.left, self.right ] )
        pass

    def write_to( self, pointer, converter, op='overwrite' ):
    res = converter.get_or_create_model( pointer )
    if self.right is None:
        if self.is_literal( self.left ):
        return cpp_ast.UnbracedBlock( [ function_from_op( self.op,  [ self.literal_value(self.left) ], res  ] )
        elif type( self.left ) == OpenCall:
        pre = self.left.write_to( pointer, converter )
        pre.contents.append( function_from_op( self.op, [ res ], res, in_place = True ) )
        return pre
        else:
        raise TypeError( 'Invalid argument to unary operation: %s' % str( self.left )
    elif self.is_literal( self.right ):
        if is_literal( self.left ):
        return cpp_ast.UnbracedBlock( [ function_from_op( self.op, [ self.literal_valie( self.left ), self.literal_value( self.right ) ], res ) ] )
        elif type( self.left ) == OpenCall:
        pre = self.left.write_to( pointer, converter )
        pre.contents.append( function_from_op( self.op, [ res, self.literal_value( self.right )], res, in_place = True ) )
        return pre
        else:
        raise TypeError( 'Invalid left argument to binary operation: %s' % str( self.left ) ) 
    elif  type( self.right ) == OpenCall:
        if self.is_literal( self.left ):
        post = self.right.write_to( pointer, converter )
        post.insert( 0, function_from_op( self.op, [ self.literal_value( self.left ), res ], res, in_place = True ) )
        elif type( self.left ) == OpenCall:
        register, regnum = converter.get_register( res.scalar_t, res.dimension() )
        
        pre = self.left.write_to( pointer, converter )
        post = self.right.write_to( register, converter )
        pre.content.extend( post.content )
        pre.content.append( function_from_op( self.op, [ res, register ], res, in_place = True ) )
        converter.release_register( res.scalar_t, res.dimension(), regnum )
        return pre
        else:
        raise TypeError( 'Invalid left argument to binary operation: %s' % str( self.left ) )
    else:
        raise TypeError( 'Invalid right argument to binary operation: %s' str( self.right ) ) 
    """