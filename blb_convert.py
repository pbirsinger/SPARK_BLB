
import ast
import asp.codegen.cpp_ast as cpp_ast
import asp.codegen.ast_tools as ast_tools 


class BLBConverter( ast_tools.ConvertAST ):
    def __init__( self, dimension, weighted ):
	self.dim = dimension
	self.weighted = weighted
	self.desired_funcs = []
	self.out_dim = self.dim

    def visit_Call(self, node):
	print node.func.id
	if node.func.id == 'linreg':
	    self.out_dim = self.dim - 1
	elif node.func.id == 'mean_norm':
	    self.out_dim = 1
 	print self.output_dim()
	params=[ 'data' ]
	if self.weighted:
	    params.append('weights')
	params.extend( [ 'size', 'result' ] )  
	func_name = self.mangle( node.func.id )
	self.desired_funcs.append( ( node.func.id, func_name, self.weighted, self.dim, self.dim ) )
	return cpp_ast.FunctionCall(cpp_ast.CName(func_name), params )

    def render( self, node ):
        model = self.visit(node)
	lines = str( model ).split('\n')
	#discard the wrongly-constructed function header and ending brace
	return "\n".join( lines[2:-1] )
	

    def mangle( self, func_name ):
	return "%s%s_%d" % ( "weighted_" if self.weighted else "", func_name, self.dim )

    def output_dim( self ):
	return self.out_dim