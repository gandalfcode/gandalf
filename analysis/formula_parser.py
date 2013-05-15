import numpy as np
import re
from pyparsing import Word, alphas, ParseException, Literal, CaselessLiteral \
, Combine, Optional, nums, Or, Forward, ZeroOrMore, StringEnd, alphanums
from data_fetcher import get_fetcher

# Debugging flag can be set to either "debug_flag=True" or "debug_flag=False"
debug_flag=False

exprStack = []
varStack  = []

def pushFirst( str, loc, toks ):
    exprStack.append( toks[0] )

def assignVar( str, loc, toks ):
    varStack.append( toks[0] )

# define grammar
point = Literal('.')
e = CaselessLiteral('E')
plusorminus = Literal('+') | Literal('-')
number = Word(nums) 
integer = Combine( Optional(plusorminus) + number )
floatnumber = Combine( integer +
                       Optional( point + Optional(number) ) +
                       Optional( e + integer )
                     )

ident = Word(alphas,alphanums + '_') 

plus  = Literal( "+" )
minus = Literal( "-" )
mult  = Literal( "*" )
div   = Literal( "/" )
lpar  = Literal( "(" ).suppress()
rpar  = Literal( ")" ).suppress()
addop  = plus | minus
multop = mult | div
expop = Literal( "^" )
assign = Literal( "=" )

expr = Forward()
atom = ( ( e | floatnumber | integer | ident ).setParseAction(pushFirst) | 
         ( lpar + expr.suppress() + rpar )
       )
        
factor = Forward()
factor << atom + ZeroOrMore( ( expop + factor ).setParseAction( pushFirst ) )
        
term = factor + ZeroOrMore( ( multop + factor ).setParseAction( pushFirst ) )
expr << term + ZeroOrMore( ( addop + term ).setParseAction( pushFirst ) )
bnf = Optional((ident + assign).setParseAction(assignVar)) + expr

pattern =  bnf + StringEnd()

# map operator symbols to corresponding arithmetic operations
opn = { "+" : ( lambda a,b: a + b ),
        "-" : ( lambda a,b: a - b ),
        "*" : ( lambda a,b: a * b ),
        "/" : ( lambda a,b: a / b ),
        "^" : ( lambda a,b: a ** b ) }

# Recursive function that evaluates the stack
def evaluateStack( s, snap ):
  op = s.pop()
  if debug_flag:
      print op
  if op in "+-*/^":
    op2 = evaluateStack( s, snap )
    op1 = evaluateStack( s, snap )
    return opn[op]( op1, op2 )
  elif op == "PI":
    return np.pi
  elif op == "E":
    return np.e
  elif re.search('^[a-zA-Z][a-zA-Z0-9_]*$',op):
      return get_fetcher(op).fetch(snap, "default")[1]
  elif re.search('^[-+]?[0-9]+$',op):
    return int( op )
  elif re.search('[-+]?[0-9]+[.][0-9]*',op):
    return float( op )
  else:
    raise ParserException("Unable to parse string: "  + op)
