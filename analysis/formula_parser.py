#==============================================================================
#  formula_parser.py
#  Contains functions and variables for parsing string-formulas.
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G Rosotti
#
#  GANDALF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  GANDALF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License (http://www.gnu.org/licenses) for more details.
#==============================================================================
import numpy as np
import re
from pyparsing import Word, alphas, ParseException, Literal, \
     CaselessLiteral, Combine, Optional, nums, Or, Forward, ZeroOrMore, \
     StringEnd, alphanums, delimitedList, Group, Suppress
from data_fetcher import UserQuantity, _KnownQuantities

# Debugging flag can be set to either "True" or "False"
debug_flag = False

exprStack = []
varStack  = []


#------------------------------------------------------------------------------
def pushFirst( str, loc, toks ):
    exprStack.append( toks[0] )


#------------------------------------------------------------------------------
def assignVar( str, loc, toks ):
    varStack.append( toks[0] )
    

#------------------------------------------------------------------------------
def arg_number(str, loc, toks):
    
    # Work out the correct number of arguments
    if toks[0] in functions1arg:
        narg = 1
    elif toks[0] in functions2arg:
        narg = 2
    else:
        raise ParseException("We do not recognize your function: " + toks[0])
    
    # Check that the number is correct
    if narg != len(toks)-1:
        raise ParseException("Error: function " + toks[0] + "takes " \
                             + narg + "arguments, you gave "+ len(toks)-1)


# Define grammar
#------------------------------------------------------------------------------
point       = Literal('.')
e           = CaselessLiteral('E')
plusorminus = Literal('+') | Literal('-')
number      = Word(nums) 
integer     = Combine( Optional(plusorminus) + number )
floatnumber = Combine( integer +
                       Optional( point + Optional(number) ) +
                       Optional( e + integer )
                     )
ident       = Word(alphas,alphanums + '_') 
plus        = Literal( "+" )
minus       = Literal( "-" )
mult        = Literal( "*" )
div         = Literal( "/" )
lpar        = Literal( "(" ).suppress()
rpar        = Literal( ")" ).suppress()
addop       = plus | minus
multop      = mult | div
expop       = Literal( "^" )
assign      = Literal( "=" )

expr = Forward()
atom = ( ( e | floatnumber | integer | ident ).setParseAction(pushFirst) |
         ( lpar + expr.suppress() + rpar )
       )
functions1arg = ['sqrt','sin','cos','tan','arcsin','arccos','arctan','log10']
functions2arg = ['arctan2']
function_tokens = map(Literal,functions1arg + functions2arg)

argument_list = (( Group(expr) + ZeroOrMore( Suppress( ',' ) + expr ) )).setResultsName('argumentlist')

function = (Or(function_tokens)+lpar+argument_list+rpar).setParseAction(arg_number,pushFirst)

primary = function | atom
factor = Forward()
factor << primary + ZeroOrMore( ( expop + factor ).setParseAction( pushFirst ) )
        
term = factor + ZeroOrMore( ( multop + factor ).setParseAction( pushFirst ) )
expr << term + ZeroOrMore( ( addop + term ).setParseAction( pushFirst ) )
bnf = Optional((ident + assign).setParseAction(assignVar)) + expr

pattern =  bnf + StringEnd()

# Map operator symbols to corresponding arithmetic operations
opn = { "+" : ( lambda a,b: a + b ),
        "-" : ( lambda a,b: a - b ),
        "*" : ( lambda a,b: a * b ),
        "/" : ( lambda a,b: a / b ),
        "^" : ( lambda a,b: a ** b ) }


# Recursive function that evaluates the stack
#------------------------------------------------------------------------------
def evaluateStack( s, type, snap ):
  op = s.pop()
  if debug_flag:
      print op
  if op in "+-*/^":
    op2 = evaluateStack( s, type, snap )
    op1 = evaluateStack( s, type, snap )
    return opn[op]( op1, op2 )
  elif op == "pi":
    return np.pi
  elif op == "e":
    return np.e
  elif op in functions1arg:
    operand = evaluateStack(s, type, snap)
    return getattr(np, op)(operand)
  elif op in functions2arg:
    operand2 = evaluateStack(s, type, snap)
    operand1 = evaluateStack(s, type, snap)
    return getattr(np, op)(operand1, operand2)
  elif op in _KnownQuantities():
    return UserQuantity(op).fetch(type, snap, "default")[1]
  elif op in snap.sim.GetIntAndFloatParameterKeys():
    return float(snap.sim.GetParam(op))
  elif re.search('^[-+]?[0-9]+$',op):
    return int( op )
  elif re.search('[-+]?[0-9]+[.][0-9]*',op):
    return float( op )
  else:
    raise ParseException("Unable to parse string: "  + op)
