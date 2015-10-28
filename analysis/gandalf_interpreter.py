#!/usr/bin/env python
#==============================================================================
# gandalf_interpreter.py
# Python program for initialising interactive interpreter for using
# gandalf inside python.
#==============================================================================
import cmd2 as cmd
import facade
import inspect
import types

toexcludefunctions = ['Manager', 'Queue', 'handle', 'init', 
                      'sigint', 'cleanup', 'get_sim_no',
                      'to_list','to_bool','Event','_load',
                      '_relative_load']

'''This module contains the implementation of the interpreter of our special
mini-language'''


#------------------------------------------------------------------------------
def excludefunctions(function):
    if function[0] in toexcludefunctions:
        return False
    return True


#------------------------------------------------------------------------------
class WrappedFunction:
    '''This class creates a wrapper around each one of the functions in facade.
The wrapper code is responsible for taking the arguments passed on the command
line, and converting them to proper arguments for a python function call.
'''
    def __init__(self, function):
        self.function = function
        self.__doc__ = function.__doc__
    
    def __call__ (self, line):
        args = line.split()
        dict = {}
        required_arguments = []
        for arg in args:
            potkeywords = arg.split('=')
            if len(potkeywords) == 2:
                kw, value = potkeywords
                dict[kw] = value
            elif len(potkeywords) >2:
                raise Exception
            else:
                required_arguments.append(arg)
        try:
            self.function(*required_arguments, **dict)
        except Exception as e:
            facade.handle(e)


#------------------------------------------------------------------------------
class Interpreter(cmd.Cmd):
    '''This class implements the interpreter proper. It gets all the functions 
defined in facade, filters out some name that we do not want (because it\'s
reserved), and defines a command for each one of them. Also imports the
documentation string for each one.
''' 
    abbreviations = {
                      'next': 'n',
                      'previous': 'p',
                      'snap': 's',
                      'help': 'h'
                      }
    
    def __init__(self):
        cmd.Cmd.__init__(self)
        
        # Retrieve the functions in facade and filter out some
        functions = inspect.getmembers(facade, inspect.isfunction)
        functions=filter(excludefunctions, functions)
        
        # Create the wrapper for each one of them
        for function in functions:
            wrapper = WrappedFunction(function[1])
            setattr(self, 'do_'+function[0], wrapper)
        
        # Generate the abbreviated versions of the commands
        for key in self.abbreviations:
            abbr = self.abbreviations[key]
            name_abbr='do_'+abbr
            name_long='do_'+key
            setattr(self, name_abbr, getattr(self, name_long))
            
        # Convince facade that we are in interactive mode
        facade.interactive = True

        # Set command ine prompt to 'gandalf'
        self.prompt = 'gandalf > '
        
    def get_names(self):
        names = dir(self)
        names.remove('do_help')
        return names
        
if __name__ == "__main__":
    interpreter = Interpreter()
    interpreter.cmdloop()
