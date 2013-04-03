#!/usr/bin/env python

import cmd
import facade
import inspect
import types

toexcludefunctions = ['Manager', 'Queue', 'handle', 'init', 
                      'sigint', 'cleanup', 'get_sim_no']

def excludefunctions(function):
    if function[0] in toexcludefunctions:
        return False
    return True

class WrappedFunction:
    def __init__(self, function):
        self.function = function
        self.__doc__ = function.__doc__
    
    def __call__ (self, line):
        args = line.split()
        dict={}
        for arg in args:
            potkeywords = arg.split('=')
            if len(potkeywords) == 2:
                kw, value = potkeywords
            elif len(potkeywords) >2:
                raise Exception
            else:
                continue
            dict[kw]=value
            args.remove(arg)
        self.function(*args, **dict)

class Interpreter(cmd.Cmd):
    def __init__(self):
        cmd.Cmd.__init__(self)
        
        #retrieve the functions in facade and filter out some
        functions = inspect.getmembers(facade, inspect.isfunction)
        functions=filter(excludefunctions, functions)
        
        #create the wrapper for each one of them
        for function in functions:
            wrapper = WrappedFunction(function[1])
            setattr(self, 'do_'+function[0], wrapper)
        
    def do_EOF(self, line):
        '''Quits the interpreter'''
        print
        return True
        
    def get_names(self):
        names = dir(self)
        names.remove('do_help')
        return names
        
if __name__ == "__main__":
    interpreter = Interpreter()
    interpreter.cmdloop('Welcome to PySeren++')