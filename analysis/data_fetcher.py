#==============================================================================
#  data_fetcher.py
#  ..
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
from swig_generated.SphSim import UnitInfo
from facade import SimBuffer

direct = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az',
          'm', 'h', 'rho', 'u', 'dudt','sma','ecc', 'mbin','period',
          'qbin']

time_fetchers={}

derived_fetchers = {}
function_fetchers = {}

#------------------------------------------------------------------------------
def _KnownQuantities():
  '''Return the list of the quantities that we know'''
  return derived_fetchers.keys()+direct+function_fetchers.keys()

#------------------------------------------------------------------------------
def UserQuantity(quantity):
    '''Given a quantity, return a fetcher that we can query to get that quantity'''
    if quantity in direct:
        return DirectDataFetcher(quantity)
    elif quantity in derived_fetchers:
        return derived_fetchers[quantity]
    elif quantity in function_fetchers:
        return function_fetchers[quantity]
    else:
        raise Exception("UserQuantity: we don't know how to compute " + quantity)
    
from formula_parser import evaluateStack, exprStack, varStack, pattern    


#------------------------------------------------------------------------------
def CreateUserQuantity(name, formula, unitlabel='', unitname='',scaling_factor=1, label=''):
    '''Given a mathematical formula, build a data fetcher from it.
The quantity is given a name, which can now be used in plots and in other 
formulae.  When you construct a quantity, you can rely on one of the units we 
provide, in which case you can just pass as the scaling_factor parameter the 
name of the unit you want inside the SimUnits class. For example, if your unit 
has dimensions of acceleration, you can pass 'a' as the scaling_factor 
parameter. Doing this allows the unit system to work seamlessly when plotting 
(i.e., you can specify the units you want the plot in).  Alternatively, you can 
build your own unit passing a numerical value for the scaling_factor, 
a unitname and a latex label.  In this case, however, no rescaling is possible, 
as the unit system does not know how to rescale your unit.    
'''
    if isinstance(formula, basestring):
        fetcher = FormulaDataFetcher(name, formula, unitlabel, unitname, scaling_factor, label)
        derived_fetchers[name]=fetcher
    else:
        fetcher = FunctionFetcher(name,formula,unitlabel,unitname,scaling_factor,label)
        function_fetchers[name] = fetcher
    return fetcher

  
#------------------------------------------------------------------------------
def CreateTimeData(name, function, *args, **kwargs):
    '''Given a function that takes a snapshot as input, construct a 
FunctionTimeDataFetcher object from it and register it.
    Return the FunctionTimeDataFetcher newly constructed.'''
    fetcher = FunctionTimeDataFetcher(function, *args, **kwargs)
    time_fetchers [name] = fetcher
    return fetcher

#------------------------------------------------------------------------------
def TimeData(quantity, id=None):
    '''Given a quantity, return the FunctionTimeDataFetcher object that we can query'''
    try:
        fetcher = time_fetchers[quantity]
    except KeyError:
        #check if the quantity exists in UserQuantities
        try:
            quantity_fetcher = UserQuantity(quantity)
        except Exception:
            raise KeyError("TimeData: we do not know how to compute " + quantity)
        #if the quantity exists, we can now build a data fetcher on the fly
        if id==None:
            raise KeyError("TimeData: you didn't specify the id of the particle for plotting quantity " + quantity)
        id=int(id)
        from compute import particle_data
        name='part_' + quantity + '_' + str(id)
        fetcher = CreateTimeData(name,particle_data,quantity=quantity,id=id)
        
    return fetcher


#------------------------------------------------------------------------------
def check_requested_quantity(quantity, snap):
    '''Check the requested quantity exists, depending on the dimensionality of the snapshot.
    Also return information about the kind of the quantity (direct or derived)'''
    
    #check dimensionality
    twod = ('y', 'vy', 'ay')
    threed = ('z', 'vz', 'az')
    minus3 = quantity in threed+('r','theta') and snap.ndim<3
    minus2 = quantity in twod+('R','phi') and snap.ndim<2
    if minus3 or minus2:
        raise Exception("Error: you requested the quantity " + quantity + ", but the simulation is only in " + str(snap.ndim) + " dims")
    
    #if it's not a live snapshot, check that we are not requesting quantities defined only for live snapshots
    if not snap.live:
        if quantity in ('ax', 'ay', 'az'):
            raise Exception ("Error: accelerations are available only for live snapshots")
        elif quantity in ('dudt',):
            raise Exception ("Error: dudt is available only for live snapshots")
    
    #check that we know how to compute the quantity
    if quantity in direct:
        return "direct"
    elif quantity in derived_fetchers:
        return "derived"
    elif quantity in function_fetchers:
        return "function"
    else:
        raise Exception("We don't know how to compute " + quantity)
    

#------------------------------------------------------------------------------
class DirectDataFetcher:
    
    quantitylabels = {'x': 'x', 'y': 'y', 'z': 'z', 'rho': '$\\rho$',
                      'vx': '$v_x$', 'vy': '$v_y$', 'vz': '$v_z$', 
                      'ax': '$a_x$', 'ay': '$a_y$', 'az': '$a_z$',
                      'm': 'm', 'h': 'h', 'u': 'u', 't': 't',
                      'sma': 'Semi-major axis', 'ecc': 'Orbital eccentricity',
                      'mbin': 'Binary mass',
                      'period': 'Orbital period', 'qbin': 'Mass ratio'
                      }
    
    def __init__(self, quantity):
        
        if quantity not in direct:
            raise Exception ("Error: the quantity" + quantity + " is not a direct quantity!")
        self._quantity = quantity
        
    def fetch(self, type="default", snap="current", unit="default"):
        
        if snap=="current":
            snap=SimBuffer.get_current_snapshot()
            
        kind = check_requested_quantity(self._quantity, snap)
        if kind != "direct":
            raise Exception ("Error: the quantity" + quantity + " is not a direct quantity!")
        
        return snap.ExtractArray(self._quantity, type, unit) + [self.quantitylabels[self._quantity]]


#------------------------------------------------------------------------------
class FormulaDataFetcher:
    
    def __init__(self, name, formula, unitlabel='', unitname='',scaling_factor=1,label=''):
        self.scaling_factor=scaling_factor
        self._name = name
        exprStack[:]=[]
        pattern.parseString(formula)
        self._stack = list(exprStack)
        self.unitinfo = UnitInfo()
        self.unitinfo.label=unitlabel
        self.unitinfo.name=unitname
        self.label=label
        
    def fetch(self, type="default", snap="current", unit="default"):
        
        if snap=="current":
            snap=SimBuffer.get_current_snapshot()
        
        result = evaluateStack(list(self._stack), type, snap)

        if isinstance(self.scaling_factor,basestring):
            try:
                unitobj=getattr(snap.sim.simunits, self.scaling_factor)
                if unit=="default":
                    unit=unitobj.outunit
                scaling_factor=unitobj.OutputScale(unit)
                self.unitinfo.name=unit
                self.unitinfo.label=unitobj.LatexLabel(unit)
            except AttributeError:
                raise AttributeError("Sorry, we do not know the unit " + self.scaling_factor)
        else:
            scaling_factor=float(self.scaling_factor)
        return self.unitinfo, result, scaling_factor, self.label


#------------------------------------------------------------------------------
class FunctionTimeDataFetcher:
    
    def __init__(self, function, unitlabel="", unitname="",scaling_factor =1,label='',**kwargs):
        self._function = function
        self._kwargs = kwargs
        self.unitlabel=unitlabel
        self.unitname=unitname
        self.unitinfo=UnitInfo()
        self.scaling_factor=scaling_factor
        self.label=label
    
    def fetch(self, sim="current", type="default", unit="default"):
        
        if sim=="current":
            sim=SimBuffer.get_current_sim()
        elif isinstance(sim,int):
            sim=SimBuffer.get_sim_no(sim)
        
        iterator = SimBuffer.get_sim_iterator(sim)
        results = map(lambda snap: self._function(snap,type=type,unit=unit,**self._kwargs),iterator)
        if isinstance(results[0],list) or isinstance(results[0],tuple):
            results_zipped = zip(*results)
        
            values = np.asarray(results_zipped[1])
        
            return results_zipped[0][0], values, results_zipped[2][0], results_zipped[3][0]
        
        else:
            values=np.asarray(results)
            if isinstance(self.scaling_factor,basestring):
                try:
                    unitobj=getattr(sim.simunits, self.scaling_factor)
                    if unit=="default":
                        unit=unitobj.outunit
                    scaling_factor=unitobj.OutputScale(unit)
                    self.unitinfo.name=unit
                    self.unitinfo.label=unitobj.LatexLabel(unit)
                except AttributeError:
                    raise AttributeError("Sorry, we do not know the unit " + self.scaling_factor)
            else:
                scaling_factor=float(self.scaling_factor)
                self.unitinfo.label=self.unitlabel
                self.unitinfo.name=self.unitname
                
        return self.unitinfo, values, scaling_factor, self.label
            
    
class FunctionFetcher:
    
    def __init__(self, name, function, unitlabel='', unitname='',scaling_factor=1,label='', **kwargs):
        self._name=name
        self._function = function
        self.unitlabel = unitlabel
        self.unitname = unitname
        self.scaling_factor = scaling_factor
        self.label = label
        self._kwargs = kwargs
        self.unitinfo = UnitInfo()
        
    def fetch(self,type="default", snap="current", unit="default"):
        
        if snap=="current":
            snap=SimBuffer.get_current_snapshot()
            
        data=self._function(snap,type=type,unit=unit,**self._kwargs)
        print data
        if isinstance(data,list) or isinstance(data,tuple):
            return data
        else:
            if isinstance(self.scaling_factor,basestring):
                try:
                    unitobj=getattr(snap.sim.simunits, self.scaling_factor)
                    if unit=="default":
                        unit=unitobj.outunit
                    scaling_factor=unitobj.OutputScale(unit)
                    self.unitinfo.name=unit
                    self.unitinfo.label=unitobj.LatexLabel(unit)
                except AttributeError:
                    raise AttributeError("Sorry, we do not know the unit " + self.scaling_factor)
            else:
                scaling_factor=float(self.scaling_factor)
                self.unitinfo.label=self.unitlabel
                self.unitinfo.name=self.unitname
        
        return self.unitinfo, data, scaling_factor, self.label

def get_time_snapshot(snap, type=None, unit="default"):
    '''Return the time of a given snapshot (plus all the needed information about units and labels)'''
    sim=snap.sim
    time_unit_obj = sim.simunits.t
    if unit=="default":
        unit = time_unit_obj.outunit
    unitinfo = UnitInfo()
    unitinfo.name = unit
    unitinfo.label = time_unit_obj.LatexLabel(unit)
    scaling_factor = time_unit_obj.OutputScale(unit)
    label = 't'
    return unitinfo, snap.t, scaling_factor, label
