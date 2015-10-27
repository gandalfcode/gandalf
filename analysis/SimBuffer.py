#==============================================================================
#  SimBuffer.py
#  Contains class definition for main simulation object buffer which is used
#  for storing and managing simulations in the python front-end.
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G. Rosotti
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
import fnmatch
import os
from swig_generated.SphSim import SimulationBase, SphSnapshotBase, Parameters, CodeTiming



#------------------------------------------------------------------------------
class SimBuffer:
    '''This class is intended both to act as a cache for the snapshots and to
provide routines to create new simulation instances.  Currently, it stores 
a list of the created simulations and a queue of the stored snapshots.
Each snapshot holds a reference to the parent simulation, while each
simulation holds a dictionary, listing all the snapshot files that have been
found on the disk, and whether they are cached.  To keep these structures
synchronized, it is better to use only the methods provided by this class to
modify its members (reading new snapshots and simulations and deleting them). 
The class is implemented as singleton, so it should never be instantiated; for
this reason, all of its methods are static.
'''

    # Initialise sim counters and empty lists
    # TODO : now hardwired to 1 GB.  Make more adaptive in future
    Nsim = int(0)
    simlist = []
    snapshots = []
    maxmemory = 1024**3 
    currentsim = -1
    

    #--------------------------------------------------------------------------
    @staticmethod
    def add_snapshot(snap, sim):
        ''' Adds a snapshot to a given simulation.  e.g. used by interactive
        run to dynamically add snapshots while the simulation is running.'''
        sim.snapshots.append(snap)
        SimBuffer.snapshots.append(snap)
        snap.sim = sim
        snap.live = False
        

    #--------------------------------------------------------------------------
    @staticmethod
    def _add_simulation(sim):
        '''Adds a new simulation object to the simulation buffer'''
        SimBuffer.simlist.append(sim)
        SimBuffer.Nsim += 1
        SimBuffer.currentsim = SimBuffer.Nsim - 1
        
    
    # TODO: improve the algorithm for making space/allocating space
    # TODO: there are still problems with the way memory is managed. Should 
    # check the reference counting before deallocating - we don't want to 
    # deallocate a snapshot that might be still in use.
    #--------------------------------------------------------------------------
    @staticmethod        
    def _findmemoryfor(snapshot):
        '''Add a new snapshot to the list of snapshots.  Checks for memory
        usage, and deallocates as many snapshots as needed for making space.
        Updates the reference to the parent simulation accordingly. Raises
        BufferFull if there isn\'t enough space for the snapshot being loaded. 
        '''
        snapshotsize = snapshot.CalculatePredictedMemoryUsage()
        if snapshotsize > SimBuffer.maxmemory:
            raise BufferFull("The requested snapshot can't fit inside the buffer memory")
        while 1:
            usedmemory = SimBuffer.total_memory_usage()
            available_memory = SimBuffer.maxmemory - usedmemory
            if snapshotsize > available_memory:
                SimBuffer._deallocateSnapshot(snapshot)
            else:
                break


    #--------------------------------------------------------------------------
    @staticmethod
    def _fillsnapshot(snapshot):
        '''Find free memory to store snapshot in buffer'''
        SimBuffer._findmemoryfor(snapshot)
        snapshot.ReadSnapshot(snapshot.sim.simparams.stringparams["out_file_form"])


    #--------------------------------------------------------------------------
    @staticmethod
    def _deallocateSnapshot(snapshottest):
        '''Deallocates a snapshot to make space.  The argument is not the
        snapshot to be deallocated, but the one that we are making space for.
        It is needed to check that we are not deallocating this snapshot
        itself.  The exact snapshot to deallocate depends on the specific
        caching algorithm; at the moment a LRU (least recently used) algorithm
        is used. The list of snapshots is sorted depending on the time at
        which the snapshot was used for the last time; the first snapshots to
        get deallocated are the ones that were used most time ago. Not that
        this technique is not scan resistant (but there are ways around that).
        '''
        for snapshot in sorted(SimBuffer.snapshots, key=lambda element: element.LastUsed):
            if snapshot.allocated:
                if snapshot != snapshottest:
                    snapshot.DeallocateBufferMemory()
                    return
        raise RuntimeError('SimBuffer._deallocateSnapshot: should never get to this line!!!!')
        

    #--------------------------------------------------------------------------
    @staticmethod
    def _destroy():
        '''Deallocates all snapshot memory in C++ memory buffer and clears
        the snapshots and simulation lists
        '''
        for snapshot in SimBuffer.snapshots:
            if snapshot.allocated:
                snapshot.DeallocateBufferMemory()
        SimBuffer.snapshots = []
        SimBuffer.simlist = []


    #--------------------------------------------------------------------------
    @staticmethod
    def newsim (paramfile=None, ndim=None, simtype=None):
        '''Creates a new simulation from specified parameter file, if given;
        otherwise, just used default values for the parameters, but needs ndim
        to be specified.  Returns the simulation created.
        '''
        
        # Create the parameter and timing object
        params = Parameters()
        timing = CodeTiming()
        
        # If a paramfile name was given, read it
        if paramfile is not None:
            params.ReadParamsFile(paramfile)
        if ndim is None:
            if paramfile is None:
                raise BufferException("You need to specify either a parameter file, or ndim and simtype")
            ndim = params.intparams["ndim"]
        if simtype is None:
            if paramfile is None:
                raise BufferException("You need to specify either a parameter file, or ndim and simtype")
            simtype = params.stringparams["sim"]
        sim = SimulationBase.SimulationFactory(ndim, simtype, params);
        sim.timing = timing
        SimBuffer._add_simulation(sim)
        sim.snapshots = []
        return sim


    #--------------------------------------------------------------------------
    @staticmethod
    def load_live_snapshot(sim):
        ''' Loads into memory the live snapshot of the simulation'''
        try:
            snap = sim.live
        except AttributeError:
            snap = SphSnapshotBase.SphSnapshotFactory("",sim,sim.simparams.intparams["ndim"])
            snap.sim = sim
            snap.live = True
        snap.CopyDataFromSimulation()
        snap.t = sim.t
        sim.live = snap
        sim.current = sim.live


    #--------------------------------------------------------------------------
    @staticmethod        
    def loadsim (run_id, fileformat=None, buffer_flag='cache'):
        '''Loads a simulation into the buffer.  Takes the run_id (can also be
        a relative path), the fileformat (read from the parameters file) and
        a string that flags how the buffer should behave with respect to
        memory management. Possible values are:
        
            cache (default)
                When simulation is loaded, only first snapshot is read. When a
                new snapshot is requested, it is saved into memory so to be
                able to quickly retrieve it.
            store
                When simulation is loaded, all snapshots are loaded in memory,
                up to filling the memory buffer available. It requires a
                longer start-up time, but then accessing the snapshot in
                memory will be very fast.
            nocache
                No caching of snapshots is done; the snapshots are always
                re-read from the disk.
        '''
        
        # Create parameters object and read parameters file using given run_id
        paramfile = run_id + '.param'
        run_id_base = os.path.basename(run_id)
        parameters = Parameters()
        parameters.ReadParamsFile(paramfile)
        ndim = parameters.intparams["ndim"]
        fileformat = parameters.stringparams["out_file_form"]
        simtype = parameters.stringparams["sim"]
        timing = CodeTiming()
        
        # Construct the simulation object and initialize it
        sim = SimulationBase.SimulationFactory(ndim, simtype, parameters);
        sim.timing = timing
        SimBuffer._add_simulation(sim)
        sim.ProcessParameters()
        
        # Get list of all files in the directory containing the parameter file
        paramfilepath = os.path.join(os.getcwd(),paramfile)
        dirname = os.path.dirname(paramfilepath)
        folderfiles = os.listdir(dirname)
        folderfiles.sort()
        
        # Search for all the files containing the given run_id string
        filetest = run_id_base + '.' + fileformat + '.?????'
        sim.snapshots = []
        for filename in folderfiles:
            if fnmatch.fnmatch(filename, filetest):
                snap = SphSnapshotBase.SphSnapshotFactory(os.path.join(dirname,filename),sim,ndim)
                snap.sim = sim
                snap.live = False
                sim.snapshots.append(snap)              
        if len(sim.snapshots) == 0:
            raise BufferException("Warning: no snapshots files found for simulation " + run_id)
            return
            
        sim.setup = True    

        print 'Found ',len(sim.snapshots),' snapshot file(s)'
        #for snapshot in sim.snapshots:
        #    print snapshot.filename      
        
        # Depending on buffer_flag, cache all the snapshots or not
        for i, snapshot in enumerate(sim.snapshots):
            snapshot.sim = sim
            SimBuffer.snapshots.append(snapshot)
            if buffer_flag == "store":
                SimBuffer._fillsnapshot(snapshot)
            elif buffer_flag == "cache":
                if i==0:
                    SimBuffer._fillsnapshot(snapshot)
                
        sim.current = sim.snapshots[0]


    #--------------------------------------------------------------------------
    @staticmethod  
    def total_memory_usage():
        '''Returns the total memory usage of the buffer (in bytes).'''
        total_memory = 0
        for snapshot in SimBuffer.snapshots:
            total_memory += snapshot.CalculateMemoryUsage()
        for sim in SimBuffer.simlist:
            try:
                total_memory += sim.live.CalculateMemoryUsage()
            except AttributeError:
                pass
        return total_memory


    #--------------------------------------------------------------------------
    @staticmethod  
    def get_current_sim():
        '''Returns the current simulation'''
        return SimBuffer.simlist[SimBuffer.currentsim]


    #--------------------------------------------------------------------------
    @staticmethod
    def get_current_sim_no():
        '''Returns the i.d. number of the current simulation'''
        return SimBuffer.currentsim


    #--------------------------------------------------------------------------
    @staticmethod
    def set_current_sim_no(simno):
        '''Sets the current simulation to the requested number'''
        
        # Checks that the simulation exists
        sim = SimBuffer.get_sim_no(simno)
        SimBuffer.currentsim = simno
        return sim
    

    #--------------------------------------------------------------------------
    @staticmethod
    def get_sim_no(no):
        '''Returns the simulation specified by the number'''
        try:
            sim = SimBuffer.simlist[no]
        except IndexError:
            raise BufferException ("The specified simulation number does not exist")
        return sim


    #--------------------------------------------------------------------------
    @staticmethod  
    def get_current_snapshot_by_sim(sim):
        '''Queries the simulation given and returns its current snapshot'''
        try:
            snap = sim.current
        except AttributeError:
            SimBuffer.load_live_snapshot(sim)
            snap=sim.current
        if (not snap.allocated):
            SimBuffer._fillsnapshot(snap)
        return snap


    #--------------------------------------------------------------------------
    @staticmethod      
    def get_current_snapshot():
        '''Returns the current snapshot of the current simulation'''
        sim = SimBuffer.get_current_sim()
        return SimBuffer.get_current_snapshot_by_sim(sim)


    #--------------------------------------------------------------------------
    @staticmethod
    def get_live_snapshot_current():
        '''Returns the live snapshot of the current simulation, if exists'''
        sim = SimBuffer.get_current_sim()
        try:
            return sim.live
        except AttributeError:
            raise BufferException("The current simulation does not have a live snapshot")


    #--------------------------------------------------------------------------
    @staticmethod
    def get_live_snapshot_sim(sim):
        '''Returns the live snapshot of the specified simulation'''
        try:
            return sim.live
        except AttributeError:
            raise BufferException("The specified simulation does not have a live snapshot")
    

    #--------------------------------------------------------------------------
    @staticmethod
    def get_snapshot_number_sim(sim, no):
        '''Queries the buffer for the given snapshot number and simulation'''
        no=int(no)
        try:
            snap = sim.snapshots[no]
        except IndexError:
            raise BufferException ("The selected snapshot does not exist")
        if (not snap.allocated):
            SimBuffer._fillsnapshot(snap)
        return snap


    #--------------------------------------------------------------------------
    @staticmethod  
    def set_current_snapshot_number_sim(sim, no):
        '''Queries the buffer for the given snapshot number and simulation
        and sets it to the current snapshot of the given simulation'''
        snap = SimBuffer.get_snapshot_number_sim(sim, no)
        sim.current = snap
        return snap


    #--------------------------------------------------------------------------
    @staticmethod
    def set_current_snapshot_number (no):
        '''Queries the buffer for the given snapshot of the current simulation
        and sets it to the current snapshot of the current simulation'''
        snap = SimBuffer.get_snapshot_number(no)
        sim = SimBuffer.get_current_sim()
        sim.current = snap
        return snap


    #--------------------------------------------------------------------------
    @staticmethod
    def get_snapshot_number (no):
        '''Returns the snapshot number from the current simulation'''
        sim = SimBuffer.get_current_sim()
        return SimBuffer.get_snapshot_number_sim(sim, no)


    #--------------------------------------------------------------------------
    @staticmethod
    def get_no_next_snapshot():
        '''Returns the number of the next snapshot, without modifying it'''
        sim = SimBuffer.get_current_sim()
        next_index = sim.snapshots.index(sim.current)+1
        if next_index >= len(sim.snapshots):
            raise BufferException ("Reached the last snapshot")
        return next_index


    #--------------------------------------------------------------------------
    @staticmethod
    def get_no_previous_snapshot():
        '''Returns the number of the previous snapshot, without modifying it'''
        sim = SimBuffer.get_current_sim()
        previous_index = sim.snapshots.index(sim.current)-1
        if previous_index < 0:
            raise BufferException ("Reached the first snapshot")
        return previous_index


    #--------------------------------------------------------------------------
    @staticmethod
    def get_snapshot_extended(sim, snapno):
        '''Returns the desired snapshot of the given simulation object.
        Handles the case where snapno is not only a number, but can be also a
        string (either "current" or "live");  also handles the case where
        sim is the "current" string.
        '''
        if sim=="current":
            sim_obj = SimBuffer.get_current_sim()
        else:
            sim_obj = sim
        
        if snapno == "current":
            snap = SimBuffer.get_current_snapshot_by_sim(sim_obj)
        elif snapno == "live":
            snap = SimBuffer.get_live_snapshot_sim(sim_obj)
        elif isinstance(snapno,int) or isinstance(snapno, basestring):
            snap = SimBuffer.get_snapshot_number_sim(sim_obj, snapno)
        elif isinstance(snapno,SphSnapshotBase):
            snap = snapno
            
        return snap


    #--------------------------------------------------------------------------
    @staticmethod
    def get_sim_iterator(sim):
        return SimIterator(sim)


    #--------------------------------------------------------------------------
    @staticmethod
    def get_next_snapshot_from_object(snap):
        '''Returns next snapshot of current simulation'''
        sim = snap.sim
        snapno = sim.snapshots.index(snap)
        next_index = snapno + 1
        if next_index >= len(sim.snapshots):
            raise BufferException("Reached the last snapshot")
        snapret = sim.snapshots[next_index]
        if (not snapret.allocated):
            SimBuffer._fillsnapshot(snapret)
        return snapret


    #--------------------------------------------------------------------------
    @staticmethod
    def get_previous_snapshot_from_object(snap):
        '''Returns previous snapshot of current simulation'''
        sim = snap.sim
        snapno = sim.snapshots.index(snap)
        previous_index = snapno-1
        if previous_index<0:
            raise BufferException("Reached the first snapshot")
        snapret = sim.snapshots[previous_index]
        if (not snapret.allocated):
            SimBuffer._fillsnapshot(snapret)
        return snapret


#------------------------------------------------------------------------------
class BufferException(Exception):   
    pass


#------------------------------------------------------------------------------
class BufferFull(BufferException):
    pass


#------------------------------------------------------------------------------
class SimIterator:
    
    def __init__(self, sim):
        self._sim = sim
        self.snapnumber = 0
    
    def __iter__(self):
        return self
    
    def next(self):
        try:
            result = SimBuffer.get_snapshot_number_sim(self._sim, self.snapnumber)
        except BufferException:
            raise StopIteration
        
        self.snapnumber += 1
        return result
