import fnmatch
import os
from swig_generated.SphSim import SimulationBase, SphSnapshotBase, Parameters

# ============================================================================
# CLASS SIMBUFFER
# ============================================================================
class SimBuffer:
    '''
    This class is intended both to act as a cache for the snapshots and to
    provide routines to create new simulation instances.
    At the moment, it stores a list of the created simulations and a queue of the
    stored snapshots. Each snapshot holds a reference to the parent simulation,
    while each simulation holds a dictionary, listing all the snapshot files that
    have been found on the disk, and whether they are cached.
    To keep these structures synchronized, it is better to use only the methods
    provided by this class to modify its members (reading new snapshots and
    simulations and deleting them). 
    The class is implemented as singleton, so it should never be instantiated; for
    this reason, all of its methods are static.
    '''

    Nsim = int(0)
    simlist = []
    snapshots = []
    maxmemory = 1024**3 #now hardwired to 1 GB
    currentsim = -1


    # ========================================================================
    # SIMBUFFER __INIT__
    # ========================================================================
    @staticmethod
    def _add_simulation(sim):

        SimBuffer.simlist.append(sim)
        SimBuffer.Nsim += 1
        SimBuffer.currentsim = SimBuffer.Nsim - 1
    
    # TODO: improve the algorithm for making space/allocating space
    # TODO: there are still problems with the way memory is managed. Should check
    #the reference counting before deallocating - we don't want to deallocate a snapshot
    #that might be still in use
    @staticmethod        
    def _findmemoryfor(snapshot):
        '''
        Add a new snapshot to the list of snapshots.
        Checks for memory usage, and deallocates as many snapshots as needed
        for making space. Updates the reference
        to the parent simulation accordingly. Raises BufferFull if there isn't enough space
        for the snapshot being loaded. 
        '''
        snapshotsize = snapshot.CalculateMemoryUsage()
        if snapshotsize > SimBuffer.maxmemory:
            raise BufferFull("The requested snapshot can't fit inside the buffer memory")
        while 1:
            usedmemory = SimBuffer.total_memory_usage()-snapshotsize
            available_memory = SimBuffer.maxmemory - usedmemory
            if snapshotsize > available_memory:
                SimBuffer._deallocateSnapshot(snapshot)
            else:
                break
    
    @staticmethod
    def _fillsnapshot(snapshot):
        snapshot.ReadSnapshot(snapshot.sim.simparams.stringparams["in_file_form"])
        SimBuffer._findmemoryfor(snapshot)
    
    @staticmethod
    def _deallocateSnapshot(snapshottest):
        ''' Deallocates a snapshot to make space.
        The argument is not the snapshot to be deallocated, but the one that we are making space for.
        It is needed to check that we are not deallocating this snapshot itself.
        The exact snapshot to deallocate depends on the specific caching algorithm; at the moment
        a LRU (least recently used) algorithm is used. The list of snapshots is sorted depending on the
        time at which the snapshot was used for the last time; the first snapshots to get deallocated are
        the ones that were used most time ago. Not that this technique is not scan resistant (but there
        are ways around that).
        '''
        for snapshot in sorted(SimBuffer.snapshots, key=lambda element: element.LastUsed):
            if snapshot.allocated:
                if snapshot != snapshottest:
                    snapshot.DeallocateBufferMemory()
                    return
        raise RuntimeError('SimBuffer._deallocateSnapshot: should never get to this line!!!!')
        
    
    @staticmethod
    def _destroy():
        '''Deallocates all the snapshots and clears the snapshots and simulation lists'''
        for snapshot in SimBuffer.snapshots:
            if snapshot.allocated:
                snapshot.DeallocateBufferMemory()
        SimBuffer.snapshots = []
        SimBuffer.simlist = []
    
    @staticmethod
    def newsim (paramfile=None, ndim=None):
        '''
        This method creates a new simulation from the specified parameter file, if given;
        otherwise, just used default values for the parameters, but needs ndim to be specified.
        Returns the simulation created.
        '''
        #create the parameter object
        params = Parameters()
        #if a paramfile name was given, read it
        if paramfile is not None:
            params.ReadParamsFile(paramfile)
        if ndim is None:
            if paramfile is None:
                raise BufferException("You need to specify either the number of dimensions, either the parameter file")
            ndim = params.intparams["ndim"]
        sim = SimulationBase.SimulationFactory(ndim, params);
        SimBuffer._add_simulation(sim)
        sim.snapshots = []
        return sim
        
    @staticmethod
    def load_live_snapshot(sim):
        ''' This method loads into memory the live snapshot of the simulation'''
        try:
            snap = sim.live
        except AttributeError:
            snap = SphSnapshotBase.SphSnapshotFactory("",sim,sim.simparams.intparams["ndim"])
            snap.sim = sim
            snap.live = True
            #TODO: check this is correct
            snap.units = sim.simunits
        snap.CopyDataFromSimulation()
        snap.t = sim.t
        sim.live = snap
        sim.current = sim.live
    
    @staticmethod        
    def loadsim (run_id, fileformat = 'ascii', buffer_flag = 'cache'):
        '''
        This method loads a simulation into the buffer.
        
        Takes the run_id (can also be a relative path), the fileformat (read from
        the parameters file) and a string that flags how the buffer should behave
        with respect to memory management. Possible values are:
            
            cache
                (the default) When the simulation is loaded, only the first snapshot is read.
                When a new snapshot is requested, it is saved into memory so to be able to quickly
                retrieve it.
            store
                When the simulation is loaded, all the snapshots are loaded in memory, up to filling
                the memory buffer available. It requires a longer start-up time, but then accessing
                the snapshot in memory will be very fast.
            nocache
                No caching of the snapshots is done; the snapshots are always re-read from the disk.
        
        '''
        
        #construct the simulation object and initialize it
        paramfile = run_id + '.param'
        run_id_base = os.path.basename(run_id)
        parameters = Parameters()
        parameters.ReadParamsFile(paramfile)
        ndim=parameters.intparams["ndim"]
        sim = SimulationBase.SimulationFactory(ndim, parameters);
        SimBuffer._add_simulation(sim)
        sim.ProcessParameters()
        fileformat = parameters.stringparams["in_file_form"]
        
        #get the list of all files in the directory where the parameter file is
        paramfilepath = os.path.join(os.getcwd(),paramfile)
        dirname = os.path.dirname(paramfilepath)
        folderfiles = os.listdir(dirname)
        folderfiles.sort()
        
        #search for all the files containing the given run_id string
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
            
        for snapshot in sim.snapshots:
            print snapshot.filename      
        
        #depending on buffer_flag, caches all the snapshots or not
        for i, snapshot in enumerate(sim.snapshots):
            snapshot.sim = sim
            SimBuffer.snapshots.append(snapshot)
            if buffer_flag == "store":
                SimBuffer._fillsnapshot(snapshot)
            elif buffer_flag == "cache":
                if i==0:
                    SimBuffer._fillsnapshot(snapshot)
                
        sim.current = sim.snapshots[0]
    
    @staticmethod  
    def total_memory_usage():
        '''
        This function return the total memory usage of the buffer (in bytes).
        '''
        total_memory = 0
        for snapshot in SimBuffer.snapshots:
            total_memory += snapshot.CalculateMemoryUsage()
        for sim in SimBuffer.simlist:
            try:
                total_memory += sim.live.CalculateMemoryUsage()
            except AttributeError:
                pass
        return total_memory
    
    @staticmethod  
    def get_current_sim():
        '''This function returns the current simulation'''
        return SimBuffer.simlist[SimBuffer.currentsim]
    
    @staticmethod
    def get_current_sim_no():
        '''This function returns the number of the current simulation'''
        return SimBuffer.currentsim
    
    @staticmethod
    def set_current_sim_no(simno):
        '''This function sets the current simulation to the given number'''
        #Checks that the simulation exists
        sim=SimBuffer.get_sim_no(simno)
        SimBuffer.currentsim=simno
        return sim
    
    
    @staticmethod
    def get_sim_no(no):
        '''This function returns the simulation specified by the number'''
        try:
            sim = SimBuffer.simlist[no]
        except IndexError:
            raise BufferException ("The specified simulation number does not exist")
        return sim
    
    @staticmethod  
    def get_current_snapshot_by_sim(sim):
        '''This function queries the simulation given and returns its current snapshot'''
        try:
            snap = sim.current
        except AttributeError:
            SimBuffer.load_live_snapshot(sim)
            snap=sim.current
        if (not snap.allocated):
            SimBuffer._fillsnapshot(snap)
        return snap
    
    @staticmethod      
    def get_current_snapshot():
        '''This function returns the current snapshot of the current simulation'''
        sim = SimBuffer.get_current_sim()
        return SimBuffer.get_current_snapshot_by_sim(sim)
    
    @staticmethod
    def get_live_snapshot_current():
        '''This function returns the live snapshot of the current simulation, if exists'''
        sim = SimBuffer.get_current_sim()
        try:
            return sim.live
        except AttributeError:
            raise BufferException("The current simulation does not have a live snapshot")
        
    @staticmethod
    def get_live_snapshot_sim(sim):
        '''This function returns the live snapshot of the specified simulation'''
        try:
            return sim.live
        except AttributeError:
            raise BufferException("The specified simulation does not have a live snapshot")
    
    
    @staticmethod
    def get_snapshot_number_sim(sim, no):
        '''This function queries the buffer for the given snapshot number of the given simulation'''
        try:
            snap = sim.snapshots[no]
        except IndexError:
            raise BufferException ("The selected snapshot does not exist")
        if (not snap.allocated):
            SimBuffer._fillsnapshot(snap)
        return snap
    
    @staticmethod  
    def set_current_snapshot_number_sim(sim, no):
        '''This function queries the buffer for the given snapshot number of the given simulation
        and sets it to the current snapshot of the given simulation'''
        snap = SimBuffer.get_snapshot_number_sim(sim, no)
        sim.current = snap
        return snap
    
    @staticmethod
    def set_current_snapshot_number (no):
        '''This function queries the buffer for the given snapshot of the current simulation
        and sets it to the current snapshot of the current simulation'''
        snap = SimBuffer.get_snapshot_number(no)
        sim = SimBuffer.get_current_sim()
        sim.current = snap
        return snap
    
    @staticmethod
    def get_snapshot_number (no):
        '''This function returns the snapshot number from the current simulation'''
        sim = SimBuffer.get_current_sim()
        return SimBuffer.get_snapshot_number_sim(sim, no)
    
    @staticmethod
    def get_no_next_snapshot():
        '''This function returns the number of the next snapshot, without modifying it'''
        sim = SimBuffer.get_current_sim()
        next_index = sim.snapshots.index(sim.current)+1
        if next_index >= len(sim.snapshots):
            raise BufferException ("Reached the last snapshot")
        return next_index
    
    @staticmethod
    def get_no_previous_snapshot():
        '''This function returns the number of the previous snapshot, without modifying it'''
        sim = SimBuffer.get_current_sim()
        previous_index = sim.snapshots.index(sim.current)-1
        if previous_index < 0:
            raise BufferException ("Reached the first snapshot")
        return previous_index
    
    @staticmethod
    def get_sim_iterator(sim):
        return SimIterator(sim)
    
         
class BufferException( Exception):   
    pass

class BufferFull (BufferException):
    pass

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