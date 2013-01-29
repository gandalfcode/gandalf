from collections import deque
import fnmatch
import os
from SphSim import SphSimulation
from SphSnap import SphSnapshot

# ============================================================================
# CLASS SIMBUFFER
# ============================================================================
class SimBuffer:


    # ========================================================================
    # SIMBUFFER __INIT__
    # ========================================================================
    def __init__ (self):

        self.Nsim = int(0)
        self.simlist = []
        self.snapshots = deque()
        self.maxmemory = 1024**2



    # ========================================================================
    # SIMBUFFER __INIT__
    # ========================================================================
    def _add_simulation(self,sim):

        self.simlist.append(sim)
        self.Nsim += 1
        
    def _addsnapshot(self, sim, snapshot):
        while 1:
            snapshotsize = snapshot.CalculateMemoryUsage()
            if snapshotsize > self.maxmemory:
                print "Warning: the requested snapshot can't fit inside the buffer memory. It will not be cached"
                raise BufferFull
            available_memory = self.maxmemory - snapshotsize
            if snapshotsize > available_memory:
                self._makespace()
            else:
                break
        self.snapshots.append(snapshot)
        snapshot.sim = sim
        sim.cachedsnapshots[snapshot.filename] = True
        
    def _makespace(self):
        ''' Deletes the first snapshot on the queue (that is, the first created)
        and deallocates its memory. Also the structure inside the cachedsnapshot dictionary
        needs to be updated
        '''
        to_be_deleted = self.snapshots.popleft()
        if to_be_deleted.allocated:
            to_be_deleted.DeallocateBufferMemory()
        #needs to update also the structure inside the cachedsnapshots dictionary
        index = self.simlist.index(to_be_deleted.sim)
        sim = self.simlist[index]
        sim.cachedsnapshots[to_be_deleted.filename] = False        
        
        
    def loadsim (self, run_id, fileformat = 'su', buffer_flag = 'cache'):
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
        sim = SphSimulation()
        self._add_simulation(sim)
        parameters = sim.simparams
        parameters.ReadParamsFile(paramfile)
        fileformat = parameters.stringparams["file_format"]
        
        #get the list of all files in the directory where the parameter file is
        paramfilepath = os.path.join(os.getcwd(),paramfile)
        dirname = os.path.dirname(paramfilepath)
        folderfiles = os.listdir(dirname)
        
        #search for all the files containing the given run_id string
        filetest = run_id_base + '.' + fileformat + '.?????'
        sim.cachedsnapshots = {}
        for filename in folderfiles:
            if fnmatch.fnmatch(filename, filetest):
                sim.cachedsnapshots[filename] = False              
        if len(sim.cachedsnapshots) == 0:
            print "Warning: no snapshots files found for simulation " + run_id        
        
        #depending on buffer_flag, caches all the snapshots or not 
        if buffer_flag == "store":
            for snapshotfile in sim.cachedsnapshots:
                #reads the snapshot from memory (we don't have the code for doing that yet!)
                snapshot = SphSnapshot()
                snapshot.Nsph = 100
                snapshot.AllocateBufferMemory()
                snapshot.filename = snapshotfile
                try:
                    self._addsnapshot(sim, snapshot)
                except BufferFull:
                    snapshot.DeallocateBufferMemory()
        elif buffer_flag == "cache":
            snapshotfile = run_id + '.' + fileformat + '.00000'
            snapshot = SphSnapshot()
            snapshot.Nsph = 100
            snapshot.AllocateBufferMemory()
            snapshot.filename = snapshotfile
            try:
                self._addsnapshot(sim, snapshot)
            except BufferFull:
                snapshot.DeallocateBufferMemory()
             
    
    def total_memory_usage(self):
        '''
        This function return the total memory usage of the buffer (in bytes).
        '''
        total_memory = 0
        for snapshot in self.snapshots:
            total_memory += snapshot.CalculateMemoryUsage()
        return total_memory
    
class BufferFull (Exception):
    pass    

from sys import getrefcount
import unittest
class Test(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        for filename in ['ciao.param', 'ciao.sf.00000', 'ciao.sf.00001']:
            file = open(filename,'w')
            file.close() 
        
    def testEmptyBuffer(self):
        buffer = SimBuffer()
        buffer.maxmemory = 0
        buffer.loadsim('ciao',buffer_flag = 'cache')
        self.assertEqual(buffer.total_memory_usage(),0)
        
    def testCacheBuffer(self):
        buffer = SimBuffer()
        buffer.maxmemory = 1024**2
        buffer.loadsim('ciao',buffer_flag = 'cache')
        self.assertEqual(buffer.total_memory_usage(),5200)
        self.assertEqual(len(buffer.snapshots), 1)
        snapshot = buffer.snapshots[0]
        self.assertGreater(getrefcount(snapshot), 2)
        filename = snapshot.filename
        self.assertEqual(buffer.simlist[0],snapshot.sim)
        self.assertEqual(buffer.simlist[0].cachedsnapshots[filename], True)
        buffer._makespace()
        self.assertEqual(buffer.simlist[0],snapshot.sim)
        self.assertEqual(buffer.simlist[0].cachedsnapshots[filename], False)
        self.assertEqual(getrefcount(snapshot), 2)
        self.assertEqual(buffer.total_memory_usage(),0)
        
    def testStoreBuffer(self):
        buffer = SimBuffer()
        buffer.maxmemory = 1024**2
        buffer.loadsim('ciao',buffer_flag = 'store')
        self.assertEqual(buffer.total_memory_usage(),10400)
        self.assertEqual(len(buffer.snapshots),2)
        
    def testRelativePath(self):
        try:
            os.mkdir('temp')
        except OSError:
            pass
        os.chdir('temp')
        buffer = SimBuffer()
        buffer.maxmemory = 1024**2
        buffer.loadsim('../ciao',buffer_flag = 'store')
        self.assertEqual(buffer.total_memory_usage(),10400)
        self.assertEqual(len(buffer.snapshots),2)
        os.chdir('..')
        os.rmdir('temp')