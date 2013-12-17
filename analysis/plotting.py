#==============================================================================
#  plotting.py
#  ..
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
from multiprocessing import Process
from Queue import Empty
import os


#------------------------------------------------------------------------------
class PlottingProcess (Process):
    '''Contains the code for the plotting process.  The run method is the main
    method that gets executed. It contains an infinite loop, where the command
    queue is inspected;  if there are commands on it, they get executed.
    The method also does some bookkeeping, removing closed figures.
    After that, it pauses for a while until the queue is inspected again.
    The only way to exit from this loop is to get a 'STOP' string on the
    command queue.  The class contains also a couple of helper functions, and
    many structures to do the bookkeeping of figures, limits, quantities, ...
    (see the comments in the constructor).
    '''

    #--------------------------------------------------------------------------
    def __init__(self, queue, commands, completedqueue, globallimits, free):
        Process.__init__(self)
        self.queue = queue #queue for receiving commands and data
        self.commands = commands #list of commands to execute
        self.completedqueue = completedqueue #queue to signal main process that we are done
        self.commandsfigures = {} #dictionary that associate to each command id the figure, axis, product triplet produced by the command
        self.quantitiesfigures = {} #dictionary that associate to figure, axis, string (either 'x', 'y' or 'render') triplet the quantity, unit name double plotted
        self.globallimits = globallimits #dictionary that associate to each quantity the global limits
        self.axesimages = {} #dictionary that for each axis associate the corresponding image
        self.lastid = 0 #id of the last command received
        self.free = free #event variable that says if we are free


    #--------------------------------------------------------------------------
    def run(self):
        self.free.clear()
        import matplotlib.pyplot as plt
        import warnings
        warnings.filterwarnings("ignore", "matplotlib is currently using a non-GUI backend, so cannot show the figure")
        self.plt = plt
        self.ppid = os.getppid()
        
        # Main loop
        while 1:
            
            self.free.clear()
            
            if self.ppid != os.getppid():
                break
            
            self.mypause(0.01)
            
            # Reads data from the queue
            try:
                jobs = []
                while 1:
                    job = self.queue.get(block = False)
                    if job[0] == "STOP":
                        return 
                    jobs.append(job)
            except Empty:
                pass
            
            if len(jobs) == 0:
                self.free.set()
                continue
            else:
                self.free.clear()
                        
            self.remove_closed_figures()
            
            for job in jobs:
                command, data = job
                command.processCommand(self, data)    


    #--------------------------------------------------------------------------
    def command_in_list (self, id):
        '''Returns true if the given id is in the list,
        false otherwise'''
        for command in self.commands:
            if id == command.id:
                return True
        return False
    

    #--------------------------------------------------------------------------
    def mypause (self, interval):
        '''Pauses for a while, allowing the event loops of the figures to run,
        so that the user can interact with them (zoom, pan, ...).'''
        plt = self.plt
        backend = plt.rcParams['backend']
        if backend in plt._interactive_bk:
            figManager = plt._pylab_helpers.Gcf.get_active()
            if figManager is not None:
                canvas = figManager.canvas
                canvas.start_event_loop(interval)
                return
    
        # No on-screen figure is active, so sleep() is all we need.
        import time
        time.sleep(interval)

    
    # TODO: for efficiency reason, refactor this routine to be a callback
    # of when a window is closed
    #--------------------------------------------------------------------------
    def remove_closed_figures(self):
        '''Remove closed figures from the dictionary'''
        for index, commanddummy in enumerate(self.commands):
            fig, ax, line = self.commandsfigures[commanddummy.id]
            try:
                if not self.plt.fignum_exists(fig.number):
                    self.commandsfigures.pop(commanddummy.id)
                    self.commands.pop(index)
            # TODO: for efficiency, should remove the figures that have been
            # closed from the command list 
            except AttributeError:
                pass       
