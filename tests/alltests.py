from seren.analysis.facade import *
import os
import sys

setupsfolder = '../setups'

if __name__ == "__main__":
    setupfiles = os.listdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),setupsfolder))
    print "Found the following setup files: "
    print setupfiles
    for setupfile in setupfiles:
        try:
            newsim(setupfile)
            run()
        except:
            print "Unexpected error while evaluating setup file ", setupfile
            type, value, traceback = sys.exc_info()
            sys.excepthook(type, value, traceback)
            print "Going to the next setup file..."
            sys.exc_clear()
            continue