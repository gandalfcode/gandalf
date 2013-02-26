from seren.analysis.facade import *
import os
import sys

if __name__ == "__main__":
    #get the files in the same directory of the script
    testfiles = [""]
    for testfile in testfiles:
        try:
            newsim(testfile)
            run()
        except:
            #if an error has happened, print info about that and goes to the next test
            print "Unexpected error while evaluating test file ", testfile
            type, value, traceback = sys.exc_info()
            sys.excepthook(type, value, traceback)
            print "Going to the next test file..."
            sys.exc_clear()
            continue