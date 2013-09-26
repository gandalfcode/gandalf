from gandalf.analysis.facade import *
import os
import sys
import subprocess

if __name__ == "__main__":
    pyfiles = ["adsod-L1error.py","adsod-energyerror.py","soundwave-L1error.py"]
    for pyfile in pyfiles:
        try:
            return_code = subprocess.call(["python",pyfile])
            print "return_code : ",return_code
        except:
            print "Unexpected error while evaluating test script ",pyfile
            type, value, traceback = sys.exc_info()
            sys.excepthook(type, value, traceback)
            print "Going to the next test file..."
            sys.exc_clear()
            continue


    sys.exit()

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
