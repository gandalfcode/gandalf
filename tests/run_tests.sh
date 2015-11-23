#!/bin/bash
set -v
ERRORS=0
for file in $(cat to_run)
do
    if [ -z "$MPICPP" ];
    then
        mpirun -np 4 ../bin/gandalf $file    
    else
        ../bin/gandalf $file
    fi
    if [ $? -ne 0 ];
    then
        ERRORS += 1
    fi
    echo $ERRORS
done
exit $ERRORS
