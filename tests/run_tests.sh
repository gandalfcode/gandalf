#!/bin/bash

CODE=0
for file in $(cat to_run)
do
    if [ -z "$VAR" ];
    then
        mpirun -np 4 ../bin/gandalf $file    
    else
        ../bin/gandalf $file
    fi
    CODE=$? || $CODE
    echo $CODE
done
exit $CODE
