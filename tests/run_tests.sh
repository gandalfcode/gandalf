#!/bin/bash

CODE=0
for file in $(cat to_run)
do
    ../bin/gandalf $file
    CODE=$? || $CODE
    echo $CODE
done
exit $CODE
