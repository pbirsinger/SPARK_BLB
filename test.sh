#!/bin/bash

PYTHONPATH=../../:$PYTHONPATH

echo PYTHONPATH
echo ${PYTHONPATH}

if [ -z "${PYTHON}" ]
then
    PYTHON=python
fi
if [ -z "${PYTHONARGS}" ]
then
    PYTHONARGS=
fi

PYTHONPATH=`pwd`:${PYTHONPATH} CLASSPATH=${CLASSPATH} ${PYTHON} ${PYTHONARGS} tests/blb_test.py $1
