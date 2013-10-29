#! /usr/bin/sh

testcases/tc_trivial/run_test.sh
return_code=$?

if [ ${return_code} -eq 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
fi

if [ ${return_code} -ne 0 ]; then
    echo "PASSED"
else
    echo "FAILED"
fi