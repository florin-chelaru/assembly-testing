#! /usr/bin/sh

timestamp=`date +%FT%R:%S`

echo "<testsuites>"
echo "    <testsuite name=\"Assembler testing\" tests=\"2\" failures=\"1\" timestamp=\"${timestamp}\">"

for tc in testcases/* ; do
    tc=`basename ${tc}`
    testcases/${tc}/run_test.sh
    return_code=$?

    if [ ${return_code} -eq 0 ]; then
        echo "<testcase name=\"${tc}\" classname=\"None\">"
        echo "</testcase>"
    else
        echo "<testcase name=\"${tc} fail\" classname=\"None\">"
        echo "     <error message=\"Failed\">Failed2.</error>"
        echo "</testcase>"
    fi

    if [ ${return_code} -ne 0 ]; then
        echo "<testcase name=\"${tc}\" classname=\"None\">"
        echo "</testcase>"
    else
        echo "<testcase name=\"${tc} fail\" classname=\"None\">"
        echo "     <error message=\"Failed\">Failed2.</error>"
        echo "</testcase>"
    fi

done

echo "    </testsuite>"
echo "</testsuites>"

