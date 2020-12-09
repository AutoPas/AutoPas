#!/bin/bash
# This script adds a namespace to the specified file. The name of the namespace is the name of the file without the extension.
# It will add an opening `namespace namespacename {` after the last include.
# The namespace is closed by adding a line at the end of the file.

f=$1

namespace=$(basename $f)
namespace="${namespace%.*}"
line=`grep -n '^#include' $f | tail -1 | cut -f1 -d:`
head -$line $f > tempfile
echo "" >> tempfile
echo "namespace $namespace {" >> tempfile
let line++
tail --lines=+$line $f >> tempfile
echo "" >> tempfile
echo "}  // end namespace $namespace" >> tempfile
mv tempfile $f
