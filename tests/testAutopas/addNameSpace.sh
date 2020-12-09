#!/bin/bash
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
echo "} // end namespace $namespace" >> tempfile
mv tempfile $f
