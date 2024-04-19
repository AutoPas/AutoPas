#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <xsd-file>"
    exit 1
fi

# The `xsd` tool is used to generate a xml-parser for the given xsd file.
# The `xsd` tool is part of the CodeSynthesis XSD package. And can be downloaded from the following link:
# https://www.codesynthesis.com/download/xsd/4.2/

xsd cxx-tree --std c++17 --hxx-suffix .h --cxx-suffix .cpp --generate-doxygen $1