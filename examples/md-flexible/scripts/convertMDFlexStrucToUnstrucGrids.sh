#!/bin/bash

# Script to convert the old StructuredGrid vts files for the Rank information to UnstructuredGrid vtu files.
# This is useful because ParaView doesn't render inner surfaces in StructuredGrids.
# For more Information see: https://github.com/AutoPas/AutoPas/pull/955

# Check if the filename is provided
if [ -z "$1" ]; then
    echo "Usage: $(basename $0) <filename.vts>"
    echo
    echo "The scripts prints the new file to stdout. To create a new file pipe the output:"
    echo "  $(basename $0) foo.vts > foo_Ranks.vtu"
    exit 1
fi

InputFile=$1

# Perform the replacements and modifications
awk '
    # Replace StructuredGrid with UnstructuredGrid (and trailing spaces)
    {
        gsub(/StructuredGrid */, "UnstructuredGrid");
    }

    # Remove WholeExtent=".*"
    {
        gsub(/WholeExtent="[^"]*"/, "");
    }

    # Replace Extent=".*" with NumberOfPoints="8" NumberOfCells="1"
    {
        gsub(/Extent="[^"]*"/, "NumberOfPoints=\"8\" NumberOfCells=\"1\"");
    }

    # Print any other line
    {
        print;
    }

    # Insert the Cells block after </Points> line
    /<\/Points>/ {
        print "      <Cells>";
        print "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
        print "          0 1 2 3 4 5 6 7";
        print "        </DataArray>";
        print "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
        print "          8";
        print "        </DataArray>";
        print "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
        print "          11";
        print "        </DataArray>";
        print "      </Cells>";
    }
' "$InputFile"

