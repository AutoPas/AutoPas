#!/bin/sh

rel="../../../build/examples/md-julia/libjulia_bindings.so"
echo $rel

pwd
# use nm to create symbols
nm $rel > ./nm_output2.txt

# use objdump to create object dump
objdump -TC $rel > obj_output2.txt

# _ZN7autopas7AutoPasINS_10MoleculeLJIdEEEC1Ev
undef="_ZN7autopas7AutoPasINS_10MoleculeLJIdEEEC1Ev"
# undef="Mole"
# cat nm_output.txt | grep $undef