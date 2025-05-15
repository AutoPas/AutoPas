for traversal in vcl_c06 vcl_sliced vcl_sliced_balanced vcl_sliced_c02;
do
    for pattern in 1xVectorLength 2xVectorLengthDiv2 VectorLengthDiv2x2 VectorLengthx1 ;
    do

        #for newton3 in enabled disabled;
        #do
            ./md-flexible --yaml-filename Gravity_Sphere.yaml --vectorizationPattern=$pattern --traversal=$traversal #--newton3=$newton3
        #done
    done
done