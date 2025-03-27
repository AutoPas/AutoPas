for traversal in lc_sliced lc_sliced_balanced lc_sliced_c02 lc_c01 lc_c01_combined_SoA lc_c04 lc_c04_HCP lc_c04_combined_SoA lc_c08 lc_c18 ;
do
    for pattern in 1xVectorLength 2xVectorLengthDiv2 VectorLengthDiv2x2 VectorLengthx1 ;
    do

        ./md-flexible --yaml-filename Exploding_Cube_hwy.yaml --vectorizationPattern=$pattern --traversal=$traversal

    done
done