# The comparator
# todo(s.schulz): use -Wall -Wextra in compiler dependant way
add_executable(TestCompare TestCompare.c)
target_include_directories(TestCompare PRIVATE ${STB_INCLUDE_DIR})
target_link_libraries(TestCompare LINK_PUBLIC m)
target_compile_options(TestCompare PRIVATE -Wall -Wextra -g)
# todo(s.schulz): should preferably be in a temporary directory that is
# cleaned after successful run, but mktemp functionality does not exist
# in CMake..
set(OUTPUTDIR .)
set(KNOWNGOOD ${CMAKE_CURRENT_SOURCE_DIR}/knowngood)
foreach( FEATURE
		ALL_Staggered
		ALL_Staggered_f
		ALL_Tensor
		ALL_Tensor_f
		)
	foreach( NP
			4 16 21 64 127 128
			)
		set(OUTPUTBINFILE ${OUTPUTDIR}/${FEATURE}_${NP}.bin)
		set(GOODBINFILE ${KNOWNGOOD}/${FEATURE}_${NP}.bin)
		set(OUTPUTFILE ${OUTPUTDIR}/out_${FEATURE}_${NP}.dat)
		set(GOODFILE ${KNOWNGOOD}/${FEATURE}_${NP}.dat)
		add_test(NAME feature_${FEATURE}_${NP}
			COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NP} ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${FEATURE}>)
        # tell CMake how many cores a test uses and allow openmpi to oversub cores via env variable
        set_tests_properties(feature_${FEATURE}_${NP} PROPERTIES PROCESSORS "${NP}" ENVIRONMENT "PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe")
		#add_test( NAME feature_${FEATURE}_${NP}_compare COMMAND $<TARGET_FILE:TestCompare> ${OUTPUTFILE} ${GOODFILE})
		add_test( NAME feature_${FEATURE}_${NP}_compare COMMAND ${CMAKE_COMMAND} -E compare_files ${OUTPUTBINFILE} ${GOODBINFILE})
        set_tests_properties(feature_${FEATURE}_${NP}_compare PROPERTIES DEPENDS feature_${FEATURE}_${NP})
	endforeach(NP)
endforeach(FEATURE)

# todo(s.schulz): generate known good results from the same set of
# parameters as above. For now run this for yourself for each NP using
# generate.sh in the knowngood directory. Make sure to run this without
# debug info, or additional output may be generated and saved

# vim: tw=72 et sw=4 ts=4
