add_executable(testEarthModel
		testOscRndm.cc
		)
list(APPEND TEST_TARGETS testEarthModel)

add_executable(testOscRndm
		testOscRndm.cc
)
list(APPEND TEST_TARGETS testOscRndm)

foreach(test_tgt ${TEST_TARGETS})
	target_link_libraries(${test_tgt} PUBLIC OscLibFunc)
endforeach()


install(TARGETS ${TEST_TARGETS} RUNTIME DESTINATION ${TARGET_BINDIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
		GROUP_READ GROUP_EXECUTE
		WORLD_READ WORLD_EXECUTE)
