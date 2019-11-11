ADD_EXECUTABLE(test_work_precision test_work_precision.F90)
TARGET_LINK_LIBRARIES(test_work_precision work_precision)
#ADD_TEST(NAME t_work_precision COMMAND test_work_precision)


#SET_TESTS_PROPERTIES(work_precision PROPERTIES PASS_REGULAR_EXPRESSION "PASSED")
