find_package(Perl REQUIRED)

add_test(
  NAME scripts_helpers
  COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/helpers.t
          ${PROJECT_SOURCE_DIR}/scripts $<TARGET_FILE_DIR:integ>
)
set_tests_properties(scripts_helpers PROPERTIES TIMEOUT 120)

add_test(
  NAME comparator_contract
  COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/comparators.t
          ${PROJECT_SOURCE_DIR}/test
)
set_tests_properties(comparator_contract PROPERTIES TIMEOUT 120)

add_test(NAME chain_backend_contract
         COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/chain-backend.t ${PROJECT_SOURCE_DIR}/test $<TARGET_FILE:nrgchain>)
set_tests_properties(chain_backend_contract PROPERTIES TIMEOUT 120 ENVIRONMENT "${TEST_BUILD_ENV}")

add_test(NAME chain_registration_contract
         COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/chain-registration.t
                 ${PROJECT_SOURCE_DIR}/test ${CMAKE_COMMAND} ${CMAKE_CTEST_COMMAND})
set_tests_properties(chain_registration_contract PROPERTIES TIMEOUT 120)
