find_package(Perl REQUIRED)

add_test(
  NAME scripts_helpers
  COMMAND ${PERL_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/helpers.t
          ${PROJECT_SOURCE_DIR}/scripts $<TARGET_FILE_DIR:integ>
)
set_tests_properties(scripts_helpers PROPERTIES TIMEOUT 120)
