find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -B -c [=[
import sys
import numpy
from numpy.lib import NumpyVersion
print(numpy.__version__)
sys.exit(not (NumpyVersion("1.26.0") <= NumpyVersion(numpy.__version__) < NumpyVersion("3.0.0")))
]=]
  RESULT_VARIABLE SCIENTIFIC_NUMPY_RESULT
  OUTPUT_VARIABLE SCIENTIFIC_NUMPY_VERSION
  ERROR_VARIABLE SCIENTIFIC_NUMPY_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT SCIENTIFIC_NUMPY_RESULT STREQUAL "0")
  message(FATAL_ERROR
    "TEST_SCIENTIFIC=ON requires NumPy >=1.26,<3 for ${Python3_EXECUTABLE}.\n"
    "Install with: \"${Python3_EXECUTABLE}\" -m pip install -r \"${CMAKE_CURRENT_SOURCE_DIR}/requirements.txt\"\n"
    "NumPy check output: ${SCIENTIFIC_NUMPY_VERSION}\n${SCIENTIFIC_NUMPY_ERROR}")
endif()
message(STATUS "Scientific tests: NumPy ${SCIENTIFIC_NUMPY_VERSION}")

set(SCIENTIFIC_TEST_ENV "${TEST_BUILD_ENV}"
  "OPENBLAS_NUM_THREADS=1"
  "MKL_NUM_THREADS=1"
  "OMP_NUM_THREADS=1"
  "BLIS_NUM_THREADS=1"
  "VECLIB_MAXIMUM_THREADS=1"
  "NUMEXPR_NUM_THREADS=1"
  "OMP_DYNAMIC=FALSE"
  "MKL_DYNAMIC=FALSE"
  "LC_ALL=C")

add_test(NAME scientific_unit
  COMMAND "${Python3_EXECUTABLE}" -B -m unittest discover
    -s "${CMAKE_CURRENT_SOURCE_DIR}" -p "test_*.py")
set(scientific_tests scientific_unit)

foreach(fixture IN ITEMS siam_qsz siam_qs siam_u0_qsz)
  foreach(temperature IN ITEMS 0.05 0.2)
    string(REPLACE "." "" temperature_tag "${temperature}")
    foreach(mode IN ITEMS rescaled absolute)
      set(test_name "scientific_${fixture}_T${temperature_tag}_${mode}")
      add_test(NAME ${test_name}
        COMMAND "${Python3_EXECUTABLE}" -B "${CMAKE_CURRENT_SOURCE_DIR}/validate_siam.py"
          --nrg $<TARGET_FILE:nrg>
          --fixture "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/${fixture}"
          --temperature ${temperature}
          --mode ${mode}
          --work-root "${CMAKE_CURRENT_BINARY_DIR}/runs/$<CONFIG>/${test_name}")
      list(APPEND scientific_tests ${test_name})
    endforeach()
  endforeach()
endforeach()

set_tests_properties(${scientific_tests} PROPERTIES
  LABELS scientific
  TIMEOUT 180
  PROCESSORS 1
  ENVIRONMENT "${SCIENTIFIC_TEST_ENV}")
