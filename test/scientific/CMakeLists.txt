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

list(JOIN CHAIN_TEST_BACKENDS "," SCIENTIFIC_CHAIN_BACKENDS)
set(SCIENTIFIC_TEST_ENV "${TEST_BUILD_ENV}"
  "NRG_TEST_CHAIN_BACKENDS=${SCIENTIFIC_CHAIN_BACKENDS}"
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

if(TEST_SCIENTIFIC AND Mathematica_FOUND AND SYM_ALL)
  foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
    foreach(fixture IN ITEMS siam_qsz siam_qs siam_u0_qsz)
      set(test_name "scientific_prepare_${fixture}_${backend}")
      add_test(NAME ${test_name}
        COMMAND "${Python3_EXECUTABLE}" -B "${CMAKE_CURRENT_SOURCE_DIR}/prepare_siam.py"
          --source "${PROJECT_SOURCE_DIR}"
          --fixture "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/${fixture}"
          --kernel "${Mathematica_KERNEL_EXECUTABLE}"
          --nrg $<TARGET_FILE:nrg>
          --backend ${backend}
          --work-root "${CMAKE_CURRENT_BINARY_DIR}/runs/$<CONFIG>/${test_name}"
          --check)
      set_tests_properties(${test_name} PROPERTIES
        LABELS "scientific;chain-${backend};chain-generation"
        TIMEOUT 1200
        PROCESSORS 1
        ENVIRONMENT "${SCIENTIFIC_TEST_ENV}")
    endforeach()
  endforeach()
endif()

# Independent nonflat finite-star qualification. No production/default backend changes.
option(TEST_SCIENTIFIC_STAR_EXTENDED "Test both star unit conventions, temperatures and solver scales" OFF)
foreach(backend IN LISTS CHAIN_TEST_BACKENDS)
  set(star_frontends nrgchain runtime instantiate)
  if(TEST_SCIENTIFIC AND Mathematica_FOUND AND SYM_ALL)
    list(APPEND star_frontends nrginit)
  endif()
  foreach(frontend IN LISTS star_frontends)
    foreach(D IN ITEMS 2 1)
      foreach(temperature IN ITEMS 0.05 0.2)
        foreach(mode IN ITEMS rescaled absolute)
          set(star_tag "")
          if(NOT (D EQUAL 2 AND temperature STREQUAL "0.05" AND mode STREQUAL "rescaled"))
            if(NOT TEST_SCIENTIFIC_STAR_EXTENDED)
              continue()
            endif()
            string(REPLACE "." "" temperature_tag "${temperature}")
            set(star_tag "_D${D}_T${temperature_tag}_${mode}")
          endif()
          set(test_name "scientific_star_${frontend}${star_tag}_${backend}")
          set(star_args)
          set(star_labels "scientific;chain-${backend};chain-qualification")
          if(frontend STREQUAL "nrgchain")
            list(APPEND star_args --nrgchain $<TARGET_FILE:nrgchain>)
          elseif(frontend STREQUAL "instantiate")
            list(APPEND star_args --instantiate $<TARGET_FILE:instantiate>)
          elseif(frontend STREQUAL "nrginit")
            list(APPEND star_args --source "${PROJECT_SOURCE_DIR}" --kernel "${Mathematica_KERNEL_EXECUTABLE}")
            list(APPEND star_labels chain-generation)
          endif()
          if(NOT star_tag STREQUAL "")
            list(APPEND star_labels scientific-extended)
          endif()
          add_test(NAME ${test_name}
            COMMAND "${Python3_EXECUTABLE}" -B "${CMAKE_CURRENT_SOURCE_DIR}/validate_star.py"
              --fixture "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/asymmetric_star"
              --nrg $<TARGET_FILE:nrg> --backend ${backend} --lane ${frontend}
              --D ${D} --temperature ${temperature} --mode ${mode} ${star_args}
              --work-root "${CMAKE_CURRENT_BINARY_DIR}/runs/$<CONFIG>/${test_name}")
          set_tests_properties(${test_name} PROPERTIES
            LABELS "${star_labels}" TIMEOUT 900 PROCESSORS 1 ENVIRONMENT "${SCIENTIFIC_TEST_ENV}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()
