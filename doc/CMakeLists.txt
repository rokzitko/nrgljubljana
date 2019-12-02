# Generate the sphinx config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# -----------------------------------------------------------------------------
# Create an optional target that allows us to regenerate the C++ doc with c++2rst
# -----------------------------------------------------------------------------
add_custom_target(docs_cpp2rst)
include(${PROJECT_SOURCE_DIR}/share/cmake/extract_flags.cmake)
extract_flags(nrgljubljana_c)
separate_arguments(nrgljubljana_c_CXXFLAGS)
macro(generate_docs header_file)
  add_custom_command(
    TARGET docs_cpp2rst
    COMMAND rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/cpp2rst_generated
    COMMAND
      c++2rst
      ${header_file}
      -N nrgljubljana
      --output_directory ${CMAKE_CURRENT_SOURCE_DIR}/cpp2rst_generated
      -I${PROJECT_SOURCE_DIR}/c++
      --cxxflags="${nrgljubljana_CXXFLAGS}"
  )
endmacro(generate_docs)

#generate_docs(${PROJECT_SOURCE_DIR}/c++/nrgljubljana/nrgljubljana.hpp)

# --------------------------------------------------------
# Build & Run the C++ doc examples and capture the output
# --------------------------------------------------------

add_custom_target(docs_example_output)
file(GLOB_RECURSE ExampleList RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach(example ${ExampleList})
  get_filename_component(f ${example} NAME_WE)
  get_filename_component(d ${example} DIRECTORY)
  add_executable(doc_${f} EXCLUDE_FROM_ALL ${example})
  set_property(TARGET doc_${f} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${d})
  add_custom_command(TARGET doc_${f}
   COMMAND doc_${f} > ${CMAKE_CURRENT_SOURCE_DIR}/${d}/${f}.output 2>/dev/null
   WORKING_DIRECTORY ${d}
  )
  add_dependencies(docs_example_output doc_${f})
endforeach()

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sphinx has internal caching, always run it
add_custom_target(docs_sphinx ALL)
add_custom_command(
  TARGET docs_sphinx
  COMMAND PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} sphinx-build -c . -j8 -b html ${CMAKE_CURRENT_SOURCE_DIR} html
)

option(Sphinx_Only "When building the documentation, skip the Python Modules and the generation of C++ Api and example outputs" OFF)
if(NOT Sphinx_Only)
  # Generation of C++ Api and Example Outputs
  add_dependencies(docs_sphinx docs_cpp2rst docs_example_output)
endif()

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/nrgljubljana
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|js|xsl|css|py|txt|inv|bib)$"
  PATTERN "_*"
)
