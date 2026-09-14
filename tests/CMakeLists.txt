add_library(IntraCrab ${CRAB_LIBS_TYPE} crab.cpp)
target_link_libraries(IntraCrab ${CRAB_LIBS})
# add_library(InterCrab ${CRAB_LIBS_TYPE} icrab.cpp)
# target_link_libraries(InterCrab ${CRAB_LIBS})
# add_library(BackwardCrab ${CRAB_LIBS_TYPE} bcrab.cpp)
# target_link_libraries(BackwardCrab ${CRAB_LIBS})

function (AddTest testSrc dirName)
  #Extract the filename without an extension (NAME_WE)
  get_filename_component(testName ${testSrc} NAME_WE)
  #Add compile target
  add_executable(${testName} ${testSrc} )
  #link with crab libraries
  if (${dirName} MATCHES "inter")
    target_link_libraries(${testName} ${Boost_PROGRAM_OPTIONS_LIBRARY} InterCrab)
  elseif (${dirName} MATCHES "backward")
    target_link_libraries(${testName} ${Boost_PROGRAM_OPTIONS_LIBRARY} BackwardCrab)
  else ()  
    target_link_libraries(${testName} ${Boost_PROGRAM_OPTIONS_LIBRARY} IntraCrab)
  endif ()
  #move testing binaries into a ${TEST_DIR} directory
  set_target_properties(${testName} PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY  ${TEST_DIR})
  #Finally add it to test execution
  add_test(NAME ${testName} 
    WORKING_DIRECTORY ${TEST_DIR}
    COMMAND ${TEST_DIR}/${testName} --disable-warnings)
endfunction()

set(EXCLUDE_FILES
  ""
 )

function (AddTestDir dirName)
  file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${dirName}/*.cc)
  foreach(testSrc ${TEST_SRCS})
    if (${testSrc} IN_LIST EXCLUDE_FILES)
      message (STATUS " ##### Excluding ${testSrc} from compilation")
    else()
      AddTest (${testSrc} ${dirName})      
    endif()
  endforeach(testSrc)
endfunction ()  

AddTestDir(assertion_crawler)
AddTestDir(domains)
AddTestDir(domains/wrapint)
AddTestDir(domains/region)
AddTestDir(domains/tvpi)
AddTestDir(transforms)
AddTestDir(inter)
AddTestDir(cg)
AddTestDir(cfg)
AddTestDir(thresholds)
AddTestDir(checkers)
AddTestDir(backward)
AddTestDir(preconditions)
AddTestDir(liveness)

# Boost.Test unit tests (self-checking). Their binaries live in
# ${CMAKE_BINARY_DIR}/unit-bin, not in ${TEST_DIR}, so that
# run_tests.sh (which diffs the stdout of every binary in ${TEST_DIR}
# against expected_results.out) never sees them.
set(UNIT_TEST_DIR ${CMAKE_BINARY_DIR}/unit-bin)
function (AddUnitTest testSrc)
  get_filename_component(testName ${testSrc} NAME_WE)
  add_executable(${testName} ${testSrc})
  target_link_libraries(${testName} Boost::unit_test_framework IntraCrab)
  target_compile_definitions(${testName} PRIVATE BOOST_TEST_DYN_LINK)
  set_target_properties(${testName} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${UNIT_TEST_DIR})
  add_test(NAME unit-${testName}
    WORKING_DIRECTORY ${UNIT_TEST_DIR}
    COMMAND ${UNIT_TEST_DIR}/${testName} --log_level=test_suite --report_level=short)
endfunction()

if (TARGET Boost::unit_test_framework)
  file(GLOB UNIT_TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} unit/*.cc)
  foreach(testSrc ${UNIT_TEST_SRCS})
    AddUnitTest(${testSrc})
  endforeach(testSrc)
else ()
  message (STATUS "Boost.Test (unit_test_framework) not found: skipping tests/unit")
endif ()
