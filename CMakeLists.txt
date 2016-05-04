cmake_minimum_required(VERSION 2.8.3)
project(soth_test)
find_package(cmake_modules REQUIRED)
find_package(PkgConfig REQUIRED)


find_package(catkin REQUIRED COMPONENTS expressiongraph)
pkg_search_module(SOTH REQUIRED soth)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  
if (SOTH_FOUND)

  message("${Red}include dir ${SOTH_INCLUDE_DIRS} ${ColourReset}")
 message("${Red}Lib ${SOTH_LIBRARIES} ${ColourReset}")
 message("${Red}Lib ${SOTH_VERSION} ${ColourReset}")
SET(SLIB ${SOTH_LIBRARIES}_${SOTH_VERSION})
 message("${Red}SLIB ${SLIB} ${ColourReset}")
else(SOTH_FOUND)
message(WARNING "not Found!!!!!!!!!!!!!!!!!!")
endif (SOTH_FOUND)

#target_link_libraries(testapp ${SOTH_LIBRARIES})
#target_include_directories(testapp PUBLIC ${SOTH_INCLUDE_DIRS})
#target_compile_options(testapp PUBLIC ${SOTH_CFLAGS_OTHER})




include_directories( include  ${SOTH_INCLUDE_DIRS} ${catkin_INCLUDE_DIRS})
add_library(rand src/RandomGenerator.cpp)
target_link_libraries(rand ${SOTH_LIBRARIES})
add_executable(test1 src/test1.cpp)
target_link_libraries(test1 rand ${SOTH_LIBRARIES})
add_executable(test2 src/test2.cpp)
target_link_libraries(test2 rand ${catkin_LIBRARIES} ${SOTH_LIBRARIES})
