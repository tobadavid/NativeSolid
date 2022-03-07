# CMake input file of the SWIG wrapper around "Native.so"

INCLUDE(${SWIG_USE_FILE})

FILE(GLOB SRCS *.h *.cpp *.inl *.swg)
FILE(GLOB ISRCS *.i)

SET(CMAKE_SWIG_FLAGS "")
SET_SOURCE_FILES_PROPERTIES(${ISRCS} PROPERTIES CPLUSPLUS ON)

SET(SWINCFLAGS 
-I${PROJECT_SOURCE_DIR}/src 
)

SET_SOURCE_FILES_PROPERTIES(${ISRCS} PROPERTIES SWIG_FLAGS "${SWINCFLAGS}")

IF(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  SWIG_ADD_MODULE(NativeSolid python ${ISRCS} ${SRCS})
ELSE()
  SWIG_ADD_LIBRARY(NativeSolid
                   LANGUAGE python
                   SOURCES ${ISRCS} ${SRCS})
ENDIF()

TARGET_INCLUDE_DIRECTORIES(_NativeSolid PRIVATE ${PROJECT_SOURCE_DIR}/_src
                                          ${PROJECT_SOURCE_DIR}/src
                                          ${PYTHON_INCLUDE_PATH})

SWIG_LINK_LIBRARIES(NativeSolid Native ${PYTHON_LIBRARIES})
