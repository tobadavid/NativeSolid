# Compilation of the source code
# Generation of a shared librairy libNative.so
# Generation of an executable NativeSolid (linked with libNative.so)

#FILE(GLOB SRCS ../include/*.hpp *.cpp ../include/*.inl ../inclue/*.h)

ADD_LIBRARY(Native SHARED ../include/config.h ../include/geometry.h ../include/integration.h ../include/MatVec.h ../include/output.h ../include/structure.h ../include/solver.h config.cpp geometry.cpp integration.cpp MatVec.cpp output.cpp structure.cpp solver.cpp ../api/NativeSolid_API.h ../api/NativeSolid_API.cpp)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/api
)

TARGET_LINK_LIBRARIES(Native lapacke cblas)

#ADD_EXECUTABLE(NativeSolid ../include/Native.h ../api/main.cpp)
ADD_EXECUTABLE(TestCVector ../include/MatVec.h MatVec.cpp testcvector.cpp)
TARGET_LINK_LIBRARIES(TestCVector lapacke cblas)
