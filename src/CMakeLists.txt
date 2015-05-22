# Example of compiled problem using "waves"

#FILE(GLOB SRCS ../include/*.hpp *.cpp ../include/*.inl ../inclue/*.h)

ADD_LIBRARY(Native SHARED ../include/config.h ../include/integration.h ../include/MatVec.h ../include/Native.h ../include/output.h ../include/structure.h config.cpp integration.cpp MatVec.cpp Native.cpp output.cpp structure.cpp)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include
)

TARGET_LINK_LIBRARIES(Native lapacke blas )

ADD_EXECUTABLE(NativeSolid ../include/Native.h main.cpp)

TARGET_LINK_LIBRARIES(NativeSolid Native lapacke blas )
