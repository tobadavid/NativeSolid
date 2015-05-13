# Example of compiled problem using "waves"

#FILE(GLOB SRCS ../include/*.hpp *.cpp ../include/*.inl ../inclue/*.h)

ADD_LIBRARY(Native SHARED ../include/config.hpp ../include/integration.hpp ../include/MatVec.hpp ../include/Native.hpp ../include/output.hpp ../include/structure.hpp config.cpp integration.cpp MatVec.cpp Native.cpp output.cpp structure.cpp)

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include
)

TARGET_LINK_LIBRARIES(Native lapacke blas )

ADD_EXECUTABLE(NativeSolid ../include/Native.hpp main.cpp)

TARGET_LINK_LIBRARIES(NativeSolid Native lapacke blas )
