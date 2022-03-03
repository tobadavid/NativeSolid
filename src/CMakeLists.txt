FILE(GLOB SRCS *.cpp)
ADD_LIBRARY(Native SHARED ${SRCS})

TARGET_LINK_LIBRARIES(Native ${LAPACKE_LIBRARIES} 
                             ${CBLAS_LIBRARIES} 
                             ${BLAS_LIBRARIES})
