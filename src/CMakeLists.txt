# Compilation of the source code
# Generation of a shared librairy libNative.so

ADD_LIBRARY(Native SHARED ../include/config.h ../include/geometry.h ../include/integration.h ../include/MatVec.h ../include/output.h ../include/structure.h ../include/solver.h config.cpp geometry.cpp integration.cpp MatVec.cpp output.cpp structure.cpp solver.cpp ../api/NativeSolid_API.h ../api/NativeSolid_API.cpp)

TARGET_LINK_LIBRARIES(Native ${LAPACKE_LIBRARIES} ${CBLAS_LIBRARIES} ${BLAS_LIBRARIES})

ADD_EXECUTABLE(TestCVector ../include/MatVec.h MatVec.cpp testcvector.cpp)
TARGET_LINK_LIBRARIES(TestCVector ${LAPACKE_LIBRARIES} ${CBLAS_LIBRARIES} ${BLAS_LIBRARIES})
