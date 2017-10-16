
SET(CBLAS_INCLUDE_SEARCH_PATHS
  /usr/include
  /usr/include/atlas
  /usr/include/openblas
  CACHE PATH ""
)

SET(CBLAS_LIB_SEARCH_PATHS
  /usr/lib
  /usr/lib/atlas-base
  /usr/lib/openblas-base
  CACHE PATH ""
)

SET(BLAS_LIB_SEARCH_PATHS
  /usr/lib
  /usr/lib/libblas
  /usr/lib/atlas-base
  /usr/lib/atlas-base/atlas
  /usr/lib/openblas-base
  CACHE PATH ""
)

SET(LAPACKE_INCLUDE_SEARCH_PATHS
  /usr/include
  CACHE PATH ""
)

SET(LAPACKE_LIB_SEARCH_PATHS
  /usr/lib
  CACHE PATH ""
)
