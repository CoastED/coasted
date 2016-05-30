! These are hard-coded array / matrix sizes for different types
! of elements and quadratures

! Included to stop preprocessor complaining about redefines
#undef NFACES
#undef NDIM
#undef NLOC
#undef FLOC
#undef P_FLOC
#undef NGI
#undef FNGI

#undef SCHEME_CDG
#undef SCHEME_IP

! Defines for element types

#ifdef DG1_TET_QUAD_3_CDG
#define NFACES 4
#define NDIM 3
#define NLOC 4
#define FLOC 3
#define P_FLOC 6

#define NGI 5
#define FNGI 6

#define SCHEME_CDG

#endif

#ifdef DG1_TET_QUAD_3_IP
#define NFACES 4
#define NDIM 3
#define NLOC 4
#define FLOC 3
#define P_FLOC 6

#define NGI 5
#define FNGI 6

#define SCHEME_IP

#endif

#ifdef DG1_TET_QUAD_3_BASSI
#define NFACES 4
#define NDIM 3
#define NLOC 4
#define FLOC 3
#define P_FLOC 6

#define NGI 5
#define FNGI 6

#define SCHEME_BASSI

#endif


#ifdef DG1_TET_QUAD_4_CDG
#define NFACES 4
#define NDIM 3
#define NLOC 4
#define FLOC 3
#define P_FLOC 6

#define NGI 11
#define FNGI 6

#define SCHEME_CDG

#endif

#ifdef DG1_TET_QUAD_4_IP
#define NFACES 4
#define NDIM 3
#define NLOC 4
#define FLOC 3
#define P_FLOC 6

#define NGI 11
#define FNGI 6

#define SCHEME_IP

#endif

#ifdef DG1_TET_QUAD_4_BASSI
#define NFACES 4
#define NDIM 3
#define NLOC 4
#define FLOC 3
#define P_FLOC 6

#define NGI 11
#define FNGI 6

#define SCHEME_BASSI

#endif


#define EFLOC ( NLOC + NFACES * FLOC )
#define NBOURS NFACES
