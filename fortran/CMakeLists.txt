CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
ENABLE_LANGUAGE(Fortran)

SET(SOURCES_PFAPACK
	dlasktrd.f dskmv.f dskr2.f dsktd2.f dsktrd.f
	zlasktrd.f zskmv.f zskr2.f zsktd2.f zsktrd.f
	dlasktrf.f dskpfa.f dsktf2.f dsktrf.f
	zlasktrf.f zskpfa.f zsktf2.f zsktrf.f
	)
if(USE_GEMMT)
    SET(SOURCES_PFAPACK
	    ${SOURCES_PFAPACK}
	    blas_dskr2k.f
	    blas_zskr2k.f
	    )
else(USE_GEMMT)
    SET(SOURCES_PFAPACK
	    ${SOURCES_PFAPACK}
	    dskr2k.f
	    zskr2k.f
	    )
endif(USE_GEMMT)
ADD_LIBRARY(pfapack STATIC ${SOURCES_PFAPACK})
TARGET_LINK_LIBRARIES(pfapack ${LAPACK_LIBRARIES})
INSTALL(TARGETS pfapack
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)

