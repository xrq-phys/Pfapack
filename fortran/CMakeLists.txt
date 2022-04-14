CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
ENABLE_LANGUAGE(Fortran)

SET(SOURCES_PFAPACK
	dlasktrd.f dskmv.f dskr2.f dsktd2.f dsktrd.f
	zlasktrd.f zskmv.f zskr2.f zsktd2.f zsktrd.f
	dlasktrf.f dskpfa.f dsktf2.f dsktrf.f
	zlasktrf.f zskpfa.f zsktf2.f zsktrf.f
	dskr2k.f
	zskr2k.f
	)
ADD_LIBRARY(pfapack STATIC ${SOURCES_PFAPACK})
TARGET_LINK_LIBRARIES(pfapack ${LAPACK_LIBRARIES})
INSTALL(TARGETS pfapack
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)

