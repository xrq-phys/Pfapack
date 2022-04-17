CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
ENABLE_LANGUAGE(Fortran)

SET(SOURCES_PFAPACK
	slasktrd.f sskmv.f sskr2.f ssktd2.f ssktrd.f
	dlasktrd.f dskmv.f dskr2.f dsktd2.f dsktrd.f
	clasktrd.f cskmv.f cskr2.f csktd2.f csktrd.f
	zlasktrd.f zskmv.f zskr2.f zsktd2.f zsktrd.f
	slasktrf.f sskpfa.f ssktf2.f ssktrf.f sskr2k.f
	dlasktrf.f dskpfa.f dsktf2.f dsktrf.f dskr2k.f
	clasktrf.f cskpfa.f csktf2.f csktrf.f cskr2k.f
	zlasktrf.f zskpfa.f zsktf2.f zsktrf.f zskr2k.f
	sskpf10.f sskbpf10.f sskbpfa.f sskbtrd.f
	dskpf10.f dskbpf10.f dskbpfa.f dskbtrd.f
	cskpf10.f cskbpf10.f cskbpfa.f cskbtrd.f
	zskpf10.f zskbpf10.f zskbpfa.f zskbtrd.f
	mul10.f
	)
ADD_LIBRARY(pfapack STATIC ${SOURCES_PFAPACK})
TARGET_LINK_LIBRARIES(pfapack ${LAPACK_LIBRARIES})
INSTALL(TARGETS pfapack
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)

