include ../make.inc

FFLAGS += -I..

ifeq (,$(BUILD_F90))
$(error Testing routine requires F95 interface be built. Set BUILD_F90 = 1 in make.inc and rebuild)
endif

ifeq (,$(FORT))
FORT = $(F90)
endif

OBJECTS = ../libpfapack.a \
	  error.o matrix_tools.o matrix_tools_s.o matrix_tools_d.o	\
	  matrix_tools_c.o matrix_tools_z.o check_decomp.o		\
	  check_decomp_s.o check_decomp_d.o check_decomp_c.o		\
	  check_decomp_z.o check_pfaffian.o check_pfaffian_s.o		\
	  check_pfaffian_d.o check_pfaffian_c.o check_pfaffian_z.o	\
	  test_ilaenv.o

test: test_pfaff.out test_skew.out
	./test_skew.out
	./test_pfaff.out

test_pfaff.out: test_pfaff.f90 $(OBJECTS) $(MODULES) ../libpfapack.a
	$(FORT) $(FFLAGS) test_pfaff.f90 -o test_pfaff.out $(OBJECTS) $(LFLAGS)

test_skew.out: test_skew.f90 $(OBJECTS) $(MODULES) ../libpfapack.a
	$(FORT) $(FFLAGS) test_skew.f90 -o test_skew.out $(OBJECTS) $(LFLAGS)

clean:
	rm -f *.out
	rm -f *.o
	rm -f *.mod

.SUFFIXES: .o .f .f90

.f90.o:
	$(FORT) $(FFLAGS) -c $< -o $@

.f.o:
	$(FORT) $(FFLAGS) -c $< -o $@
