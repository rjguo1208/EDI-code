include ../../make.inc

MODFLAGS= $(BASEMOD_FLAGS) \
    	  $(MOD_FLAG)../../PW/src \
    	  $(MOD_FLAG)../../external/wannier90/src/obj \

PWOBJS = ../../PW/src/libpw.a
W90LIB = ../../external/wannier90/libwannier.a
QEMODS = ../../Modules/libqemod.a\
    	 ../../upflib/libupf.a\
    	 ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/src/libqefft.a\
    	 ../../UtilXlib/libutil.a\
    	 ../../XClib/xc_lib.a\
    	 ../../LAXlib/libqela.a\
    	 ../../MBD/libmbd.a\
    	 ../../dft-d3/libdftd3qe.a\

MODULES = $(PWOBJS) $(QEMODS) $(W90LIB)

F90FLAGS_W90 = $(F90FLAGS) -D__WANLIB

TLDEPS= pwlibs

# EDI core modules
EDIMODS = edic_mod.f90 delta_weights.f90 bz_symmetry.f90 \
          wigner_seitz.f90 wan2bloch.f90 wannierize.f90 edi_read_hr.f90

# Wannier/transport modules
W90MODS = ep_constants.f90 io_var_edi.f90 input_edi.f90 wann_common.f90 \
          global_var.f90 low_lvl_edi.f90 parallelism_edi.f90 \
          openfilepw_edi.f90 io_edi.f90 kfold_edi.f90 edi_input.f90 \
          edbloch2wan.f90 ed_coarse.f90 transport.f90

# Subroutines (matrix elements, potentials)
SUBS = get_betavkb.f90 get_vloc_onthefly.f90

OBJ = $(EDIMODS:.f90=.o) $(W90MODS:.f90=.o) $(SUBS:.f90=.o)

EXE = edi.x

.DEFAULT_GOAL := $(EXE)

pw2wan_edi.o: pw2wan_edi.f90
	$(MPIF90) $(F90FLAGS_W90) -c pw2wan_edi.f90 -o pw2wan_edi.o

edi_pw2wan.o: edi_pw2wan.f90
	$(MPIF90) $(F90FLAGS_W90) -c edi_pw2wan.f90 -o edi_pw2wan.o

edi.o: edi.f90
	$(MPIF90) $(F90FLAGS_W90) -c edi.f90 -o edi.o

$(EXE): $(OBJ) pw2wan_edi.o edi_pw2wan.o edi.o
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) pw2wan_edi.o edi_pw2wan.o edi.o $(MODULES) $(LIBOBJS) $(QELIBS)

clean:
	rm -f *.o *~ *.d *.mod *.i $(EXE)
