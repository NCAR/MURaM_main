include ../Make_defs
############################################

MAKEFILE      = Makefile
SHELL         = /bin/sh

SUBDIRS	      = rt
LIBS	     += rt/librt.a

OBJS = ACCH.o \
       io_xysl.o \
       slice_write.o \
       slice_write_rebin.o \
       Iout.o \
       Qthin.o\
       tau_slice.o \
       xy_slice.o \
       yz_slice.o \
       xz_slice.o \
       corona_emission_dem_xyz.o \
       comm_split.o \
       driver.o \
       grid.o \
       init.o \
       run.o \
       physics.o \
       analysis_light.o \
       analysis_hmean_Corona.o \
       dfparser.o \
       Adjust_Valf_Max.o \
       Add_Sources_Integrate_SR.o \
       mhdres_SR.o \
       tvdlimit_SR.o \
       clean_divB.o \
       exchange.o \
       potential_sd_fftw3.o \
       potential_sd_fftw3_par.o \
       potential_sd_heffte_kernels.o \
       boundary_pdmp_1_fftw3.o \
       solver.o \
       eos_init.o \
       cons_to_prim.o 

all:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) all -C   $$subdir ; done
	$(MAKE) $(OBJS) $(HEADERS)
	@echo "Linking $(PROGRAM) ..."
	$(LD) -o $(PROGRAM) $(OBJS) $(LIBS) $(LDFLAGS)
	@echo "done"

clean:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) clean -C   $$subdir ; done
	rm -f $(OBJS) core

distclean: 
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) distclean -C   $$subdir ; done
	rm -f $(OBJS) core $(PROGRAM) 

install:
	cp -v $(PROGRAM) /usr/local/bin/

dep:
	@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) dep -C   $$subdir ; done
	$(MAKE) makedepend

makedepend:
	@makedepend -f$(MAKEFILE) $(OBJS:.o=.cc) $(INCLUDE)

update:
	$(FILEPATH)/$(PROGRAM) 

# DO NOT DELETE THIS LINE -- make depend depends on it.
