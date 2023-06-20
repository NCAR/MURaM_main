include ../../Make_defs
############################################

MAKEFILE      = Makefile

DEST	      = .
PROGRAM       = mhd3d

OBJS          = rt.o mem.o

SUBDIRS       = 
                
LIBRARY       = librt.a

all:            
		@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) all -C   $$subdir; done
		$(MAKE) $(OBJS) $(HEADERS)
		@echo "Creating $(LIBRARY) ..."
		@$(AR) cru $(LIBRARY) $(OBJS)
		@echo -n "$(AR) cru $(LIBRARY) $(OBJS) "
		@list='$(LIB_LIBS)' ; for lib in $$list ; do\
                  files=`$(AR) t $$lib | grep -v SYMDEF`;\
                  echo -n $$files" ";\
                  $(AR) x $$lib;\
                  $(AR) ru $(LIBRARY) $$files;\
                  $(RM) -f $$files;\
                done
		@echo ""
		$(RANLIB) $(LIBRARY)
		@echo "done"

clean:
		@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) clean -C   $$subdir; done
		$(RM) -f $(OBJS) core $(LIBRARY)

distclean:
		@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) distclean -C   $$subdir; done
		$(RM) -f $(OBJS) core $(LIBRARY) 

dep:		depend

depend:
		@list='$(SUBDIRS)' ; for subdir in $$list ; do $(MAKE) depend -C   $$subdir; done
		@makedepend -f $(MAKEFILE) $(SRCS)

update:		$(DEST)/$(LIBRARY) 

# DO NOT DELETE THIS LINE -- make depend depends on it.
