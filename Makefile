default:
	make -C src
clean: 
	make clean -C src
        
distclean: 
	make distclean -C src
	rm -rf lib/{lib,include,share}
        

