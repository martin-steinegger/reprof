# This can be overridden e.g.: make install INSTALL_DIR=...
INSTALL_DIR?=$(PWD)

	
install:
	cd fannorg/src && make
	cd reprof && make
	mkdir -p $(INSTALL_DIR)/bin
	cp reprof/reprof $(INSTALL_DIR)/bin

clean:
	cd reprof && make clean
	cd fannorg/src && make clean
