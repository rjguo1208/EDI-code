
default: all

all: edi

edi: pw
	cd src && make && make extract_pot.x

pw: 
	cd .. && make pw

clean:
	cd src && make clean
