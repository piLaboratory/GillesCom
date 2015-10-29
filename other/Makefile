CC=g++
CFLAGS=-I/usr/share/R/include -DNDEBUG   -I"/usr/local/lib/R/site-library/Rcpp/include" -I"/usr/local/lib/R/site-library/RcppArmadillo/include"   -fpic  -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -g -pg
LDFLAGS= -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -L/usr/lib/R/lib -pg
all: ../src/main.cpp profile.cpp
	$(CC) $(CFLAGS) -c profile.cpp -o profile.o
	$(CC)  $(LDFLAGS) profile.o -o profile -lR -lblas
