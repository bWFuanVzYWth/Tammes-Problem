tammes :	dSFMT.c dSFMT.h dSFMT-common.h dSFMT-params.h dSFMT-params19937.h\
			vec.c vec.h\
			main.c
	gcc -o tammes -Ofast -flto -pipe -fopenmp -march=native -mtune=native -fno-strict-aliasing -DDSFMT_MEXP=19937 -DHAVE_SSE2 dSFMT.c vec.c main.c

clean :
	rm tammes*