tammes :	dSFMT.c dSFMT.h dSFMT-common.h dSFMT-params.h dSFMT-params19937.h\
			vec.h\
			avlmini.c avlmini.h\
			main.c
	gcc -Wall -o tammes -Ofast -flto -pipe -s -static -fopenmp -march=native -mtune=native -fno-strict-aliasing -DDSFMT_MEXP=19937 -DHAVE_SSE2 dSFMT.c avlmini.c main.c

clean :
	rm tammes.exe