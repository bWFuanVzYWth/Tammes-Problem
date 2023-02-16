src = $(wildcard *.c)

release: $(src)
	gcc $(src) -o tammes -Wall -s -static -fexec-charset=GBK -Ofast -flto -pipe -mavx2 -fopt-info -fopenmp

release_sapphirerapids: $(src)
	gcc $(src) -o tammes_sapphirerapids -Wall -s -static -fexec-charset=GBK -Ofast -flto -pipe -march=sapphirerapids -mtune=sapphirerapids -fopt-info -fopenmp

debug: $(src)
	gcc $(src) -o test -Wall -fexec-charset=GBK -fopt-info -fopenmp -DDEBUG_TREESIZE -DDEBUG_TREE_CMP

clean:
	rm tammes.exe test.exe