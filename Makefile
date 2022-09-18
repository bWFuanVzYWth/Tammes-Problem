src = $(wildcard *.c)

release: $(src)
	gcc $(src) -o tammes -Wall -s -static -fexec-charset=GBK -Ofast -flto -pipe -march=native -fopt-info -fopenmp

debug: $(src)
	gcc $(src) -o test -Wall -fexec-charset=GBK -fopt-info -fopenmp -DDEBUG_TREESIZE -DDEBUG_TREE_CMP

clean:
	rm tammes.exe test.exe