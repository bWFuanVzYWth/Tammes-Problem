src = $(wildcard *.c)

debug: $(src)
	gcc $(src) -o test -Wall -fexec-charset=GBK -fopt-info -fopenmp -DDEBUG_TREESIZE -DDEBUG_TREE_CMP

release: $(src)
	gcc $(src) -o tammes -Wall -s -static -fexec-charset=GBK -Ofast -flto -pipe -march=native -fopt-info -fopenmp

clean:
	rm tammes.exe test.exe