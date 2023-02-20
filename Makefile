SRC = main.c vec.c random_point.c tammes.c

debug: $(SRC)
	gcc -Wall -fexec-charset=GBK -o test -DDEBUG $(SRC)

clear:
	rm *.exe