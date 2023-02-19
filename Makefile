SRC = main.c vec.c random_point.c

debug: $(SRC)
	gcc -Wall -fexec-charset=GBK -o test $(SRC)

clear:
	rm *.exe