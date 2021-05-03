


hapsimul:	hash.o hapsimul.c
	$(CC) -fsanitize=address  -O0 -g -o hapsimul -std=c99 hapsimul.c  hash.o

hash.o: hash.c
	$(CC) -fsanitize=address -O0 -c -g -std=c99 hash.c


clean:
	/bin/rm -f *.o hapsimul


debug:
	 /opt/local/bin/gcc-mp-9  -O0 -g -o hapsimul -std=c99 hapsimul.c 
