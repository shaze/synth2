


hapsimul:	hash.o hapsimul.c
	$(CC)  -lbsd -O2 -g -o hapsimul -std=c99 hapsimul.c  hash.o

hash.o: hash.c
	$(CC)  -lbsd -O2 -c -g -std=c99 hash.c

ubuntu:
	$(CC) -fno-omit-frame-pointer -fsanitize=address  -O2 -g -o hapsimul -std=c99 hapsimul.c  hash.c -lbsd 


clean:
	/bin/rm -f *.o hapsimul


debug:
	 /opt/local/bin/gcc-mp-9  -O0 -g -o hapsimul -std=c99 hapsimul.c 
