CC = g++ $(CFLAGS)
CFLAGS = -O0 -g

matrixTst: Matrix.h Matrix.cpp matrixTst.cpp
	$(CC) -o matrixTst matrixTst.cpp Matrix.cpp

clean:
	rm -f matrixTst
