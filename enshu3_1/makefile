CC=g++
CFLAGS=-Wall -O2

simulator: simulator.cpp
	$(CC) simulator.cpp -o simulator $(CFLAGS)

clean:
	rm -f simulator result result.gif

test: simulator
	./simulator
	python visualize.py
