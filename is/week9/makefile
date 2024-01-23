CC=g++
CFLAGS=-Wall -O2
TARGET=training

all: $(TARGET)

.PHONY: clean all

$(TARGET): $(TARGET).cpp
	$(CC) $(TARGET).cpp -o $(TARGET) $(CFLAGS)

clean:
	rm -rf $(TARGET)
