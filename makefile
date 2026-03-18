CC = gcc
CFLAGS = -O3 -march=native -mtune=native -fopenmp -pipe
LDFLAGS = -lm

TARGET = lantern_ga
SRC = lantern_ga.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
