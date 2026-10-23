# Makefile for crispr_scan — builds the binary in place (./crispr_scan).
# Flags per spec: -O3 -march=native -pthread -Wall -Wextra -std=c11 (zero warnings allowed).

CC      = gcc
CFLAGS  = -O3 -march=native -pthread -Wall -Wextra -std=c11
TARGET  = crispr_scan

all: $(TARGET)

$(TARGET): crispr_scan.c
	$(CC) $(CFLAGS) -o $(TARGET) crispr_scan.c

test: all
	python3 ../tests/test_scanner.py

clean:
	rm -f $(TARGET)

.PHONY: all test clean