# C compiler
CC = gcc

# compiler flags:
CFLAGS = -O3

# build target
TARGET = disperse

all: $(TARGET)

$(TARGET): ./ffi_demo/ctypes/$(TARGET).c
	$(CC) $(CFLAGS) -o ./build/ctypes/$(TARGET).o -c ./ffi_demo/ctypes/$(TARGET).c
	$(CC) -shared -o ./build/ctypes/$(TARGET).so ./build/ctypes/$(TARGET).o

test/ctypes:
	python3 ./test/ctypes_test.py
