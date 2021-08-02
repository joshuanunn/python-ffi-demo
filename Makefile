# C compiler
CC = gcc

# compiler flags:
CFLAGS = -O3

# build target
TARGET = disperse

all: $(TARGET)

$(TARGET): ./ffi_demo_ctypes/src/$(TARGET).c
	$(CC) $(CFLAGS) -o ./ffi_demo_ctypes/src/$(TARGET).o -c ./ffi_demo_ctypes/src/$(TARGET).c
	$(CC) -shared -o ./ffi_demo_ctypes/core/$(TARGET).so ./ffi_demo_ctypes/src/$(TARGET).o

test/ctypes:
	python3 -m ffi_demo_ctypes.tests.runtests
