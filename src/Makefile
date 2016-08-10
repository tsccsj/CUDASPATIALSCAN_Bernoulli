# CUDA
CUDA_PATH ?= "/sw/cuda/7.0"
CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDALIB   ?= -L$(CUDA_PATH)/lib64 -lcudart

GCC             ?= g++
NVCC		?= nvcc

# Gdal
GDALDIR = /sw/geosoft/gdal-1.11.2-fgdb
GDALLIB := -L$(GDALDIR)/lib -lgdal

# source files
TARGETS := io scan
OBJS    := $(TARGETS:=.o)
SRCS    := $(TARGETS:=.cu)
HDRS    := $(TARGETS:=.cuh)

# Common includes and paths for CUDA and GDAL
INCLUDES      := -I$(CUDA_INC_PATH) -I. -I$(GDALDIR)/include

# Target rules
all: cudaScan

io.o: io.c
	$(GCC) -I. -I$(GDALDIR)/include -o $@ -c $<
scan.o: scan.cu
	$(NVCC) -I$(CUDA_INC_PATH) -I. -o $@ -c $<
main.o: main.c
	$(GCC) -I. -o $@ -c $<
cudaScan: main.o io.o scan.o
	$(GCC) -o ../$@ $+ $(LDFLAGS) $(EXTRA_LDFLAGS) $(GDALLIB) $(CUDALIB)
clean:
	rm -f ../cudaScan $(OBJS) main.o
