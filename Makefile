#Platform detection
UNAME_S     := $(shell uname -s)

#Compiler and Linker
# Using c++17 unlocks parallel algorithms in mesodyn when you:
#   Windows: MSVC uses its own builtin thread pool (no TBB needed)
#   Linux:   install TBB (apt install libtbb-dev)
#   MacOS:   install gcc and tbb via Homebrew (brew install gcc tbb)
#            Apple clang does not support the <execution> header.
ifeq ($(UNAME_S),Darwin)
	# MacOS specific settings
    BREW_PREFIX := $(shell brew --prefix)
    GCC_VERSION := $(shell ls $(BREW_PREFIX)/bin/g++-* | sort -t- -k2 -n | tail -1)
    CC          := $(if $(GCC_VERSION),$(GCC_VERSION),g++)
else
    CC          := g++
endif

NVCC        := $(shell which nvcc)
CUDA_DIR    := $(if $(NVCC),$(realpath $(dir $(NVCC))/..))
# c++ 17 is required if you'd like to use the parallel algorithms in mesodyn
CXX_STD     := c++14

#The Target Binary Program
TARGET      := namics

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := inc
BUILDDIR    := obj
TARGETDIR   := bin
RESDIR      := res
SRCEXT      := cpp
CUDAEXT		:= cu
DEPEXT      := d
OBJEXT      := o

#flat DOUBLE

#Flags, Libraries and Includes
CFLAGS      := -Wall -O3 -ffast-math -std=$(CXX_STD) -march=native
LIB         := -lm -lpthread
INC         := -I$(SRCDIR) -I/usr/local/include -I/usr/include -I/usr/include/eigen3

# MacOS: add Homebrew paths for headers and libraries
ifeq ($(UNAME_S),Darwin)
ifdef BREW_PREFIX
    INC     += -I$(BREW_PREFIX)/include
    LIB     += -L$(BREW_PREFIX)/lib
endif
endif

ifdef CUDA_DIR
INC         += -I$(CUDA_DIR)/include
endif

ifdef CUDA
	LIB        += -L$(CUDA_DIR)/lib64 -lcuda -lcudart -lcurand
	CFLAGS     += -DCUDA
	CUDA_ARCH  := native
	NVCCFLAGS  := -g -arch=$(CUDA_ARCH) -std=$(CXX_STD) -DCUDA -diag-suppress 20011,20012,20013,20014,20015,2809
	ifdef PAR_MESODYN_THRUST
		CFLAGS += -DPAR_MESODYN_THRUST
		NVCCFLAGS += --expt-relaxed-constexpr --expt-extended-lambda -DPAR_MESODYN_THRUST
	endif
else
    ifdef PAR_MESODYN_STL
		CFLAGS += -DPAR_MESODYN_STL
		LIB += -ltbb
    endif
endif

# %.o: %.cu $(NVCC) $(NVCCFLAGS) -c $< -o $@

#Build configuration summary
$(info )
$(info Platform:    $(UNAME_S))
$(info C++ standard: $(CXX_STD))
ifdef CUDA
$(info Compiler:    nvcc ($(NVCC)))
$(info CUDA dir:    $(CUDA_DIR))
$(info CUDA arch:   $(CUDA_ARCH))
ifdef PAR_MESODYN_THRUST
$(info Mesodyn:     Thrust (GPU parallel))
else
$(info Mesodyn:     serial)
endif
else
$(info Compiler:    $(CC))
ifdef PAR_MESODYN_STL
$(info Mesodyn:     C++17 parallel STL (Intel TBB))
else
$(info Mesodyn:     serial)
endif
endif
$(info )

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
CUDASOURCES := $(shell find $(SRCDIR) -type f -name *.$(CUDAEXT))
OBJECTS     := $(filter $(BUILDDIR)/%, $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT))) \
	$(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(CUDASOURCES:.$(CUDAEXT)=.$(OBJEXT))))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories
#	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
ifdef CUDA
# use nvcc
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(INC) -x cu -c -o $@ $<

$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(CUDAEXT)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $@ $<
else
# use regular cpp compiler
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INC) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp
endif


#Non-File Targets
.PHONY: all remake clean cleaner resources
