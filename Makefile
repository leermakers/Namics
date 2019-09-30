#Compiler and Linker
ifdef GCC_8
CC          :=g++-8
else
CC			:=g++
endif

NVCC        :=nvcc

#The Target Binary Program
TARGET      := namics

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := inc
BUILDDIR    := obj
TARGETDIR   := bin
RESDIR      := res
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CFLAGS      := -Wall -Ofast -g -std=c++14 -fopenmp -march=native -ftree-parallelize-loops=12
LIB         := -lm -lpthread -lgomp
INC         := -I/usr/local/include -I/usr/include

ifdef CLENG_EXPERIMENTAL
CFLAGS      := -Wall -Ofast -std=c++14 -DCLENG_EXPERIMENTAL
LIB         += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5
INC         += -I/usr/include/hdf5/serial
endif

#INCDEP      := -I$(INCDIR)
ifdef CUDA
	LIB        += -lcuda -lcudart
	CFLAGS     += -DCUDA
	NVCCFLAGS  := -std=c++14 -DCUDA
	ifdef PAR_MESODYN
		CFLAGS += -DPAR_MESODYN
		NVCCFLAGS += -ccbin gcc-5 --expt-relaxed-constexpr -DPAR_MESODYN
	endif
endif

# %.o: %.cu $(NVCC) $(NVCCFLAGS) -c $< -o $@
#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
ifdef PAR_MESODYN
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.$(OBJEXT)) $(BUILDDIR)/tools.o $(BUILDDIR)/mesodyn.o)
else
ifdef CUDA
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.$(OBJEXT)) $(BUILDDIR)/tools.o)
else
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.$(OBJEXT)))
endif
endif

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
#$(TARGET): $(OBJECTS)
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)


#Compile
ifdef PAR_MESODYN
$(BUILDDIR)/tools.o:
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $(BUILDDIR)/tools.o $(SRCDIR)/tools.cu
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $(BUILDDIR)/mesodyn.o $(SRCDIR)/mesodyn.cu
else
ifdef CUDA
$(BUILDDIR)/tools.o:
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $(BUILDDIR)/tools.o $(SRCDIR)/tools.cu
endif
endif

$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all remake clean cleaner resources
