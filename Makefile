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
CFLAGS      := -Wall -Ofast -g -std=c++14 -march=native
LIB         := -lm -lpthread
INC         := -I/usr/local/cuda-10.0/include -I/usr/local/include -I/usr/include
#INCDEP      := -I$(INCDIR)
ifdef CUDA
	LIB        += -L/usr/local/cuda-10.0/lib64 -lcuda -lcudart
	CFLAGS     += -DCUDA
	NVCCFLAGS  := -ccbin gcc-7 -arch=sm_61 -std=c++14 -DCUDA
	ifdef PAR_MESODYN
		CFLAGS += -DPAR_MESODYN
		NVCCFLAGS += --expt-relaxed-constexpr -DPAR_MESODYN
	endif
endif

# %.o: %.cu $(NVCC) $(NVCCFLAGS) -c $< -o $@
#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.$(OBJEXT)))
ifdef CUDA
OBJECTS     += $(BUILDDIR)/tools.o
ifdef PAR_MESODYN
OBJECTS     += $(BUILDDIR)/mesodyn.o $(BUILDDIR)/neighborlist.o $(BUILDDIR)/boundary_conditions.o
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
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $(BUILDDIR)/neighborlist.o $(SRCDIR)/mesodyn/neighborlist.cu
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $(BUILDDIR)/boundary_conditions.o $(SRCDIR)/mesodyn/boundary_conditions.cu

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