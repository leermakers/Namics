#Compiler and Linker
CC          := g++
<<<<<<< HEAD
NVCC        :=nvcc
=======
>>>>>>> 6ca03aaec3f8a4a84b48ae8f4717a6e4b3acd76b

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
CFLAGS      := -Wall -g -O3 -std=c++11
LIB         := -lm -lpthread
INC         := -I/usr/local/include
#INCDEP      := -I$(INCDIR)
<<<<<<< HEAD
ifdef CUDA
	LIB        += -lcuda -lcudart
	NVCCFLAGS   := -DCUDA -arch sm_13
endif

# %.o: %.cu $(NVCC) $(NVCCFLAGS) -c $< -o $@
#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT)) 
ifdef CUDA
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.$(OBJEXT)) $(BUILDDIR)/tools.o)
else
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.$(OBJEXT)))
endif
=======

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))
>>>>>>> 6ca03aaec3f8a4a84b48ae8f4717a6e4b3acd76b

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

<<<<<<< HEAD
#Clean only Objects
=======
#Clean only Objecst
>>>>>>> 6ca03aaec3f8a4a84b48ae8f4717a6e4b3acd76b
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
<<<<<<< HEAD
#$(TARGET): $(OBJECTS)
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)


#Compile
ifdef CUDA
$(BUILDDIR)/tools.o: 
	$(NVCC) $(NVCCFLAGS) $(INC) -c -o $(BUILDDIR)/tools.o $(SRCDIR)/tools.cu
endif


=======
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
>>>>>>> 6ca03aaec3f8a4a84b48ae8f4717a6e4b3acd76b
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
