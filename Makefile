INCLUDE = -I./include
LIBS=$(wildcard lib/*)
MACOS_VER=10.15


cflags += -ftree-vectorize -mmacosx-version-min=$(MACOS_VER) $(INCLUDE)
ldflags += -lm $(LIBS)

lib.name = psl

class.sources := psl.c

datafiles = help-psl.pd

PDLIBBUILDER_DIR=pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder