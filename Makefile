INCLUDE = -I./include
LIBS=$(wildcard lib/*)
MACOS_VER=10.15


cflags += -ftree-vectorize -mmacosx-version-min=$(MACOS_VER) $(INCLUDE)
ldflags += -lm $(LIBS)

lib.name = psl

psl.class.sources := psl.c

datafiles = help-psl.pd


include Makefile.pdlibbuilder


render:
	@scripts/render.py && make

hash:
	@gcc -o ./scripts/hash scripts/hash-functions.c