
INCLUDE = -I../include
LIBS=$(wildcard ../lib/*)
ldflags += -lm

all: test1 test2 test3 test4


test1: test1.c
	@gcc -o $@ -lm $(INCLUDE) $(LIBS) $@.c

test2: test2.c
	@gcc -o $@ -lm $(INCLUDE) $(LIBS) $@.c

test3: test3.c
	@gcc -o $@ -lm $(INCLUDE) $(LIBS) $@.c

test4: test4.c
	@gcc -o $@ -lm $(INCLUDE) $(LIBS) $@.c

clean:
	@rm -f test1 test2 test3 test4

