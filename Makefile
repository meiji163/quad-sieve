##
# quad-sieve
#

CC=g++
CFLAGS=-std=c++20 -I./include -O3
LIBS=-lgmp -lgmpxx

SRC := $(wildcard src/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRC))

build: $(OBJS)
	$(CC) $(SRC) -o quad-sieve $(CFLAGS) $(LIBS)

$(OBJS): %.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@ $(LIBS)

.PHONY: clean
clean:
	rm -f src/*.o quad-sieve

# end
