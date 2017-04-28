vpath %.cpp src/
vpath %.h inc/

OBJS= main.o

CC=g++ -std=c++11
CPPFLAGS=-Wall -g3 -O0
LDLIBS=-lgsl -lgslcblas -lm -lpython2.7 -lsndfile

default: test

test: $(OBJS)
	$(CC) $(OBJS) $(CPPFLAGS) $(LDLIBS)  -o $@

$(OBJS): %.o : %.cpp
	$(CC) $(CPPFLAGS)  -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.o test

