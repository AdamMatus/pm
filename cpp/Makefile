vpath %.cpp src/
vpath %.h inc/

OBJS= mel_frame_generator.o main.o dsp_utils.o vq.o 
CC=g++ -std=c++17
CPPFLAGS=-Wall -g3 -O0
INC=-Iinc
LDLIBS=-lgsl -lgslcblas -lm -lpython2.7 -lsndfile

default: test

test: $(OBJS)
	$(CC) $(OBJS) $(CPPFLAGS) $(LDLIBS)  -o $@

$(OBJS): %.o : %.cpp %.h
	$(CC) $(INC) $(CPPFLAGS)  -c $< -o $@

.PHONY: clean

clean:
	rm -rf *.o test

