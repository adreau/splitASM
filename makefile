SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

all: molsplit

%.o: %.cpp 
	g++ -g -fsanitize=address -std=c++11 -O3 -Wall -c $< -o $@

molsplit: $(OBJS)
	g++ -g -fsanitize=address -std=c++11 -O3 -Wall -o splitASM $(OBJS)

clean:
	rm -f *~ *.o
