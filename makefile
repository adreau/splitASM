SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o )

all: molsplit

%.o: %.cpp 
	g++ -g -fsanitize=address -std=c++17 -O3 -c $< -o $@

molsplit: $(OBJS)
	g++ -g -fsanitize=address -std=c++17 -o splitASM $(OBJS)

clean:
	rm -f *~ *.o
