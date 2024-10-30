# Simple makefile

## -fext-numeric-literals
## Needed for float suffix 'q' (float128)
CC := g++ -std=c++23 -O2 -fext-numeric-literals

#CCW := -Wall -Werror -Wextra

.DEFAULT_GOAL := main

main:
	$(CC) $(CCW) -o ./bin/main ./src/main.cpp
	./bin/main
	
all: clean run

clean:
	rm -rf ./bin/main

run: main.cpp
	./bin/main
