all: marble

marble: marble.c
	$(CC) -std=c99 -Wall -Wextra -pedantic-errors -g -o marble marble.c -lSDL2 -lm
