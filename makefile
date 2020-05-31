# Makefile 
# Variables
CC = gcc
EJECUTABLE = main
# Programa Principal
all:
	$(CC) main.c -o $(EJECUTABLE) -lm
# Librerías

# Borrar los Archivos Objeto y el Ejecutable
clean:
	rm -rf *.o $(EJECUTABLE) 
