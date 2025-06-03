#!/bin/bash

# Comando 1
echo "Compilando pthread..."
gcc -lpthread -c pthread_source.c -o pthread_source.o

# Comando 2
echo "Compilando mpi..."
mpicc -pthread -c MPI_pthread.c -o MPI_pthread.o

# Comando 3
echo "Creando salida..."
mpicc -o salidaMPI_P MPI_pthread.o pthread_source.o -lm -pthread