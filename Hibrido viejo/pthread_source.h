#ifndef PTHREAD_SOURCE_H
#define PTHREAD_SOURCE_H

#include <pthread.h>

#ifndef G
#define G 6.673e-11
#endif

typedef struct cuerpo cuerpo_t;
struct cuerpo{
	float masa;
	float px;
	float py;
	float pz;
	float vx;
	float vy;
	float vz;
	float r;
	float g;
	float b;
	int cuerpo;
};

// Variables globales
extern float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
extern float *matriz_fuerzaX_l, *matriz_fuerzaY_l, *matriz_fuerzaZ_l;
extern float *matriz_fuerzaX_v, *matriz_fuerzaY_v, *matriz_fuerzaZ_v;
extern cuerpo_t *cuerpos;
extern int N, T, delta_tiempo, pasos, proceso;
extern pthread_mutex_t mutex1, mutex2;
extern pthread_barrier_t barrera;

// Declaraciones de funciones
void* gravitacion(void *arg);
void crear_hilos(int N_p, int T_p, int delta_tiempo_p, int pasos_p, int proceso_p,pthread_t* hilo_p);
void calcularFuerzas(int id);
void moverCuerpos(int id);
void cerrar_hilos();

#endif // PTHREAD_SOURCE_H