#ifndef PTHREAD_SOURCE_H
#define PTHREAD_SOURCE_H

#include <pthread.h>

#ifndef G
#define G 6.673e-11
#endif

typedef struct {
    float px, py, pz;
    float vx, vy, vz;
    float masa;
} cuerpo_t;

// Variables globales
extern float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
extern float *matriz_fuerzaX_l, *matriz_fuerzaY_l, *matriz_fuerzaZ_l;
extern float *matriz_fuerzaX_v, *matriz_fuerzaY_v, *matriz_fuerzaZ_v;
extern cuerpo_t *cuerpos;
extern int N, T, delta_tiempo, pasos, proceso;
extern pthread_mutex_t mutex1, mutex2;
extern pthread_barrier_t barrera;

// Declaraciones de funciones
void gravitacion(int *arg);
void crear_hilos(int N_p, int T_p, int delta_tiempo_p, int pasos_p, int proceso_p,
                 cuerpo_t *cuerpos_p,
                 float *matriz_fuerzaX_l_p, float *matriz_fuerzaY_l_p, float *matriz_fuerzaZ_l_p,
                 float *matriz_fuerzaX_v_p, float *matriz_fuerzaY_v_p, float *matriz_fuerzaZ_v_p,
                 pthread_mutex_t mutex1_p, pthread_mutex_t mutex2_p, pthread_barrier_t barrera_p);
void calcularFuerzas(int id);
void moverCuerpos(int id);

#endif // PTHREAD_SOURCE_H
