#include <math.h>
#include <stdlib.h>
#include "pthread_source.h"

float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
float *matriz_fuerzaX_l, *matriz_fuerzaY_l, *matriz_fuerzaZ_l;
float *matriz_fuerzaX_v, *matriz_fuerzaY_v, *matriz_fuerzaZ_v;
cuerpo_t *cuerpos;
int N, T, delta_tiempo, pasos, proceso;
pthread_mutex_t mutex1, mutex2;
pthread_barrier_t barrera;

void gravitacion(int id){
    int paso;
    for(paso=0; paso<pasos;paso++){
        //CALCULAR LAS FUERZAS Q LE TOCARON A LOS HILOS
        calcularFuerzas(id);

        pthreads_barrier_wait(&barrera);
        if ((id == 0) || (id == 1)){
            pthreads_mutex_unlock(&mutex1);//semaforo (aviso a mpi que termine)
            pthreads_mutex_lock(&mutex2);//semaforo (espera a mpi)
        }

        moverCuerpos(id);

        pthreads_barrier_wait(&barrera); //barrera
        if ((id == 0) || (id == 1)) {
            pthreads_mutex_unlock(&mutex1);//semaforo (aviso a mpi que termine)
            pthreads_mutex_lock(&mutex2);//semaforo (espera a mpi)
        }
    }
}

void crear_hilos(int N_p,int T_p,int delta_tiempo_p,int pasos_p,int proceso_p,cuerpo_t *cuerpos_p,float *matriz_fuerzaX_l_p,float *matriz_fuerzaY_l_p,float *matriz_fuerzaZ_l_p,float *matriz_fuerzaX_v_p,float *matriz_fuerzaY_v_p,float *matriz_fuerzaZ_v_p,pthread_mutex_t mutex1_p,pthread_mutex_t mutex2_p,pthread_barrier_t barrera_p){
    //INICIALIZAR VARIABLES RECIBIDAS POR PARAMETRO A GLOBALES
    int i;
    N=N_p;
    T=T_p;
    delta_tiempo=delta_tiempo_p;
    pasos=pasos_p;
    proceso=proceso_p;
    cuerpos=cuerpos_p;
    matriz_fuerzaX_l = matriz_fuerzaX_l_p;
    matriz_fuerzaY_l = matriz_fuerzaY_l_p;
    matriz_fuerzaZ_l = matriz_fuerzaZ_l_p;
    matriz_fuerzaX_v = matriz_fuerzaX_v_p;
    matriz_fuerzaY_v = matriz_fuerzaY_v_p;
    matriz_fuerzaZ_v = matriz_fuerzaZ_v_p;
    mutex1=mutex1_p;
    mutex2=mutex2_p;
    barrera=barrera_p;

    pthread_t hilo[T];
    for(int i = 0; id<T; i++){
        pthread_create(&hilo[i], NULL, gravitacion, ((void*)i*2)+proceso);
    }

	//
    for(int i = 0; i<T; i++){
        pthread_join(hilo[i], NULL);
    }
}

void calcularFuerzas(int id){
int cuerpo1, cuerpo2;
float dif_X, dif_Y, dif_Z;
float r;
float F;
int t2=T*2;

	for(cuerpo1 = id; cuerpo1<N-1 ; cuerpo1+=t2){
		for(cuerpo2 = cuerpo1 + 1; cuerpo2<N ; cuerpo2++){
			if ( (cuerpos[cuerpo1].px == cuerpos[cuerpo2].px) && (cuerpos[cuerpo1].py == cuerpos[cuerpo2].py) && (cuerpos[cuerpo1].pz == cuerpos[cuerpo2].pz))
                continue;

	            	dif_X = cuerpos[cuerpo2].px - cuerpos[cuerpo1].px; //X2-X1
			        dif_Y = cuerpos[cuerpo2].py - cuerpos[cuerpo1].py; //Y2-Y1
			        dif_Z = cuerpos[cuerpo2].pz - cuerpos[cuerpo1].pz; //Z2-Z3
                
			        r = sqrt(dif_X*dif_X + dif_Y*dif_Y + dif_Z*dif_Z);

	                F = (G*cuerpos[cuerpo1].masa*cuerpos[cuerpo2].masa)/(r*r);

	                dif_X *= F; 
			        dif_Y *= F;
			        dif_Z *= F;
			
					matriz_fuerzaX[id*N+cuerpo1] += dif_X;
	                matriz_fuerzaY[id*N+cuerpo1] += dif_Y;
	                matriz_fuerzaZ[id*N+cuerpo1] += dif_Z;

	                matriz_fuerzaX[id*N+cuerpo2] -= dif_X;
	                matriz_fuerzaY[id*N+cuerpo2] -= dif_Y;
	                matriz_fuerzaZ[id*N+cuerpo2] -= dif_Z;
	                
		}
	}
}

void moverCuerpos(int id){
 	int cuerpo,i,j;
	int t2=2*T;
	for(cuerpo = id; cuerpo<N ; cuerpo+=t2){

		for (i=0;i<T;i++){
			fuerza_totalX[i] += matriz_fuerzaX_l[i*N+cuerpo]+matriz_fuerzaX_v[i*N+cuerpo];
			fuerza_totalY[i] += matriz_fuerzaY_l[i*N+cuerpo]+matriz_fuerzaY_v[i*N+cuerpo];
			//fuerza_totalZ[i] += matriz_fuerzaZ[i*N+cuerpo];
			matriz_fuerzaX_l[i*N+cuerpo] = 0.0;
			matriz_fuerzaY_l[i*N+cuerpo] = 0.0;
            matriz_fuerzaX_v[i*N+cuerpo] = 0.0;
			matriz_fuerzaY_v[i*N+cuerpo] = 0.0;
			//matriz_fuerzaZ[i*N+cuerpo] = 0.0;
		}
		
        fuerza_totalX[cuerpo] *= 1/cuerpos[cuerpo].masa;
        fuerza_totalY[cuerpo] *= 1/cuerpos[cuerpo].masa;
        //fuerza_totalZ[cuerpo] *= 1/cuerpos[cuerpo].masa;

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo]*dt;
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo]*dt;
        //cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx *dt;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy *dt;
        //cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *dt;

		fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

	}
}