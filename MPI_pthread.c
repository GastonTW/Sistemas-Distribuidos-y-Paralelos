// Compilar:
//		gcc -o n_body_simple_NOGL n_body_simple.c -lm
// Ejecutar:
//		./n_body_simple_NOGL <nro de cuerpos> <DT> <Pasos>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "pthread_source.h"

//
// Para tiempo de ejecucion
//

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

double tIni, tFin, tTotal;

//
// Constantes para Algoritmo de gravitacion
//
#define PI (3.141592653589793)
#define ESTRELLA 0
#define POLVO 1
#define H2 2 //Hidrogeno molecular

// ===============
// ===== CPU =====
// ===============

//
// Estructuras y variables para Algoritmo de gravitacion
//
/*typedef struct cuerpo cuerpo_t;
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
};*/
void inicializarCuerpos(cuerpo_t *cuerpos, int N,
                        float *toroide_alfa, float *toroide_theta,
                        float *toroide_incremento, float *toroide_lado,
                        float *toroide_r, float *toroide_R);

void finalizar(MPI_Datatype MPI_CUERPO);

void crear_tipo_cuerpo(MPI_Datatype *tipo) {
    struct cuerpo temp;

    int bloques = 2;
    int tamanios[] = {10, 1}; // 10 floats + 1 int
    MPI_Aint desplazamientos[2];
    MPI_Datatype tipos[] = {MPI_FLOAT, MPI_INT};

    // Calcular desplazamientos
    MPI_Aint base;
    MPI_Get_address(&temp, &base);
    MPI_Get_address(&temp.masa, &desplazamientos[0]);
    MPI_Get_address(&temp.cuerpo, &desplazamientos[1]);
    desplazamientos[0] -= base;
    desplazamientos[1] -= base;

    // Crear el tipo derivado
    MPI_Type_create_struct(bloques, tamanios, desplazamientos, tipos, tipo);
    MPI_Type_commit(tipo);
}
//
//
void procesoA(int N,int dt,int pasos, int T,MPI_Datatype MPI_CUERPO,pthread_t* hilo){
	int i,j,c,paso;
	int nt=N*T;
	pthread_mutex_lock(&mutex1);
	pthread_mutex_lock(&mutex2);
	
	crear_hilos(N,dt,pasos,T,0,hilo);
	
	for(paso=0; paso<pasos;paso++){
		
		pthread_mutex_lock(&mutex1);
		
		MPI_Ssend(matriz_fuerzaX_l, nt, MPI_FLOAT, 1, paso, MPI_COMM_WORLD);
		MPI_Ssend(matriz_fuerzaY_l, nt, MPI_FLOAT, 1, paso, MPI_COMM_WORLD);
		MPI_Ssend(matriz_fuerzaZ_l, nt, MPI_FLOAT, 1, paso, MPI_COMM_WORLD);
        	MPI_Recv(matriz_fuerzaX_v, nt, MPI_FLOAT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaY_v, nt, MPI_FLOAT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaZ_v, nt, MPI_FLOAT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		pthread_mutex_unlock(&mutex2);
		pthread_mutex_lock(&mutex1); // espero a mover cuerpos
		MPI_Ssend(cuerpos, N, MPI_CUERPO, 1, paso, MPI_COMM_WORLD);
		MPI_Recv(cuerpos, N, MPI_CUERPO, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		pthread_mutex_unlock(&mutex2); //semaforo unlock
	}
	printf("cerrando hilos\n");
	cerrar_hilos();

	printf("hilos cerrados\n");
}

void procesoB(int N,int dt,int pasos, int T,MPI_Datatype MPI_CUERPO, pthread_t* hilo){
	int i,j,c,paso;
	int nt=N*T;

	pthread_mutex_lock(&mutex1);
	pthread_mutex_lock(&mutex2);
	crear_hilos(N,dt,pasos,T,1,hilo);
	for(paso=0; paso<pasos;paso++){
		
		pthread_mutex_lock(&mutex1);
		
		MPI_Recv(matriz_fuerzaX_v, nt, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaY_v, nt, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaZ_v, nt, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Ssend(matriz_fuerzaX_l, nt, MPI_FLOAT, 0, paso, MPI_COMM_WORLD);
		MPI_Ssend(matriz_fuerzaY_l, nt, MPI_FLOAT, 0, paso, MPI_COMM_WORLD);
		MPI_Ssend(matriz_fuerzaZ_l, nt, MPI_FLOAT, 0, paso, MPI_COMM_WORLD);
		pthread_mutex_unlock(&mutex2);
		pthread_mutex_lock(&mutex1); // espero a mover cuerpos
		MPI_Recv(cuerpos, N, MPI_CUERPO, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Ssend(cuerpos, N, MPI_CUERPO, 0, paso, MPI_COMM_WORLD);
		pthread_mutex_unlock(&mutex2); //semaforo unlock
	}
	cerrar_hilos();
}

int funcion_mpi(int rank,int N,int delta_tiempo,int pasos,int T,float *toroide_alfa,float *toroide_theta,
				float *toroide_incremento,float *toroide_lado,float *toroide_r,float *toroide_R){
	int i,j;
	double tIni, tFin;
	pthread_t hilo[T];

	cuerpos = (cuerpo_t*)malloc(sizeof(cuerpo_t)*N);
	fuerza_totalX = (float*)malloc(sizeof(float)*N);
	fuerza_totalY = (float*)malloc(sizeof(float)*N);
	fuerza_totalZ = (float*)malloc(sizeof(float)*N);
	matriz_fuerzaX_l = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaY_l = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaZ_l = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaX_v = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaY_v = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaZ_v = (float*)malloc(sizeof(float)*N*T);
	
	inicializarCuerpos(cuerpos,N,toroide_alfa,toroide_theta,toroide_incremento,toroide_lado,toroide_r,toroide_R); 
	
	MPI_Datatype MPI_CUERPO;
	crear_tipo_cuerpo(&MPI_CUERPO);
	
	for(i=0;i<N;i++){
		//inicializo matriz fuerza total
		fuerza_totalX[i] = 0.0;
		fuerza_totalY[i] = 0.0;
		fuerza_totalZ[i] = 0.0;
		for(j=0;j<T;j++){ //inicializo matriz fuerza local y visitante
				matriz_fuerzaX_l[j*N + i]=0.0;
				matriz_fuerzaZ_l[j*N + i]=0.0;
				matriz_fuerzaY_l[j*N + i]=0.0;
				matriz_fuerzaX_v[j*N + i]=0.0;
				matriz_fuerzaY_v[j*N + i]=0.0;
				matriz_fuerzaZ_v[j*N + i]=0.0;
		}
	}

	pthread_barrier_init(&barrera,NULL,T);

	tIni = dwalltime(); 
	
	if (rank==0){
		procesoA(N,delta_tiempo,pasos,T,MPI_CUERPO,hilo);
	}else if (rank==1){
		procesoB(N,delta_tiempo,pasos,T,MPI_CUERPO,hilo);
	}
   	tFin =	dwalltime();
	printf("tiempo %f \n", (tFin - tIni));
	
	finalizar(MPI_CUERPO);
	return tFin - tIni;
}

void inicializarEstrella(cuerpo_t *cuerpo,int i,double n,float *toroide_alfa,float *toroide_theta,
						float *toroide_incremento,float *toroide_lado,float *toroide_r,float *toroide_R){
    cuerpo->masa = 0.001*8;

        if ((*toroide_alfa + *toroide_incremento) >=2*M_PI){
            *toroide_alfa = 0;
            *toroide_theta += *toroide_incremento;
        }else{
            *toroide_alfa+=*toroide_incremento;
        }

	cuerpo->px = (*toroide_R + *toroide_r*cos(*toroide_alfa))*cos(*toroide_theta);
	cuerpo->py = (*toroide_R + *toroide_r*cos(*toroide_alfa))*sin(*toroide_theta);
	cuerpo->pz = *toroide_r*sin(*toroide_alfa);

    cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarPolvo(cuerpo_t *cuerpo,int i,double n,float *toroide_alfa,float *toroide_theta,
						float *toroide_incremento,float *toroide_lado,float *toroide_r,float *toroide_R){

    cuerpo->masa = 0.001*4;
	
        if ((*toroide_alfa + *toroide_incremento) >=2*M_PI){
            *toroide_alfa = 0;
            *toroide_theta += *toroide_incremento;
        }else{
            *toroide_alfa+=*toroide_incremento;
        }

	cuerpo->px = (*toroide_R + *toroide_r*cos(*toroide_alfa))*cos(*toroide_theta);
	cuerpo->py = (*toroide_R + *toroide_r*cos(*toroide_alfa))*sin(*toroide_theta);
	cuerpo->pz = *toroide_r*sin(*toroide_alfa);
	
	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;
    
	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 0.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 0.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarH2(cuerpo_t *cuerpo,int i,double n,float *toroide_alfa,float *toroide_theta,
					float *toroide_incremento,float *toroide_lado,float *toroide_r,float *toroide_R){

    cuerpo->masa = 0.001;

	if ((*toroide_alfa + *toroide_incremento) >=2*M_PI){
            *toroide_alfa = 0;
            *toroide_theta += *toroide_incremento;
	}else{
            *toroide_alfa+=*toroide_incremento;
	}

	cuerpo->px = (*toroide_R + *toroide_r*cos(*toroide_alfa))*cos(*toroide_theta);
	cuerpo->py = (*toroide_R + *toroide_r*cos(*toroide_alfa))*sin(*toroide_theta);
	cuerpo->pz = *toroide_r*sin(*toroide_alfa);

	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}



int main(int argc, char * argv[]) {
	int N,delta_tiempo,pasos,T;
	float toroide_alfa;
	float toroide_theta;
	float toroide_incremento;
	float toroide_lado;
	float toroide_r;
	float toroide_R;
	cuerpo_t *cuerpos;
	int rank;
	double tTotal;

	if (argc < 4){
		printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos>\n",argv[0]);
		return -1;
	}
        
	
	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	N = atoi(argv[1]);
	delta_tiempo = atof(argv[2]);
	pasos = atoi(argv[3]);
   	T = atoi(argv[4]);
	

    	tTotal = funcion_mpi(rank,N,delta_tiempo,pasos,T,&toroide_alfa,&toroide_theta,
		&toroide_incremento,&toroide_lado,&toroide_r,&toroide_R);

	MPI_Finalize();

	printf("Tiempo en segundos: %f\n",tTotal);

	
    return(0);
}
void inicializarCuerpos(cuerpo_t *cuerpos,int N,float *toroide_alfa,float *toroide_theta,
						float *toroide_incremento,float *toroide_lado,float *toroide_r,float *toroide_R){
 int cuerpo,i;
 double n = N;

	*toroide_alfa = 0.0;
	*toroide_theta = 0.0;
	*toroide_lado = sqrt(N);
	*toroide_incremento = 2*M_PI / *toroide_lado;
	*toroide_r = 1.0;
	*toroide_R = 2*(*toroide_r);
	
	srand(1);

	for(cuerpo = 0; cuerpo < N; cuerpo++){
	
        /*fuerza_totalX[cuerpo] = 0.0; TODO ESTO LO HAGO EN FUNCION_MPI
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

		for(i=0;i<T;i++){ 
			matriz_fuerzaX[i*N + cuerpo]=0;
			matriz_fuerzaZ[i*N + cuerpo]=0;
			matriz_fuerzaY[i*N + cuerpo]=0;
		}*/

		cuerpos[cuerpo].cuerpo = (rand() %3);

		if (cuerpos[cuerpo].cuerpo == ESTRELLA){
			inicializarEstrella(&cuerpos[cuerpo],cuerpo,n,toroide_alfa,toroide_theta,
				toroide_incremento,toroide_lado,toroide_r,toroide_R);
		}else if (cuerpos[cuerpo].cuerpo == POLVO){
			inicializarPolvo(&cuerpos[cuerpo],cuerpo,n,toroide_alfa,toroide_theta,
				toroide_incremento,toroide_lado,toroide_r,toroide_R);
		}else if (cuerpos[cuerpo].cuerpo == H2){
			inicializarH2(&cuerpos[cuerpo],cuerpo,n,toroide_alfa,toroide_theta,
				toroide_incremento,toroide_lado,toroide_r,toroide_R);
		}

	}

		cuerpos[0].masa = 2.0e2;
	    cuerpos[0].px = 0.0;
		cuerpos[0].py = 0.0;
		cuerpos[0].pz = 0.0;
		cuerpos[0].vx = -0.000001;
		cuerpos[0].vy = -0.000001;
		cuerpos[0].vz = 0.0;

		cuerpos[1].masa = 1.0e1;
	    cuerpos[1].px = -1.0;
		cuerpos[1].py = 0.0;
		cuerpos[1].pz = 0.0;
		cuerpos[1].vx = 0.0;
		cuerpos[1].vy = 0.0001;
		cuerpos[1].vz = 0.0;
}

void finalizar(MPI_Datatype MPI_CUERPO){
	printf("cerrando cuerpos\n");
	free(cuerpos);
	printf("cerrando matrices\n");
	free(fuerza_totalX);
	free(fuerza_totalY);
	free(fuerza_totalZ);
	free(matriz_fuerzaX_l);
	free(matriz_fuerzaY_l);
	free(matriz_fuerzaZ_l);
	free(matriz_fuerzaX_v);
	free(matriz_fuerzaY_v);
	free(matriz_fuerzaZ_v);
	printf("cerrando mutex\n");
	pthread_mutex_destroy(&mutex1);
	pthread_mutex_destroy(&mutex2);
	printf("cerrando barrera\n");
	pthread_barrier_destroy(&barrera);
	printf("cerrando mpi_cuerpo\n");
	MPI_Type_free(&MPI_CUERPO);
}
