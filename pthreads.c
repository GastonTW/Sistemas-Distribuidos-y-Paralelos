// Compilar:
//		gcc -o n_body_simple_NOGL n_body_simple.c -lm
// Ejecutar:
//		./n_body_simple_NOGL <nro de cuerpos> <DT> <Pasos>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

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
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 //Hidrogeno molecular

// ===============
// ===== CPU =====
// ===============

//
// Estructuras y variables para Algoritmo de gravitacion
//
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

float *fuerza_totalX,*fuerza_totalY, *fuerza_totalZ;
float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

cuerpo_t *cuerpos;
int delta_tiempo = 1.0f; //Intervalo de tiempo, longitud de un paso
int pasos;
int N;
//nnuevo
int T;
pthread_mutex_t mutex;
pthread_barrier_t barrera;
float *matriz_fuerzaX,*matriz_fuerzaY,*matriz_fuerzaZ;

//
// Funciones para Algoritmo de gravitacion
//

void calcularFuerzas(cuerpo_t *cuerpos, int N, int dt,int id){
int cuerpo1, cuerpo2;
float dif_X, dif_Y, dif_Z;
float r;
float F;

	for(cuerpo1 = id; cuerpo1<N-1 ; cuerpo1+=T){
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
					/*fuerza_totalX[cuerpo1] += dif_X;
	                fuerza_totalY[cuerpo1] += dif_Y;
	                fuerza_totalZ[cuerpo1] += dif_Z;

	                fuerza_totalX[cuerpo2] -= dif_X;
	                fuerza_totalY[cuerpo2] -= dif_Y;
	                fuerza_totalZ[cuerpo2] -= dif_Z;*/
					matriz_fuerzaX[id*N+cuerpo1] += dif_X;
	                matriz_fuerzaY[id*N+cuerpo1] += dif_Y;
	                matriz_fuerzaZ[id*N+cuerpo1] += dif_Z;

	                matriz_fuerzaX[id*N+cuerpo2] -= dif_X;
	                matriz_fuerzaY[id*N+cuerpo2] -= dif_Y;
	                matriz_fuerzaZ[id*N+cuerpo2] -= dif_Z;
	                
		}
	}
}

void moverCuerpos(cuerpo_t *cuerpos, int N, int dt,int id){
 	int cuerpo,i;
	
	for(cuerpo = id; cuerpo<N ; cuerpo+=T){

		for (i=0;i<T;i++){
			fuerza_totalX[i] += matriz_fuerzaX[i*N+cuerpo];
			fuerza_totalY[i] += matriz_fuerzaY[i*N+cuerpo];
			fuerza_totalZ[i] += matriz_fuerzaZ[i*N+cuerpo];
			//PONER MATRIZ FUERZA EN 0

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

void gravitacionCPU(cuerpo_t *cuerpos, int N, int dt){
	calcularFuerzas(cuerpos,N,dt);
	moverCuerpos(cuerpos,N,dt);
}
void gravitacion(void *arg){
	int id=(int*)arg;
	int paso;
	for(paso=0; paso<pasos;(paso+delta_tiempo)){
		calcularFuerzas(cuerpos,N,delta_tiempo,id);
		pthread_barrier_wait(&barrera);
		moverCuerpos(cuerpos,N,delta_tiempo,id);	
		pthread_barrier_wait(&barrera);
	}
	
}

void inicializarEstrella(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*8;

        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

    	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

		cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
		cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
		cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarPolvo(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001*4;
	
        if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
        }else{
            toroide_alfa+=toroide_incremento;
        }

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);
	
	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;
    
	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 0.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 0.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarH2(cuerpo_t *cuerpo,int i,double n){

    cuerpo->masa = 0.001;

	if ((toroide_alfa + toroide_incremento) >=2*M_PI){
            toroide_alfa = 0;
            toroide_theta += toroide_incremento;
	}else{
            toroide_alfa+=toroide_incremento;
	}

	cuerpo->px = (toroide_R + toroide_r*cos(toroide_alfa))*cos(toroide_theta);
	cuerpo->py = (toroide_R + toroide_r*cos(toroide_alfa))*sin(toroide_theta);
	cuerpo->pz = toroide_r*sin(toroide_alfa);

	cuerpo->vx = 0.0;
	cuerpo->vy = 0.0;
	cuerpo->vz = 0.0;

	cuerpo->r = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->g = 1.0; //(double )rand()/(RAND_MAX+1.0);
	cuerpo->b = 1.0; //(double )rand()/(RAND_MAX+1.0);
}

void inicializarCuerpos(cuerpo_t *cuerpos,int N){
 int cuerpo,i;
 double n = N;

	toroide_alfa = 0.0;
	toroide_theta = 0.0;
	toroide_lado = sqrt(N);
	toroide_incremento = 2*M_PI / toroide_lado;
	toroide_r = 1.0;
	toroide_R = 2*toroide_r;
	
	srand(time(NULL));

	for(cuerpo = 0; cuerpo < N; cuerpo++){
	
        fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

		for(i;i<T;i++){ //inicializo matriz fuerza total
			matriz_fuerzaX[i*N + cuerpo]=0;
			matriz_fuerzaZ[i*N + cuerpo]=0;
			matriz_fuerzaY[i*N + cuerpo]=0;
		}

		cuerpos[cuerpo].cuerpo = (rand() %3);

		if (cuerpos[cuerpo].cuerpo == ESTRELLA){
			inicializarEstrella(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == POLVO){
			inicializarPolvo(&cuerpos[cuerpo],cuerpo,n);
		}else if (cuerpos[cuerpo].cuerpo == H2){
			inicializarH2(&cuerpos[cuerpo],cuerpo,n);
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

void finalizar(void){
	free(cuerpos);
	free(fuerza_totalX);
	free(fuerza_totalY);
	free(fuerza_totalZ);
	free(matriz_fuerzaX);
	free(matriz_fuerzaY);
	free(matriz_fuerzaZ);
	pthread_barrier_destroy(&barrera);
	pthread_mutex_destroy(&mutex);
}

int main(int argc, char * argv[]) {

	if (argc < 4){
		printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos>\n",argv[0]);
		return -1;
	}
	
	N = atoi(argv[1]);
	delta_tiempo = atof(argv[2]);
	pasos = atoi(argv[3]);
    T = atoi(argv[4]);
	pthread_t hilo[T];
	int hilo_id[T];

	

	pthread_mutex_init(&mutex,NULL);
    pthread_barrier_init(&barrera,NULL,T);

	cuerpos = (cuerpo_t*)malloc(sizeof(cuerpo_t)*N);
	fuerza_totalX = (float*)malloc(sizeof(float)*N);
	fuerza_totalY = (float*)malloc(sizeof(float)*N);
	fuerza_totalZ = (float*)malloc(sizeof(float)*N);
	matriz_fuerzaX = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaY = (float*)malloc(sizeof(float)*N*T);
	matriz_fuerzaZ = (float*)malloc(sizeof(float)*N*T);

	inicializarCuerpos(cuerpos,N);
	
	tIni = dwalltime(); 
	
	for(int id = 0; id<T; id++){
        pthread_create(&hilo[id], NULL, &gravitacion, (void*)id);
    }

	//
    for(int id = 0; id<T; id++){
        pthread_join(hilo[id], NULL);
    }

	tFin =	dwalltime();
	tTotal = tFin - tIni;
	
	printf("Tiempo en segundos: %f\n",tTotal);

	finalizar();
    return(0);

}