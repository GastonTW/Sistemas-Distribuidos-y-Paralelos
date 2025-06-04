#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

// Define
#define PI (3.141592653589793)
#define G 6.673e-11
#define ESTRELLA 0
#define POLVO 1
#define H2 2 //Hidrogeno molecular

// Definicion de tipos
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

// Variables
double tIni, tFin, tTotal; // Para mostrar tiempo de ejecucion

float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

float *fuerza_totalX, *fuerza_totalY, *fuerza_totalZ;
float *matriz_fuerzaX_l, *matriz_fuerzaY_l, *matriz_fuerzaZ_l;
float *matriz_fuerzaX_v, *matriz_fuerzaY_v, *matriz_fuerzaZ_v;
cuerpo_t *cuerpos;
int N, T, pasos, proceso;
int delta_tiempo = 1.0f;
pthread_mutex_t mutex1= PTHREAD_MUTEX_INITIALIZER, mutex2= PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t barrera;
pthread_t *hilo;

// Prototipos
double dwalltime();
void crear_tipo_cuerpo(MPI_Datatype *tipo);
void inicializarEstrella(cuerpo_t *cuerpo,int i,double n);
void inicializarPolvo(cuerpo_t *cuerpo,int i,double n);
void inicializarH2(cuerpo_t *cuerpo,int i,double n);
void inicializarCuerpos(cuerpo_t *cuerpos,int N);
void finalizar();
double funcion_mpi(int rank,int N,int delta_tiempo,int pasos,int T);
void procesoA();
void procesoB();
void crear_hilos(int proceso_p);
void cerrar_hilos();
void *gravitacion(void *arg);
void calcularFuerzas(int id);
void moverCuerpos(int id);


//Main
int main(int argc, char * argv[]) {
	int rank;
	double tTotal;

	if (argc < 4){
		printf("Ejecutar: %s <nro. de cuerpos> <DT> <pasos>\n",argv[0]);
		return -1;
	}
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
    tTotal = funcion_mpi(rank,atoi(argv[1]),atof(argv[2]),atoi(argv[3]),atoi(argv[4]));

	MPI_Finalize();
	

	printf("Tiempo en segundos: %f\n",tTotal);

	
    return(0);
}
//Funciones
double dwalltime() {
	double sec;
	struct timeval tv;

	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
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
	
	srand(1);

	for(cuerpo = 0; cuerpo < N; cuerpo++){
	
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

void finalizar(){
	
	free(fuerza_totalX);
	free(fuerza_totalY);
	free(fuerza_totalZ);
	free(matriz_fuerzaX_l);
	free(matriz_fuerzaY_l);
	free(matriz_fuerzaZ_l);
	free(matriz_fuerzaX_v);
	free(matriz_fuerzaY_v);
	free(matriz_fuerzaZ_v);
	free(hilo);
	pthread_mutex_destroy(&mutex1);
	pthread_mutex_destroy(&mutex2);
	pthread_barrier_destroy(&barrera);
	free(cuerpos);
}

double funcion_mpi(int rank,int N_p,int dt_p,int pasos_p,int T_p){
	int i,j;
	double tIni, tFin;
	
    N=N_p;
    delta_tiempo=dt_p;
    pasos=pasos_p;
    T=T_p;

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
	hilo= (pthread_t*)malloc(sizeof(pthread_t)*T);
	
    
    if (rank == 0) inicializarCuerpos(cuerpos,N); 
	
	for(i=0;i<N;i++){
		fuerza_totalX[i] = 0.0;
		fuerza_totalY[i] = 0.0;
		fuerza_totalZ[i] = 0.0;
		for(j=0;j<T;j++){ 
				matriz_fuerzaX_l[j*N + i]=0.0;
				matriz_fuerzaZ_l[j*N + i]=0.0;
				matriz_fuerzaY_l[j*N + i]=0.0;
				matriz_fuerzaX_v[j*N + i]=0.0;
				matriz_fuerzaY_v[j*N + i]=0.0;
				matriz_fuerzaZ_v[j*N + i]=0.0;
		}
	}

	pthread_barrier_init(&barrera,NULL,T+1);
    // Fin de la inicializacion

	tIni = dwalltime();

	MPI_Bcast(cuerpos, N*(sizeof(cuerpo_t)), MPI_BYTE, 0, MPI_COMM_WORLD); 

	if (rank==0){
		procesoA();
	}else if (rank==1){
		procesoB();
	}

   	tFin =	dwalltime();
	if (rank==0) printf("- Tiempo en segundos: %f\n",tFin-tIni);
	finalizar();
	return tFin - tIni;
}

void procesoA(){
	int i,j,c,paso;
    	int nt = N * T;
	int n2=N/2;

	crear_hilos(0);

	for(paso=0; paso<pasos;paso++){
		pthread_barrier_wait(&barrera);
		MPI_Send(matriz_fuerzaX_l, nt, MPI_FLOAT, 1, paso, MPI_COMM_WORLD);
		MPI_Send(matriz_fuerzaY_l, nt, MPI_FLOAT, 1, paso, MPI_COMM_WORLD);
		MPI_Send(matriz_fuerzaZ_l, nt, MPI_FLOAT, 1, paso, MPI_COMM_WORLD);
        	MPI_Recv(matriz_fuerzaX_v, nt, MPI_FLOAT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaY_v, nt, MPI_FLOAT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaZ_v, nt, MPI_FLOAT, 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		pthread_barrier_wait(&barrera);
		pthread_barrier_wait(&barrera);
		MPI_Allgather(cuerpos, n2 * sizeof(cuerpo_t), MPI_BYTE,cuerpos, n2 * sizeof(cuerpo_t), MPI_BYTE, MPI_COMM_WORLD);

		
		pthread_barrier_wait(&barrera);
		
	}

	cerrar_hilos();

}

void procesoB(){
	int i,j,c,paso;
   	int nt = N * T;
	unsigned n2=N/2;
	pthread_mutex_lock(&mutex1);
	pthread_mutex_lock(&mutex2);
	crear_hilos(1);
	for(paso=0; paso<pasos;paso++){
		
		pthread_barrier_wait(&barrera);
		MPI_Recv(matriz_fuerzaX_v, nt, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaY_v, nt, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(matriz_fuerzaZ_v, nt, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(matriz_fuerzaX_l, nt, MPI_FLOAT, 0, paso, MPI_COMM_WORLD);
		MPI_Send(matriz_fuerzaY_l, nt, MPI_FLOAT, 0, paso, MPI_COMM_WORLD);
		MPI_Send(matriz_fuerzaZ_l, nt, MPI_FLOAT, 0, paso, MPI_COMM_WORLD);
		pthread_barrier_wait(&barrera);
		pthread_barrier_wait(&barrera);
		MPI_Allgather((cuerpos+n2), n2 * sizeof(cuerpo_t), MPI_BYTE,cuerpos, n2 * sizeof(cuerpo_t), MPI_BYTE, MPI_COMM_WORLD);
		pthread_barrier_wait(&barrera);
    }
    cerrar_hilos();

}

void crear_hilos(int proceso_p){
    int i;
    proceso=proceso_p;
    for(int i = 0; i<T; i++){
        pthread_create(&hilo[i], NULL, gravitacion, (void*)(i*2 + proceso));
    }
}

void cerrar_hilos(){
    int i;	
    for(int i = 0; i<T; i++){
        pthread_join(hilo[i], NULL);
    }
}

void *gravitacion(void *arg){
    //int id=(int*)arg; //cambie esto
    int id = (int*)arg;
    int paso;
    for(paso=0; paso<pasos;paso++){
        //CALCULAR LAS FUERZAS Q LE TOCARON A LOS HILOS
        calcularFuerzas(id);
        pthread_barrier_wait(&barrera);
        pthread_barrier_wait(&barrera);
        moverCuerpos(id);
        pthread_barrier_wait(&barrera); //barrera
        pthread_barrier_wait(&barrera);
    }
}

void calcularFuerzas(int id){
    int cuerpo1, cuerpo2;
    float dif_X, dif_Y, dif_Z;
    float r;
    float F;
    int t2=T*2;
    int idN=id*N;
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
		
			matriz_fuerzaX_l[idN+cuerpo1] += dif_X;
	                matriz_fuerzaY_l[idN+cuerpo1] += dif_Y;
	                matriz_fuerzaZ_l[idN+cuerpo1] += dif_Z;

	                matriz_fuerzaX_l[idN+cuerpo2] -= dif_X;
	                matriz_fuerzaY_l[idN+cuerpo2] -= dif_Y;
	                matriz_fuerzaZ_l[idN+cuerpo2] -= dif_Z;
	                
		}
	}
}

void moverCuerpos(int id){
 	int cuerpo,i,j;
	int n2t=N/(2*T);
	int inicio=((id%2)*(N/2))+(((id%(2*T))-(id%2))*(N/(4*T)));
	int fin=inicio+(N/(4*T));
	for(cuerpo = inicio; cuerpo<fin ; cuerpo++){
		for (i=0;i<T;i++){
			fuerza_totalX[cuerpo] += (matriz_fuerzaX_l[i*N+cuerpo]+matriz_fuerzaX_v[i*N+cuerpo]);
			fuerza_totalY[cuerpo] += (matriz_fuerzaY_l[i*N+cuerpo]+matriz_fuerzaY_v[i*N+cuerpo]);
			//fuerza_totalZ[i] += (matriz_fuerzaZ_l[i*N+cuerpo]+matriz_fuerzaZ_v[i*N+cuerpo]);
			
			matriz_fuerzaX_l[i*N+cuerpo] = 0.0f;
			matriz_fuerzaY_l[i*N+cuerpo] = 0.0f;
			matriz_fuerzaZ_l[i*N+cuerpo] = 0.0f;
            		matriz_fuerzaX_v[i*N+cuerpo] = 0.0f;
			matriz_fuerzaY_v[i*N+cuerpo] = 0.0f;
			matriz_fuerzaZ_v[i*N+cuerpo] = 0.0f;
			//matriz_fuerzaZ[i*N+cuerpo] = 0.0f;
		}
		
		
                fuerza_totalX[cuerpo] *= 1/cuerpos[cuerpo].masa;
                fuerza_totalY[cuerpo] *= 1/cuerpos[cuerpo].masa;

                cuerpos[cuerpo].vx += fuerza_totalX[cuerpo]*delta_tiempo;
                cuerpos[cuerpo].vy += fuerza_totalY[cuerpo]*delta_tiempo;
                //cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*delta_tiempo;

                cuerpos[cuerpo].px += cuerpos[cuerpo].vx *delta_tiempo;
                cuerpos[cuerpo].py += cuerpos[cuerpo].vy *delta_tiempo;
                //cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *delta_tiempo;

	        fuerza_totalX[cuerpo] = 0.0f;
	        fuerza_totalY[cuerpo] = 0.0f;
	        fuerza_totalZ[cuerpo] = 0.0f;

	}
}
