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
int N, T, delta_tiempo, pasos, proceso;
pthread_mutex_t mutex1, mutex2;
pthread_barrier_t barrera;
MPI_Datatype MPI_CUERPO;

// Prototipos
double dwalltime();
void crear_tipo_cuerpo(MPI_Datatype *tipo);
void inicializarEstrella(cuerpo_t *cuerpo,int i,double n);
void inicializarPolvo(cuerpo_t *cuerpo,int i,double n);
void inicializarH2(cuerpo_t *cuerpo,int i,double n);
void inicializarCuerpos(cuerpo_t *cuerpos,int N);
void finalizar();
int funcion_mpi(int rank,int N,int delta_tiempo,int pasos,int T);
void procesoA(pthread_t* hilo);
void procesoB(pthread_t* hilo);
void crear_hilos(int proceso_p,pthread_t* hilo_p);
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
	free(cuerpos);
	free(fuerza_totalX);
	free(fuerza_totalY);
	free(fuerza_totalZ);
	free(matriz_fuerzaX_l);
	free(matriz_fuerzaY_l);
	free(matriz_fuerzaZ_l);
	free(matriz_fuerzaX_v);
	free(matriz_fuerzaY_v);
	free(matriz_fuerzaZ_v);
	pthread_mutex_destroy(&mutex1);
	pthread_mutex_destroy(&mutex2);
	pthread_barrier_destroy(&barrera);
	MPI_Type_free(&MPI_CUERPO);
}

int funcion_mpi(int rank,int N_p,int dt_p,int pasos_p,int T_p){
	int i,j;
	double tIni, tFin;
	pthread_t hilo[T];
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
	
	crear_tipo_cuerpo(&MPI_CUERPO);
    
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

	pthread_barrier_init(&barrera,NULL,T);
    // Fin de la inicializacion

	tIni = dwalltime();

	MPI_Bcast(cuerpos, N, MPI_CUERPO, 0, MPI_COMM_WORLD); 

	if (rank==0){
		procesoA(hilo);
	}else if (rank==1){
		procesoB(hilo);
	}

   	tFin =	dwalltime();
	printf("tiempo %f \n", (tFin - tIni));
	
	finalizar();
	return tFin - tIni;
}

void procesoA(pthread_t* hilo){
	int i,j,c,paso;
    int nt = N * T;
	pthread_mutex_lock(&mutex1);
	pthread_mutex_lock(&mutex2);

	crear_hilos(0,hilo);

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
		pthread_mutex_unlock(&mutex2);
	}
	printf("cerrando hilos\n");
	cerrar_hilos();
	printf("hilos cerrados\n");
}

void procesoB(pthread_t* hilo){
	int i,j,c,paso;
    int nt = N * T;
	pthread_mutex_lock(&mutex1);
	pthread_mutex_lock(&mutex2);
	crear_hilos(1,hilo);
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
		MPI_Recv(cuerpos, N, MPI_CUERPO, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // Revisar
		MPI_Ssend(cuerpos, N, MPI_CUERPO, 0, paso, MPI_COMM_WORLD);
		pthread_mutex_unlock(&mutex2);
    }
    cerrar_hilos();
}

void crear_hilos(int proceso_p,pthread_t* hilo_p){
    int i;
    proceso=proceso_p;
    hilo=hilo_p;
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
	    //if(id==0) printf("paso: %d \n",paso);
        calcularFuerzas(id);
        pthread_barrier_wait(&barrera);
        if ((id == 0) || (id == 1)){	
            pthread_mutex_unlock(&mutex1);//semaforo (aviso a mpi que termine
            pthread_mutex_lock(&mutex2);//semaforo (espera a mpi)
        }
	    pthread_barrier_wait(&barrera);
        moverCuerpos(id);
        pthread_barrier_wait(&barrera); //barrera
        if ((id == 0) || (id == 1)) {
            pthread_mutex_unlock(&mutex1);//semaforo (aviso a mpi que termine)
            pthread_mutex_lock(&mutex2);//semaforo (espera a mpi)
        }
	    pthread_barrier_wait(&barrera);
        if (((paso == 0) || (paso == 999)) && id == 0){
			printf("Pos cuerpo 1: X(%f) Y(%f) Z(%f)\n",cuerpos[0].px,cuerpos[0].py,cuerpos[0].pz);
			printf("Pos cuerpo 2: X(%f) Y(%f) Z(%f)\n",cuerpos[1].px,cuerpos[1].py,cuerpos[1].pz);
        }
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
			//printf ("F(%f) = (G(%f)*cuerpos[cuerpo1].masa(%f)*cuerpos[cuerpo2].masa(%f))/(r(%f)*r(%f)) \n",F,G,cuerpos[cuerpo1].masa,cuerpos[cuerpo2].masa,r,r);
			//printf("dif cuerpos l: X(%f) Y(%f) Z(%f)\n",dif_X,dif_Y,dif_Z);
			matriz_fuerzaX_l[id*N+cuerpo1] += dif_X;
	                matriz_fuerzaY_l[id*N+cuerpo1] += dif_Y;
	                matriz_fuerzaZ_l[id*N+cuerpo1] += dif_Z;
			//printf("fuerzas l: X(%f) Y(%f) Z(%f)\n",matriz_fuerzaX_l[id*N+cuerpo1],matriz_fuerzaY_l[id*N+cuerpo1],matriz_fuerzaZ_l[id*N+cuerpo1]);
		

	                matriz_fuerzaX_l[id*N+cuerpo2] -= dif_X;
	                matriz_fuerzaY_l[id*N+cuerpo2] -= dif_Y;
	                matriz_fuerzaZ_l[id*N+cuerpo2] -= dif_Z;
	                
		}
	}
}

void moverCuerpos(int id){
 	int cuerpo,i,j;
	int t2=2*T;
	for(cuerpo = id; cuerpo<N ; cuerpo+=t2){
		for (i=0;i<T;i++){
			//printf("matrizfuerza_local : %f ",matriz_fuerzaX_l[i*N+cuerpo]);
			//printf("matrizfuerza_visitante : %f \n",matriz_fuerzaX_v[i*N+cuerpo]);
			fuerza_totalX[cuerpo] += (matriz_fuerzaX_l[i*N+cuerpo]+matriz_fuerzaX_v[i*N+cuerpo]);
			
			fuerza_totalY[cuerpo] += (matriz_fuerzaY_l[i*N+cuerpo]+matriz_fuerzaY_v[i*N+cuerpo]);
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

        cuerpos[cuerpo].vx += fuerza_totalX[cuerpo]*delta_tiempo;
        cuerpos[cuerpo].vy += fuerza_totalY[cuerpo]*delta_tiempo;
        //cuerpos[cuerpo].vz += fuerza_totalZ[cuerpo]*dt;

        cuerpos[cuerpo].px += cuerpos[cuerpo].vx *delta_tiempo;
        cuerpos[cuerpo].py += cuerpos[cuerpo].vy *delta_tiempo;
        //cuerpos[cuerpo].pz += cuerpos[cuerpo].vz *dt;

		fuerza_totalX[cuerpo] = 0.0;
		fuerza_totalY[cuerpo] = 0.0;
		fuerza_totalZ[cuerpo] = 0.0;

	}
}