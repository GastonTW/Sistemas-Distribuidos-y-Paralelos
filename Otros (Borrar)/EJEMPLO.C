#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

typedef double tipo;

double dwalltime() {
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

tipo* mm(tipo* parteA, tipo* B,tipo* parteC,int N,int filas){
	int i,j,k;
	int iN;
	int jN;
	tipo suma;
	for (i = 0; i < filas; i++) {
		iN = i*N;
		for (j = 0; j < N; j++) {
			jN = j*N;
			suma = 0;
			for (k = 0; k < N; k++) {
				suma += parteA[iN + k] * B[k + jN];
			}
			parteC[iN+j] = suma;
		}
	}
	return parteC;
}

void proceso(int ID, int P, int N){
	tipo *A, *B, *C; 
	tipo *parteA, *parteC;
	double timetick;
   	int check = 1;
	int total = N*N, parte = total/P, filas = N/P, root = 0;

	parteA = (tipo *)malloc(sizeof(tipo)*parte);
	parteC = (tipo *)malloc(sizeof(tipo)*parte);
	B = (tipo *)malloc(sizeof(tipo)*total);
	
	if (ID==root){
		
		A = (tipo*)malloc(sizeof(tipo)*total);
		C = (tipo*)malloc(sizeof(tipo)*total);
		for (int i = 0; i < total; i++) {
		    A[i] = 1;
		    B[i] = 1;
		}
		timetick = dwalltime();
	}
	
	MPI_Scatter(A, parte, MPI_DOUBLE, parteA, parte, MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Bcast(B, total, MPI_DOUBLE, root, MPI_COMM_WORLD);
	parteC = mm(parteA,B,parteC,N,filas);
	MPI_Gather(parteC, parte, MPI_DOUBLE, C, parte, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	free(parteA);
	free(parteC);
	if (ID==root){
		free(A);
		printf("Tiempo en segundos: %f\n", dwalltime() - timetick);
		// Verificación
		/*for (int i = 0; i < N; i++) {
		    for (int j = 0; j < N; j++) {
			check = check && (C[i * N + j] == N);
		    }
		}
		if (check) {
		    printf("Multiplicacion de matrices resultado correcto\n");
		}else {
		    printf("Multiplicacion de matrices resultado erroneo\n");}	*/
		free(C);
	}
	free(B);
	
}



int main( int argc, char *argv[]){
	int P;
	int ID;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &ID );
	MPI_Comm_size(MPI_COMM_WORLD,&P);
	int N = atoi(argv[1]);
	
	proceso(ID,P,N);
	
	MPI_Finalize();
	return 0;
}

