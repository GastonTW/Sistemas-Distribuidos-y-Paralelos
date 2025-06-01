float *fuerza_totalX,*fuerza_totalY,*fuerza_totalZ;
float *matriz_fuerzaX_l,*matriz_fuerzaY_l,*matriz_fuerzaZ_l;
float *matriz_fuerzaX_v,*matriz_fuerzaY_v,*matriz_fuerzaZ_v;
cuerpo_t *cuerpos;
int N,T,delta_tiempo,pasos,proceso;

void gravitacion(int id){
    for(i=0;i<){ //CALCULAR LAS FUERZAS Q LE TOCARON A LOS HILOS
        calcularFuerzas();
    }
    //barrera
    //semaforo (aviso a mpi que termine)
    //semaforo (espera a mpi)
     for(i=0;i<){ //MOVER LOS CUERPOS Q LE TOCARON A LOS HILOS
        moverCuerpos();
    }
    //barrera
    //semaforo (aviso a mpi que termine)
    //semaforo (espera a mpi)
}

void crear_hilos(N,T,delta,pasos,fuerzasTODAS,proceso,){
    //INICIALIZAR VARIABLES RECIBIDAS POR PARAMETRO A GLOBALES
    int i;
    for(int i = 0; id<T; i++){
        pthread_create(&hilo[i], NULL, gravitacion, ((void*)i*2)+proceso);
    }

	//
    for(int i = 0; i<T; i++){
        pthread_join(hilo[i], NULL);
    }

}