# N-Body Parallel Simulation

Este proyecto implementa una simulación del problema de los N-cuerpos gravitacionales utilizando programación paralela. Se desarrollaron tres versiones del algoritmo:

- 🧮 **Versión secuencial**
- 🧵 **Versión paralela con Pthreads**
- 🌐 **Versión híbrida con MPI + Pthreads**

## 📁 Estructura del repositorio

```
├── secuencial.c            # Implementación base secuencial
├── pthreads.c              # Implementación con hilos (Pthreads)
├── MPI+Pthread_source.c    # Implementación híbrida (MPI + Pthreads)
├── Informe-Grupo-10.pdf    # Informe final del proyecto (análisis y resultados)
├── README.md               # Este archivo
```

## 🛠️ Compilación

### Secuencial
```bash
gcc -o nbody_secuencial secuencial.c -lm
```

### Pthreads
```bash
gcc -o nbody_pthreads pthreads.c -lpthread -lm
```

### MPI + Pthreads
```bash
mpicc -o nbody_mpi MPI+Pthread_source.c -lpthread -lm
```

## 🚀 Ejecución

### Secuencial
```bash
./nbody_secuencial <N> <DT> <Pasos>
```

### Pthreads
```bash
./nbody_pthreads <N> <DT> <Pasos> <HILOS>
```

### MPI + Pthreads
```bash
mpirun -np <PROCESOS> ./nbody_mpi <N> <DT> <Pasos> <HILOS_POR_PROCESO>
```

> 📌 Parámetros:
- `N`: cantidad de cuerpos
- `DT`: delta de tiempo
- `Pasos`: cantidad de pasos de simulación
- `HILOS`: cantidad de hilos por proceso

## 📊 Informe

El archivo `Informe-Grupo-10.pdf` documenta:
- Estrategia de paralelización
- Análisis de escalabilidad
- Cálculo de speedup, eficiencia y overhead
- Conclusiones y comparativa de soluciones

## 👨‍💻 Autores

- Tobias Garcia Iacovelli
- Gastón Triviño  
UNLP - Facultad de Informática - Sistemas Distribuidos y Paralelos 2025