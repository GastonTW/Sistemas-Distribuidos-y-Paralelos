# El problema de los N-Cuerpos

Este proyecto implementa una simulación del problema de los N-cuerpos gravitacionales utilizando programación paralela. Se desarrollaron tres versiones del algoritmo:

- 🧮 **Versión secuencial**
- 🧵 **Versión paralela con Pthreads**
- 🌐 **Versión híbrida con MPI + Pthreads**

## 📁 Estructura del repositorio

```
├── Secuencial # Implementación base secuencial
|   ├─ secuencial.c
|   ├─ sec.sh   
├── Pthreads # Implementación con hilos (Pthreads)
|   ├─ pthreads.c
|   ├─ pth.sh           
├── MPI+Pthread # Implementación híbrida (MPI + Pthreads)
|   ├─ MPI+Pthread_source.c
|   ├─ mpi.sh  
├── Informe-Grupo-10.pdf # Informe final del proyecto
```

## 🛠️ Compilación

### Secuencial
```bash
gcc -o secuencial_g10 secuencial.c -lm
```

### Pthreads
```bash
gcc -o pthread_g10 pthreads.c -lpthread -lm
```

### MPI + Pthreads
```bash
mpicc -o mpi_g10 MPI+Pthread_source.c -lm -lpthread
```

## 🚀 Ejecución

### Secuencial
```bash
./secuencial_g10 <N> <DT> <Pasos>
```

### Pthreads
```bash
./pthread_g10 <N> <DT> <Pasos> <HILOS>
```

### MPI + Pthreads
```bash
mpirun -np 2 ./mpi_g10 <N> <DT> <Pasos> <HILOS_POR_PROCESO>
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