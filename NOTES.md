% `AmgX` working notes on the **A**lgebraic **Mu**ulti**G**rid Library.
% Alexander Samoilov
% March 2021

build **AMGX**
==============

+ there is an incompability between `CUDA` version 11 and `main` branch of the librarry, so switch to the branch `v2.1.x`

```sh
    git checkout v2.1.x
    mkdir -p build
    cd build
    cmake . -DCUDA_ARCH="35 52 60"
    make -j8
```

Getting started
===============

run examples:

+ **Running single GPU example from the build directory:**

```c++
$ examples/amgx_capi -m ../examples/matrix.mtx -c ../core/configs/FGMRES_AGGREGATION.json
```

+ **Running multi GPU example from the build directory:**

```c++
$ mpirun -n 2 examples/amgx_mpi_capi.exe -m ../examples/matrix.mtx -c ../core/configs/FGMRES_AGGREGATION.json
```

**`AMGX_Mode`**
===============

+ **`AMGX_Mode`** is a solver parameter, it has a four letter coding scheme, the sample is below:

```c++
typedef enum
{
...
    AMGX_mode_hDDI = AMGX_ASSEMBLE_MODE(
                            AMGX_host,
                            AMGX_vecDouble,
                            AMGX_matDouble,
                            AMGX_indInt)
...
```

+ a $1^{st}$ letter either **h**ost or **d**evice.
+ the latter stand for precision either for **D**ouble/**F**loat/**I**nteger/**C**omplex/**Z**DoubleComplex
+ the last 3 letters denote precision for:
    + vector
    + matrix
    + integer index


**AMGX** solver initialization
==============================

+ the solver initialization consists of the following API calls:

```c++
    AMGX_resources_create_simple(&rsrc, cfg);
    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_vector_create(&x, rsrc, mode);
    AMGX_vector_create(&b, rsrc, mode);
    AMGX_solver_create(&solver, rsrc, mode, cfg);
```

**AMGX** read input data
========================

+ the following API calls are used for reading input data:

```c++
    AMGX_read_system(A, b, x, argv[pidx + 1]);
    AMGX_matrix_get_size(A, &n, &bsize_x, &bsize_y);
    AMGX_vector_get_size(x, &sol_size, &sol_bsize);

    if (sol_size == 0 || sol_bsize == 0)
    {
        AMGX_vector_set_zero(x, n, bsize_x);
    }
```

**AMGX** solver invocation
==========================

```c++
    AMGX_solver_setup(solver, A);
    AMGX_config_add_parameters(&cfg,
          "config_version=2, default:tolerance=1e-12");
    AMGX_solver_solve(solver, b, x);
    AMGX_solver_get_status(solver, &status);
```

**howto** print residual history
================================

```c++
    /* example of how to print the residual history */
    int nit;
    double res;
    AMGX_solver_get_iterations_number(solver, &nit);
    for (int i=0; i<nit; i++) {
      printf("residual from iteration %d=", i);
      for (int j=0; j<bsize_y; j++) {
        AMGX_solver_get_iteration_residual(solver, i, j, &res);
        printf("%f ", (float)(res));
      }
      printf("\n");
    }
```

write linear system and destroy
===============================

```c++
    AMGX_write_system(A, b, x, "output.system.mtx");
    /* destroy resources, matrix, vector and solver */
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(x);
    AMGX_vector_destroy(b);
    AMGX_matrix_destroy(A);
    AMGX_resources_destroy(rsrc);
```
