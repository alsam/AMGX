% `AmgX` working notes on the **A**lgebraic **M**ulti**G**rid Library.
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

**Running single GPU example from the build directory:**
========================================================

```sh
cd build_debug
examples/amgx_capi -m ../examples/matrix.mtx -c ../core/configs/FGMRES_AGGREGATION.json
AMGX version 2.1.0.131-opensource
Built on Feb 22 2021, 12:11:39
Compiled with CUDA Runtime 11.2, using CUDA driver 11.2
Warning: No mode specified, using dDDI by default.
Reading data...
RHS vector was not found. Using RHS b=[1,…,1]^T
Solution vector was not found. Setting initial solution to x=[0,…,0]^T
Finished reading
```

**Running single GPU example:** contd.
======================================

```sh
AMG Grid:
 Number of Levels: 1
    LVL         ROWS               NNZ    SPRSTY       Mem (GB)
 --------------------------------------------------------------
   0(D)           12                61     0.424       8.75e-07
 --------------------------------------------------------------
 Grid Complexity: 1
 Operator Complexity: 1
 Total Memory Usage: 8.75443e-07 GB
```

**Running single GPU example:** end
===================================
```sh
 --------------------------------------------------------------
   iter      Mem Usage (GB)       residual           rate
 --------------------------------------------------------------
    Ini             1.73499   3.464102e+00
      0             1.73499   2.183701e-14         0.0000
 --------------------------------------------------------------
 Total Iterations: 1
 Avg Convergence Rate: 		         0.0000
 Final Residual: 		   2.183701e-14
 Total Reduction in Residual: 	   6.303801e-15
 Maximum Memory Usage: 		          1.735 GB
 --------------------------------------------------------------
Total Time: 0.283992
    setup: 0.10545 s
    solve: 0.178542 s
    solve(per iteration): 0.178542 s
```

Debug session
=============

```sh
cgdb examples/amgx_capi
[?2004hGNU gdb (GDB) 10.1
...
Reading symbols from examples/amgx_capi...
[?2004h(gdb) b main
line 182.[?2004h(gdb) r -m ../examples/matrix.mtx -c ../core/configs/FGMRES_AGGREGATION.json
...
434│     /* solver setup */
435│     AMGX_solver_setup(solver, A);
436│     /* solver solve */
437├───> AMGX_solver_solve(solver, b, x);
```

Inside `AMGX_solver_solve`
==========================

```cpp
2780│     AMGX_RC AMGX_API AMGX_solver_solve(AMGX_solver_handle slv, AMGX_vector_handle rhs, AMGX_vector_handle sol)
2781│     {
2782│         nvtxRange nvrf(__func__);
2783│
2784│         AMGX_CPU_PROFILER( "AMGX_solver_solve " );
2785│         Resources *resources;
2786│         AMGX_CHECK_API_ERROR(getAMGXerror(getResourcesFromSolverHandle(slv, &resources)), NULL);
2787│         AMGX_ERROR rc = AMGX_OK;
2788│
2789│         try
2790│         {
2791├───────────> AMGX_Mode mode = get_mode_from<AMGX_solver_handle>(slv);
...
```

Inside `AMGX_solver_solve` contd.
=================================

```cpp
2795│ #define AMGX_CASE_LINE(CASE) case CASE: { \
2796│       AMGX_ERROR rcs = solve_with<CASE,AMG_Solver,Vector>(slv, rhs, sol, resources, false); \
2797│       AMGX_CHECK_API_ERROR(rcs, resources); break;\
2798│           }
2799├───────────────────> AMGX_FORALL_BUILDS(AMGX_CASE_LINE)
2800│                     AMGX_FORCOMPLEX_BUILDS(AMGX_CASE_LINE)
2801│ #undef AMGX_CASE_LINE
```
```c++
 709│ template<AMGX_Mode CASE,
 710│          template<typename> class SolverType,
 711│          template<typename> class VectorType>
 712├>inline AMGX_ERROR solve_with(AMGX_solver_handle slv,
 713│                              AMGX_vector_handle rhs,
 714│                              AMGX_vector_handle sol,
 715│                              Resources *resources,
 716│                              bool xIsZero = false)
 717│ {
 718│     typedef SolverType<typename TemplateMode<CASE>::Type> SolverLetterT;
 ...
 746│     cudaSetDevice(solver.getResources()->getDevice(0));
 747├───> AMGX_ERROR ret = solver.solve(b, x, wrapSolver.last_solve_status(), xIsZero);
 748│     return ret;
 749│ }
.../AMGX/base/src/amgx_c.cu
```

Inside `AMGX_solver_solve` contd.
=================================

```cpp
282│ /****************************************************
283│ * Solves the AMG system Ax=b
284│ ***************************************************/
285│ template<class T_Config>
286│ AMGX_ERROR AMG_Solver<T_Config>::solve( Vector<T_Config> &b, Vector<T_Config> &x, AMGX_STATUS &status, bool xIsZero )
287│ {
288├───> profilePhaseSolve();
...
296├───> AMGX_ERROR e = solver->solve_no_throw( b, x, status, xIsZero );
297│     thrust::global_thread_handle::cudaFreeWait();
...
.../AMGX/base/src/amg_solver.cu
```


Inside `AMGX_solver_solve` contd.
=================================
```cpp
 974│ template<class TConfig>
 975├>AMGX_ERROR Solver<TConfig>::solve_no_throw(VVector &b, VVector &x,
 976│         AMGX_STATUS &status, bool xIsZero)
 977│ {
...
 991│
 992├───────────> status = this->solve(b, x, xIsZero);
.../AMGX/base/src/solvers/solver.cu
```

Inside `AMGX_solver_solve` contd.
=================================
```cpp
 588│ template<class TConfig>
 589│ AMGX_STATUS Solver<TConfig>::solve(Vector<TConfig> &b, Vector<TConfig> &x,
 590│                                    bool xIsZero)
 591│ {
 592│     PODValueB eps = (sizeof(PODValueB) == 4) ? AMGX_NUMERICAL_SZERO : AMGX_NUMERICAL_DZERO;
 593│     AMGX_CPU_PROFILER("Solver::solve ");
 594│
 595│     if (!m_is_solver_setup)
 596│     {
 597│         FatalError("Error, setup must be called before calling solve",
 598│                    AMGX_ERR_CONFIGURATION);
 599│     }
 600│
 601├───> if (b.get_block_size() != m_A->get_block_dimy())
 602│     {
 603│         FatalError("Block sizes do not match", AMGX_ERR_BAD_PARAMETERS);
 604│     }
```
