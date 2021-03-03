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

Method to compute residual
==========================
```c++
 192│ // Method to compute residual
 193│ template<class TConfig>
 194├>void Solver<TConfig>::compute_residual(const VVector &b, VVector &x)
 195│ {
 196│     AMGX_CPU_PROFILER( "Solver::compute_residual_bx " );
 197│     assert(m_A);
 198│     assert(m_r); //r and b/x are not the same size
 199│     int size, offset;
 200│     m_A->getOffsetAndSizeForView(OWNED, &offset, &size);
 201│     // r = b - Ax.
 202│     m_A->apply(x, *m_r);
 203│     axpby(b, *m_r, *m_r, types::util<ValueTypeB>::get_one(), types::util<ValueTypeB>::get_minus_one(), offset, size);
 204│ }
```

**`thrust`** is used for `axpby`
================================
```c++
 321│ //out=a*x+b*y
 322│ template<class Vector, class Scalar>
 323├>void axpby(const Vector &x, const Vector &y, Vector &out, Scalar a, Scalar b, int offset, int size)
 324│ {
 325│     if (size == -1) { size = x.size() / x.get_block_size(); }
 326│
 327│ #ifndef NDEBUG
 328│
 329│     if (x.get_block_dimx() == -1) { FatalError("x block dims not set", AMGX_ERR_NOT_IMPLEMENTED); }
 330│
 331│     if (y.get_block_dimx() == -1) { FatalError("y block dims not set", AMGX_ERR_NOT_IMPLEMENTED); }
 332│
 333│     if (out.get_block_dimx() == -1) { FatalError("out block dims not set", AMGX_ERR_NOT_IMPLEMENTED); }
 334│
 335│ #endif
 336│     thrust_axpby(x.begin() + offset * x.get_block_size(),
 337│                  x.begin() + (offset + size) * x.get_block_size(),
 338│                  y.begin() + offset * x.get_block_size(),
 339│                  out.begin() + offset * x.get_block_size(),
 340│                  a, b);
 341│     out.dirtybit = 1;
 342│     cudaCheckError();
 343│ }
.../AMGX/base/src/blas.cu
```

callstack snapshot
==================
```c++
[?2004h[?2004l[?2004h(gdb) where
[?2004l#0  amgx::axpby<amgx::Vector<amgx::TemplateConfig<(AMGX_MemorySpace)1, (AMGX_VecPrecision)0, (AMGX_MatPrecision)0, (AMGX_IndPrecision)2> >, double> (x=..., y=..., out=..., a=1, b=-1, offset=0, size=12) at .../AMGX/base/src/blas.cu:323
#1  0x00007fff8430ed15 in amgx::Solver<amgx::TemplateConfig<(AMGX_MemorySpace)1, (AMGX_VecPrecision)0, (AMGX_MatPrecision)0, (AMGX_IndPrecision)2> >::compute_residual (this=0x5555f2bb3ba0, b=..., x=...) at .../AMGX/base/src/solvers/solver.cu:203
#2  0x00007fff843120a6 in amgx::Solver<amgx::TemplateConfig<(AMGX_MemorySpace)1, (AMGX_VecPrecision)0, (AMGX_MatPrecision)0, (AMGX_IndPrecision)2> >::solve (this=0x5555f2bb3ba0, b=..., x=..., xIsZero=false) at .../AMGX/base/src/solvers/solver.cu:691
#3  0x00007fff84314523 in amgx::Solver<amgx::TemplateConfig<(AMGX_MemorySpace)1, (AMGX_VecPrecision)0, (AMGX_MatPrecision)0, (AMGX_IndPrecision)2> >::solve_no_throw (this=0x5555f2bb3ba0, b=..., x=..., status=@0x5555fc8934c0: amgx::AMGX_ST_ERROR, xIsZero=false) at .../AMGX/base/src/solvers/solver.cu:992
#4  0x00007fff82e1d798 in amgx::AMG_Solver<amgx::TemplateConfig<(AMGX_MemorySpace)1, (AMGX_VecPrecision)0, (AMGX_MatPrecision)0, (AMGX_IndPrecision)2> >::solve (this=0x5555f7cc32e0, b=..., x=..., status=@0x5555fc8934c0: amgx::AMGX_ST_ERROR, xIsZero=false) at .../AMGX/base/src/amg_solver.cu:296
#5  0x00007fff82e93bbe in amgx::(anonymous namespace)::solve_with<(AMGX_Mode)8193, amgx::AMG_Solver, amgx::Vector> (slv=0x5555fc8934a0, rhs=0x555555d02540, sol=0x5555e3b7ff40, resources=0x5555558bd000, xIsZero=false) at .../AMGX/base/src/amgx_c.cu:747
#6  0x00007fff82e3bb03 in AMGX_solver_solve (slv=0x5555fc8934a0, rhs=0x555555d02540, sol=0x5555e3b7ff40) at .../AMGX/base/src/amgx_c.cu:2799
#7  0x0000555555557890 in main (argc=5, argv=0x7fffffffd7f8) at .../AMGX/examples/amgx_capi.c:437
```

`solve_init`
===========

```c++
790│     if (!done)
 791│     {
 792├───────> solve_init(b, x, xIsZero);
 793│     }
 794│
 795│     // Run the iterations
 796│     std::stringstream ss;
 797│
 798│     for (m_curr_iter = 0; m_curr_iter < m_max_iters && !done; ++m_curr_iter)
 799│     {
 800│         // Run one iteration. Compute residual and its norm and decide convergence
 801│         bool has_converged = solve_iteration(b, x, xIsZero);
 802│         // Make sure x is not zero anymore.
.../AMGX/base/src/solvers/solver.cu
```

inside `solve_init`
===================

```c++
71│ template<class T_Config>
372│ void
373│ FGMRES_Solver<T_Config>::solve_init( VVector &b, VVector &x, bool xIsZero )
374│ {
375│     //init residual, even if we don't plan to use it, we might need it, so make sure we have enough memory to store it now
376├───> residual.resize( b.size() );
377│     residual.set_block_dimx( 1 );
378│     residual.set_block_dimy( this->m_A->get_block_dimy() );
379│     residual.dirtybit = 1;
380│     residual.delayed_send = 1;
381│ }
382│
383│ //check for convergence
384│ //al the complicated stuff happens here
385│ template <class TConfig>
386│ bool FGMRES_Solver<TConfig>::checkConvergenceGMRES(bool check_V_0)
.../AMGX/core/src/solvers/fgmres_solver.cu
```

iteration loop
==============

```c++
798│     for (m_curr_iter = 0; m_curr_iter < m_max_iters && !done; ++m_curr_iter)
 799│     {
 800│         // Run one iteration. Compute residual and its norm and decide convergence
 801│         bool has_converged = solve_iteration(b, x, xIsZero);
 802│         // Make sure x is not zero anymore.
 803├───────> xIsZero = false;
 804│         // Is it done ?
 805│         done = m_monitor_convergence && has_converged;
```


