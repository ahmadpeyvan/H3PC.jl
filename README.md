# H3PC.jl — Hypersonic, High-Order, High-Performance Code in Julia

H3PC.jl is a high-order computational framework for simulating compressible Euler and Navier–Stokes equations, with a strong emphasis on **hypersonic**, **high-speed**, and **multi-physics** flow regimes.  
It is built on top of the high-order solver **Trixi.jl**, extending its capabilities with:

- adaptive mesh refinement (AMR) specialized for hypersonics  
- advanced viscous and inviscid flux formulations  
- axisymmetric extensions  
- high-Mach-number stability improvements  
- Mutation++ real-gas thermochemistry integration  
- production-scale MPI parallelization for HPC systems  

> **Note**  
> The H3PC.jl **source code is not publicly visible** in this repository.  
> To obtain access, please contact **Ahmad Peyvan** or **Prof. George Em Karniadakis** (see “Access to Source Code” below).

---

## 🔗 Base Code — Trixi.jl

H3PC.jl is built upon and extends the functionality of:

**Trixi.jl**  
https://github.com/trixi-framework/Trixi.jl

We gratefully acknowledge the Trixi development team for providing a powerful and flexible foundation for high-order hyperbolic PDE solvers in Julia.

---

## 🔗 Real-Gas Chemistry — Mutation.jl

H3PC.jl supports **real-gas thermochemistry** through **Mutation.jl**, which wraps the Mutation++ Fortran library into Julia.

We acknowledge **Khemraj Shukla** for creating Mutation.jl and enabling robust thermochemical modeling directly in Julia.

**Mutation.jl:**  

---

# 🚀 Running H3PC.jl

This repository contains documentation only. Once you are granted access to the full source code, you may run H3PC.jl as described below.

---

## 🧪 1. Running H3PC.jl in Serial Mode

```bash
julia --project H3PC.jl run_examples/your_case.jl
```

Step-by-step:

1. Activate project:
   ```julia
   using Pkg
   Pkg.activate(".")
   ```
2. Instantiate dependencies:
   ```julia
   Pkg.instantiate()
   ```
3. Run any example:
   ```bash
   julia --project run_examples/hypersonic_plate.jl
   ```

For better JIT performance:
```bash
julia --compile=all --threads=auto --project run_examples/...jl
```

---

## 🧵 2. Running with Julia Multithreading

To use all available CPU threads:

```bash
JULIA_NUM_THREADS=auto julia --project run_examples/...jl
```

Or manually specify:

```bash
JULIA_NUM_THREADS=16 julia --project run_examples/...jl
```

H3PC.jl parallelizes reconstruction, flux computations, and AMR operations.

---

## 🖥️ 3. Running in MPI Parallel Mode (Distributed Memory)

H3PC.jl supports MPI-based distributed memory parallelism via `MPI.jl`.

### Example: 4 MPI ranks

```bash
mpiexec -n 4 julia --project run_examples/...jl
```

Hybrid mode (MPI + threads):

```bash
export JULIA_NUM_THREADS=8
mpiexec -n 4 julia --project run_examples/...jl
```

### Slurm HPC example

```bash
srun -n 64 julia --project run_examples/capsule_3D.jl
```

Verify MPI installation:

```julia
using MPI
MPI.version()
```

---

# 🔒 Access to the H3PC.jl Source Code

The full H3PC.jl source is **not included** in this public repository.

To request access:

### Contact:
**Ahmad Peyvan**  
Email: ahmadpeyvan@gmail.com
Division of Applied Mathematics, Brown University

or

**Prof. George Em Karniadakis**
Email: george_karniadakis@brown.edu  
Division of Applied Mathematics, Brown University

Please include:

- Full name & affiliation  
- Intended use (research, academic, benchmarking, etc.)  
- Confirmation of **non-commercial** use  

Access is granted on a case-by-case basis.

---

# 📝 Citation

If you use H3PC.jl in your research, please cite:

> A. Peyvan, G. E. Karniadakis, et al., **H3PC.jl: High-Order Hypersonics & Hyperbolic PDE Code** (2025).  
> Additional references will be provided with source access.

---

# 🙏 Acknowledgments

- The **Trixi.jl** team for the base solver architecture.  
- **Khemraj Shukla** for Mutation.jl and the Julia Mutation++ wrapper.  
- Collaborators and supporting institutions at Brown University DAM.
