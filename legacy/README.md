# Legacy (2020 PHYD57 CUDA code)

The original, as-submitted CUDA C source from the 2020 PHYD57 project
*"Galaxy Collisions With CUDA and FFT"*. **Preserved for reference only** — it is
not built or run by the modernized package and contains the known bugs documented
in [`../AGENT.md`](../AGENT.md) §3.6.

| File | Role |
|---|---|
| `CUDAfft2.0.cu` | FFT-Poisson solver test harness (no particles, no time loop) |
| `final_draft1.cu` | Full simulation, primary (OpenMP pragmas mostly inert) |
| `Multi-Parallel-CUDA-FFT.cu` | Full simulation, later OpenMP variant |

Original build (CUDA 9.2 + cuFFT + DISLIN, hard-coded 2020 lab paths):

```
nvcc final_draft1.cu -Xcompiler -fopenmp \
  -I/home/phyd57/N_Body1/9.2/include -L/home/phyd57/N_Body1/9.2/lib64 -lcufft \
  -o CUDAfftcu2.out -I/usr/local/dislin -ldislin
```
