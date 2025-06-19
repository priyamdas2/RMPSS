# RMPSS: Recursive Modified Pattern Search on Simplex

This repository provides MATLAB implementations of the Recursive Modified Pattern Search on the Simplex (RMPSS) algorithm.

---

## 📚 Citation

If you use this method in your research, please cite:

> Das, Priyam (2021).  
> *Recursive Modified Pattern Search on High-Dimensional Simplex: A Blackbox Optimization Technique*.  
> Sankhya B, 83 (Suppl 2), 440–483.  
> https://doi.org/10.1007/s13571-020-00236-9

---

## Files

- **`RMPSS.m`**  
  Implements Recursive Modified Pattern Search on the simplex **without parallel threading**.

- **`RMPSS_parallel.m`**  
  Implements Recursive Modified Pattern Search on the simplex **with parallel threading**.

---

## Notes

In this particular example, the objective function is **not computationally expensive**.  
As a result, the parallel version (`RMPSS_parallel.m`) **may take more time** than the single-threaded version due to the overhead of parallelization.




