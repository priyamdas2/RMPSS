# RMPSS

**Recursive Modified Pattern Search on the Simplex (RMPSS)** is a C++-based optimization routine implemented for high-dimensional, non-convex objective functions defined on the probability simplex. This package provides a fast, derivative-free optimization method accessible in R via `Rcpp`.

### 📚 Citation

If you use this method in your research, please cite:

> Das, Priyam (2021).  
> *Recursive Modified Pattern Search on High-Dimensional Simplex: A Blackbox Optimization Technique*.  
> Sankhya B, 83 (Suppl 2), 440–483.  
> https://doi.org/10.1007/s13571-020-00236-9

## 🔧 Installation

You can install the development version of the package directly from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("priyamdas2/RMPSS")
```

## ⚙️ Function Overview

The main function in this package is:

**`RMPSS_opt()`**

### Usage

```r
RMPSS_opt(x0, func,
          s_init = 1,
          no_runs = 1000,
          max_iter = 10000,
          rho_1 = 2,
          rho_2 = 2,
          phi = 1e-20,
          lambda = 1e-10,
          tol_fun = 1e-6,
          tol_fun_2 = 1e-20,
          print_output = 0)
```

# Example 1: Ackley-type function
```r
g <- function(y)
  return(-20 * exp(-0.2 * sqrt(0.5 * (y[1] ^ 2 + (y[2]-1) ^ 2))) 
         - exp(0.5 * (cos(2 * pi * y[1]) + cos(2 * pi * (y[2]-1)))) 
         + exp(1) + 20)

# Global min value is 0, achieved at c(0, 1)
starting_point <- c(0.4, 0.6)
g(starting_point)
solution <- RMPSS_opt(starting_point, g)
g(solution)
```
# Example 2: Infeasible starting point (not on simplex)

```r
g <- function(y) return(-y[1])  # Min value is 1 if y[1] = 1
RMPSS_opt(c(1, 0.2, 56, 0.4), g)  # Starting point is invalid; resets internally
```

# Example 3: 1000-dimensional problem

```r
g <- function(y) return(-sum(y^10))
# Min value is -1, achieved if only one coordinate is 1, others are 0
RMPSS_opt(rep(1/1000, 1000), g, print_output = 1)
```
