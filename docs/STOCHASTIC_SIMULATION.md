# Stochastic simulation layer

This document describes the random-number, Brownian and Volterra layer used by
the Rough Bergomi model. Pricing and calibration are covered in the dedicated
documents referenced from the repository README.

## Units and grids

All time points are in years. A Brownian increment over an interval of length
\(\Delta t_i\) is a normal variate with variance \(\Delta t_i\) and units
year\(^{1/2}\). A `TimeGrid` starts at zero and is strictly increasing; it may
be uniform or nonuniform. At least one positive-length interval is required.

Path arrays are contiguous and path-major: all intervals/points for one path
are adjacent. Generated paths use stable per-path substreams, so a path's
random values do not depend on execution ordering or thread scheduling.

## RNG contract

The generator is Philox4x32-10 with the Random123 round constants and a
64-bit seed/key, 64-bit stream identifier, and 64-bit block counter. The
uniform mapping uses 52 random bits and returns strictly interior values in
\((0,1)\). Standard normals use the Box–Muller transform with one cached
normal. Philox's known-answer vector is tested. For a fixed build/runtime,
identical seed, stream and draw sequence replay exactly. Transcendental
functions in Box–Muller may differ across platforms; cross-platform bitwise
normal equality is not promised.

`RandomStream` is explicit value-owned state; there is no global RNG. A caller
selects seed and stream ID. `substream(id)` deterministically derives another
stream. The hybrid method additionally takes an explicit bridge RNG stream and
domain-separates its per-path bridge draws from Brownian path streams.

## Brownian simulation

For each path and interval, generate
\[
\Delta W_i=\sqrt{\Delta t_i}\,Z_i,\qquad Z_i\sim N(0,1),
\]
with independent standard normals. The cumulative path is \(W_0=0\) and
\(W_{i+1}=W_i+\Delta W_i\). For two processes with target correlation
\(\rho\in[-1,1]\), use independent standard normals \(Z_i,Z_i^\perp\):
\[
\Delta W_i^1=\sqrt{\Delta t_i}Z_i,\qquad
\Delta W_i^2=\sqrt{\Delta t_i}
(\rho Z_i+\sqrt{(1-\rho)(1+\rho)}Z_i^\perp).
\]
The factored expression for \(1-\rho^2\) is used for numerical stability.
The Brownian module validates the finite correlation domain and rejects
zero paths, malformed arrays, and invalid grids.

## Volterra reference discretization

The implemented process is the Riemann–Liouville (RL) Volterra process:
\[
Y_t=\sqrt{2H}\int_0^t(t-s)^{H-\frac12}\,dW_s,\qquad 0<H<\frac12.
\]
At zero, define \(Y_0=0\) directly; do not evaluate the singular kernel at
the endpoint. The singularity has exponent \(H-\tfrac12>-1/2\), so its square
is integrable and the Itô integral is well-defined. No kernel clipping is used.

On \(0=t_0<\cdots<t_N\), approximate the kernel by its cell average on each
Brownian interval:
\[
\bar K_{i,j}=\frac{\sqrt{2H}}{t_{j+1}-t_j}
\int_{t_j}^{t_{j+1}}(t_i-s)^{H-\frac12}\,ds
=\frac{\sqrt{2H}}{(H+\frac12)(t_{j+1}-t_j)}
\left[(t_i-t_j)^{H+\frac12}
-(t_i-t_{j+1})^{H+\frac12}\right],\quad j<i,
\]
\[
Y_{t_i}^{\rm ref}=\sum_{j=0}^{i-1}\bar K_{i,j}\Delta W_j,\qquad
Y_{t_0}^{\rm ref}=0.
\]
This is the orthogonal cell-average approximation of the deterministic kernel
using Brownian increments, not the exact stochastic integral on each cell.
The continuous limit has variance \(t^{2H}\); the finite-grid reference
variance is \(\sum_{j<i}\bar K_{i,j}^2\Delta t_j\) and generally differs from
that limit. Uniform-grid lag coefficients are precomputed in \(O(N)\) memory.
The path loop is \(O(PN^2)\) arithmetic and retains \(O(PN)\) input/output
values for \(P\) paths and \(N\) intervals. Nonuniform grids use the explicit
cell endpoints in the same formula.

Coefficient differences are evaluated with `log1p`/`expm1` when intervals are
small relative to their distance from the singular endpoint, avoiding
catastrophic cancellation. The final cell has a direct non-cancelling
expression. No time-zero kernel is evaluated.

## Hybrid comparison scheme

`hybrid_exact_diagonal` treats the most recent Brownian interval separately.
For \(h=\Delta t_{i-1}\), its exact diagonal integral \(I_i\) has
\[
\operatorname{Var}(I_i)=h^{2H},\qquad
\operatorname{Cov}(I_i,\Delta W_{i-1})
=\frac{\sqrt{2H}}{H+\frac12}h^{H+\frac12}.
\]
Conditioning on the interval increment gives
\[
I_i=\frac{\sqrt{2H}}{H+\frac12}h^{H-\frac12}\Delta W_{i-1}
+\sqrt{h^{2H}\left(1-\frac{2H}{(H+\frac12)^2}\right)}Z_i^*,
\]
where \(Z_i^*\) is independent standard normal. Older intervals use kernel
midpoint weights. This removes the cell-average approximation at the singular
diagonal and matches its variance and covariance with that same interval's
Brownian increment.

The residual factor is evaluated as
\(h^H(\frac12-H)/(H+\frac12)\), using
\((H+\frac12)^2-2H=(\frac12-H)^2\) rather than subtracting nearly equal
floating-point numbers near \(H=\frac12\).

The per-output residual normals are not a full joint factorization of all
Brownian-bridge integrals: their cross-covariance with later outputs' older
interval contributions is approximated. Therefore the hybrid method is a
comparison candidate, not an exact process sampler. Both schemes have
\(O(PN^2)\) path work and \(O(PN)\) retained path storage in this initial
implementation; the hybrid also draws \(O(PN)\) bridge normals. No asymptotic
speed-up is presumed. The benchmark reports measured runtime and the
discretization variance error as grid size changes.

## Covariance and an important specification inconsistency

For the RL kernel above, the covariance is
\[
C_{\rm RL}(t,s)=2H\int_0^{\min(t,s)}
(t-u)^{H-\frac12}(s-u)^{H-\frac12}\,du.
\]
In particular \(C_{\rm RL}(t,t)=t^{2H}\), as required. The covariance
\[
\frac{t^{2H}+s^{2H}-|t-s|^{2H}}{2}
\]
is the covariance of standard fractional Brownian motion, not this RL process
for \(H\ne\frac12\). The two requested formulas cannot both describe the
same process. This implementation follows the explicitly specified
zero-start Volterra integral, uses its covariance integral for off-diagonal
tests, and tests the diagonal identity \(t^{2H}\). It does not claim that the
standard-fBm covariance identity holds. The covariance helper evaluates the
integral by fixed-order Gauss–Legendre quadrature after a variable transform
that regularizes the integrable endpoint singularity.

## Statistical validation and benchmarks

Brownian tests estimate sample means, variances, step cross-covariance, and
correlation. Test bounds are multiples of the corresponding Gaussian
large-sample standard errors, with fixed seeds; they are not arbitrary fixed
decimal tolerances.

Volterra tests compare empirical variances/covariances with the RL covariance,
including an explicit discretization allowance, validate \(Y_0=0\), exercise
nonuniform grids, and verify decreasing deterministic terminal variance
error for both schemes on refined grids. Sampling confidence and discretization
error are kept separate.

Configure with `-DRVRE_BUILD_BENCHMARKS=ON` to build
`rvre_volterra_benchmark`. It compares the methods on common Brownian input at
quarter, half, and full requested time resolutions (duplicate levels are
omitted), and prints seed, path/grid/H, runtime, terminal empirical variance
error and its sampling standard error, empirical covariance and covariance
error/standard error at the sampled times, an approximate kernel-term count,
and estimated retained double-payload bytes. Example:

```sh
cmake --preset release -DRVRE_BUILD_BENCHMARKS=ON
cmake --build --preset release --target rvre_volterra_benchmark
./build/release/rvre_volterra_benchmark 1000 64 0.1
```

The memory estimate includes Brownian input, Volterra output and coefficient
arrays, but excludes allocator metadata, executable/code, and RNG state. Time
and memory are observations for the specified machine/configuration, not
portable performance claims.
