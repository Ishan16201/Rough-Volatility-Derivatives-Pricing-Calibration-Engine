# Mathematical foundations and conventions

This note records the equations implemented in RVRE. All time inputs are
year-fractions. Rates are continuously compounded annualized rates; variance
is annualized; Brownian increments have units \(\sqrt{\text{year}}\).

## Brownian drivers and correlation

For a grid interval \(\Delta t_i\), the simulator draws independent standard
normals \(Z_i,Z_i^\perp\) and sets
\[
\Delta W_i=\sqrt{\Delta t_i}Z_i,\qquad
\Delta W_i^S=\sqrt{\Delta t_i}
(\rho Z_i+\sqrt{1-\rho^2}Z_i^\perp).
\]
Thus each increment has variance \(\Delta t_i\) and covariance
\(\rho\Delta t_i\). The correlation input is checked against its model
domain. Explicit independent streams prevent shared mutable RNG state.

## The specified Volterra process and its covariance

The implemented rough driver is the Riemann--Liouville (RL) process
\[
Y_t=\sqrt{2H}\int_0^t(t-s)^{H-\frac12}\,dW_s,\quad 0<H<\frac12.
\]
The exponent exceeds \(-1/2\), so the kernel's square is integrable. Itô
isometry gives
\[
\operatorname{Var}(Y_t)=2H\int_0^t(t-s)^{2H-1}\,ds=t^{2H}.
\]
For \(s,t\ge0\), its covariance is
\[
C_{\rm RL}(t,s)=2H\int_0^{\min(t,s)}
(t-u)^{H-\frac12}(s-u)^{H-\frac12}\,du.
\]

The sometimes-quoted expression
\[
\tfrac12\left(t^{2H}+s^{2H}-|t-s|^{2H}\right)
\]
is the covariance of standard fractional Brownian motion, not this RL process
for general \(H\). Accordingly, the implementation validates the RL integral
off the diagonal and \(t^{2H}\) on the diagonal. It does not claim both
covariance formulas hold for one process. A future fBm-based specification
would be a separate model/discretization decision.

## Discrete Volterra driver

For grid \(0=t_0<\dots<t_N\), the reference scheme projects the kernel onto
Brownian increments using cell-average coefficients
\[
\bar K_{i,j}=\frac{\sqrt{2H}}{(H+\frac12)\Delta t_j}
\left((t_i-t_j)^{H+\frac12}-(t_i-t_{j+1})^{H+\frac12}\right),
\quad j<i,
\]
\[
Y_i^{\rm ref}=\sum_{j<i}\bar K_{i,j}\Delta W_j,\qquad Y_0=0.
\]
The approximation variance is
\(\nu_i^{\rm ref}=\sum_{j<i}\bar K_{i,j}^2\Delta t_j\), which need not equal
\(t_i^{2H}\) at finite resolution. The hybrid candidate uses midpoint weights
on older cells and a conditionally exact diagonal Brownian integral. It
matches the diagonal integral's marginal variance and covariance with the
last Brownian increment, but does not factorize all bridge cross-covariances
exactly.

## Rough Bergomi discrete model

The model is
\[
V_t=\xi_0(t)\exp\left(\eta Y_t-\tfrac12\eta^2t^{2H}\right),\qquad
\frac{dS_t}{S_t}=(r_t-q_t)dt+\sqrt{V_t}\,dW_t^S,
\quad d\langle W^S,W\rangle_t=\rho dt.
\]
For the discretized driver \(Y_i\), the code uses its actual marginal
\(\nu_i=\operatorname{Var}(Y_i)\) and sets
\[
V_i=\xi_0(t_i)\exp(\eta Y_i-\tfrac12\eta^2\nu_i).
\]
This yields \(\mathbb E[V_i]=\xi_0(t_i)\) under that discrete Gaussian
marginal. It does not remove joint-covariance or asset-payoff discretization
bias.

Let \(D_r(0,t),D_q(0,t)\) be positive normalized discount curves. The
left-variance log step is
\[
\log\frac{S_{i+1}}{S_i}
=\log\frac{D_r(0,t_i)}{D_r(0,t_{i+1})}
-\log\frac{D_q(0,t_i)}{D_q(0,t_{i+1})}
-\tfrac12V_i\Delta t_i+\sqrt{V_i}\Delta W_i^S.
\]
Conditionally on time-\(t_i\) information, \(V_i\) is fixed and the Gaussian
exponential has expectation one after the variance drift term. Hence the
discrete conditional forward relation is preserved under the stated
correlated-driver construction. Empirical finite-sample diagnostics do not
equal exact identities.

## Heston and option payoffs

Heston follows
\[
\frac{dS_t}{S_t}=(r-q)dt+\sqrt{v_t}\,dW_t^S,\quad
dv_t=\kappa(\theta-v_t)dt+\sigma_v\sqrt{v_t}\,dW_t^v,\quad
d\langle W^S,W^v\rangle_t=\rho dt.
\]
Its Monte Carlo variance scheme stores a raw Euler state and uses
\(v^+=\max(v,0)\) in drift, diffusion and stock variance; the next state is
not projected. This is full-truncation Euler, not exact CIR simulation.
Characteristic-function prices target the continuous-time model instead.

For a European call or put, the discounted estimator is
\[
\widehat C=D_r(0,T)\frac1N\sum_{n=1}^N(S_T^{(n)}-K)^+,\qquad
\widehat P=D_r(0,T)\frac1N\sum_{n=1}^N(K-S_T^{(n)})^+.
\]
Confidence intervals quantify sampling error under the specified estimator;
they do not quantify grid bias, model error or market-data uncertainty.
