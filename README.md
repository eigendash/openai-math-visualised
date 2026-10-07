# Numerical probes of seven results in the `openai/math` collection

The repository [`openai/math`](https://github.com/openai/math) contains 722
manuscripts organised into 372 result families, produced by an internal OpenAI model
and released at varying stages of verification. The bulk of the collection concerns
objects for which no direct numerical evaluation is available. The seven families
treated here are those, among the 372 surveyed, whose central quantities can be
computed directly.

For each family we state the theorem in the form given by the manuscript, identify
the quantity the theorem constrains, implement it, and compare the result against the
predicted law. Figures are in `results/figures/` and per-run statistics in
`results/`. The animated presentation of the same computations is in `docs/`.

## Status of the computations

A numerical agreement with a stated asymptotic law is evidence about the statement
and not about the proof; a disagreement admits two readings, of which an error in the
manuscript is the less likely. In six of the seven cases below the computed quantity
approaches the predicted value over the accessible range. In one case, family 011,
the approach is measurably incomplete at the largest limit we could reach, and we
quantify the shortfall rather than treat it as noise. One computation, family 003,
does not test the theorem at all, for reasons given in that section.

## Summary

| Family | Quantity computed | Predicted | Result |
|---|---|---|---|
| 011 | `V_j(p)` for primes `p <= 3*10^7` | `PD(1)` | components 1–2 agree; 3–8 deviate, see text |
| 012 | density of `P+(n) < P+(n+1)`, `n <= 6*10^7` | `1/2` | `0.499983` |
| 021 | `j(P_k)`, `k <= 10` | `<= C k^2/(log log 3k)^2` | ratio rises to `0.71`, then flat |
| 026 | density of `d_n > C log p_n`, `p <= 2*10^8` | positive for each `C` | `0.60, 0.34, 0.11` at `C = 0.5, 1, 2` |
| 025 | `N(b)`, `b <= 220` | `Theta(log log b)` | `N(b)/log log b` in `[2.43, 4.80]` |
| 017 | `q^2 abs(pi - p/q)`, 40 convergents | bounded, `mu(pi) = 2` | max `0.94`, exponent settles at `2.0` |
| 003 | 600 zeros of `zeta`, with real parts | no zero in `Re s > 7/8` | 0 violations; claim not tested |

---

## 011. Prime-factor statistics of `p - 1`

**Theorem.** For a prime `p`, let `q_1(p) >= q_2(p) >= ...` be the prime factors of
`p - 1` listed with multiplicity, and set `q_j(p) = 1` once the list is exhausted.
Put `V_j(p) = log q_j(p) / log(p - 1)`. Then for every fixed `k` and every bounded
continuous `F: [0,1]^k -> R`,

```
(1/(pi(x) - 1)) * sum_{3 <= p <= x} F(V_1(p), ..., V_k(p))  ->  E F(L_1, ..., L_k),
```

where `(L_1, L_2, ...)` is the decreasing rearrangement of the stick-breaking sequence
`B_1 = 1 - U_1`, `B_j = (U_1 ... U_{j-1})(1 - U_j)`, whose law is `PD(1)`. This
resolves a conjecture of Ford, Konyagin and Luca.

**Computation.** The largest prime factor of every integer up to `3*10^7` is obtained
by sieving, which factors `p - 1` for all `1,857,858` primes in the range. The
reference law is generated in two ways: `L_1` has the exact distribution
`Pr[L_1 <= x] = rho(1/x)`, with `rho` obtained by fourth-order integration of the delay
equation `u rho'(u) = -rho(u-1)`; the remaining components are drawn by stick breaking
with the leading twelve fragments retained, a tail of mean `2^-12`.

**Results.** Write `e_j` for the empirical mean of `V_j` over all primes in range and
`lambda_j` for the corresponding `PD(1)` mean from `4*10^5` draws.

| `j` | `e_j` | `lambda_j` | `e_j / lambda_j` | KS distance |
|---:|---:|---:|---:|---:|
| 1 | 0.5976 | 0.6245 | 0.957 | 0.045 |
| 2 | 0.1968 | 0.2094 | 0.940 | 0.075 |
| 3 | 0.0891 | 0.0882 | 1.010 | 0.265 |
| 4 | 0.0503 | 0.0403 | 1.248 | 0.426 |
| 5 | 0.0301 | 0.0190 | 1.584 | 0.433 |
| 6 | 0.0173 | 0.0090 | 1.922 | 0.430 |

Two features require comment. The first two components agree with `PD(1)` to within a
few per cent, and `e_1` increases monotonically across six equal windows of the prime
range, from `0.5932` to `0.5996`, so the discrepancy has the sign and the monotonicity
of an unconverged limit rather than of a different one. The components from `j = 3`
onward, however, exceed their predicted means by factors that grow with `j`, and the
Kolmogorov--Smirnov distances against the simulated law are correspondingly large. The
mass accounted for by the first eight components is `0.99542`, against `0.99607` for
`PD(1)`.

The deficit is consistent with the known slow convergence of the largest-part law:
the theorem's regime is `x -> infinity`, and the largest prime factor of `p - 1` is
governed by a shifted Dickman law whose approach is logarithmic in `x`. The behaviour
we observe is that the empirical partition of logarithmic mass is flatter than
`PD(1)` at this scale, carrying more mass in intermediate components and less in the
leading one. We do not claim this identifies the mechanism, and the range
`x <= 3*10^7` is in any case too short to establish the rate. What can be said is that
the first two components are consistent with the theorem while the remainder are not,
and that no attempt was made to extend the computation far enough to decide whether
the deviation decays.

## 012. Independence of the largest prime factors of consecutive integers

**Theorem.** For fixed `a, b` in `(0, 1)`,

```
(1/X) # {2 <= n <= X : P+(n) <= n^a, P+(n+1) <= n^b}  ->  rho(1/a) rho(1/b),
```

so that `x_n = log P+(n)/log n` and `y_n = log P+(n+1)/log n` have independent
limiting distributions with common distribution function `a -> rho(1/a)`. It follows
that `P+(n) < P+(n+1)` has natural density `1/2`.

**Computation.** `P+(n)` for all `n <= 6*10^7`, by the same largest-prime-factor sieve.
The comparison density is evaluated directly, and the joint distribution of
`(x_n, y_n)` is binned on a uniform grid to compare the joint distribution function
against the product of the two marginals.

**Results.** Over `59,999,998` consecutive pairs the density of `P+(n) < P+(n+1)` is
`0.499983`. Evaluated in six equal windows the values are `0.50000`, `0.50005`,
`0.49995`, `0.50000`, `0.50001`, `0.49989`; there is no monotone drift. The supremum
distance between the binned joint distribution function and the product of marginals
is `0.0051`, which is of the order of the binning bias and does not decrease further
with sample size over the range examined. The empirical marginal of `x_n` lies below
`rho(1/a)` for `a` near `1`, by `0.039` at `a = 0.8` and `a = 0.9`; this is the same
direction of finite-range effect as in family 011, and we note it rather than
attribute it.

## 021. A quadratic bound for Jacobsthal's function

**Theorem.** Let `j(n)` be the least `m` such that every interval of `m` consecutive
integers contains an integer coprime to `n`, and `h(k) = sup { j(n) : omega(n) <= k }`.
Then `h(k) <= C k^2 / (log log 3k)^2` for every `k >= 1`, with `C` absolute. Since
replacing `n` by its radical does not change `j(n)`, `h(k) - 1` is the greatest length
of a run of consecutive integers covered by the divisibility classes of `k` primes.
The bound answers Jacobsthal's question whether `h(k) << k^2` and improves Iwaniec's
`h(k) << k^2 log^2 k`.

**Computation.** `h` is a supremum over all sets of at most `k` primes, so only lower
bounds are accessible. For a fixed prime set the multiples are periodic with period
equal to the product of the primes, and `j` is the longest covered run in one period,
plus one. That product grows as `e^{k log k}`, which limits the computation to `k = 10`,
where the period is `6,469,693,230`. The values obtained,

```
j(P_k) = 2, 4, 6, 10, 14, 22, 26, 34, 40, 46    (k = 1, ..., 10),
```

reproduce the published primorial sequence.

**Results.** Writing `E(k) = k^2/(log log 3k)^2`, the ratios `j(P_k)/E(k)` for
`k = 3, ..., 10` are

```
0.413, 0.518, 0.556, 0.688, 0.658, 0.710, 0.702, 0.689,
```

and the corresponding ratios against Iwaniec's envelope `k^2 log^2 k` fall from
`0.552` to `0.087`. Both sequences are consistent with the respective theorems, a lower
bound for `h` being necessarily below an upper bound. The first sequence rises and then
flattens near `0.70`; we do not interpret this as an estimate of `C`, since the range in
`k` is short and the primorial case is not extremal in general. On the latter point, the
manuscript notes that `j(P_24) = 234 < 236 = h(24)`, so the computed values bound `h(k)`
strictly from below.

## 026. Positive lower density of large prime gaps

**Theorem.** With `d_n = p_{n+1} - p_n`, for every fixed real `C > 0` there are
`c(C) > 0` and `N_0(C)` such that

```
# {1 <= n <= N : d_n > C log p_n} >= c(C) N      for all N >= N_0(C).
```

The density is taken over prime indices. Consequently the set
`{n : p_n/n < p_{n+1}/(n+1)}` has positive lower density; the equivalence with
`d_n > p_n/n` is exact.

**Computation.** All primes to `2*10^8`, giving `11,078,937` gaps, and the exceedance
proportion at thresholds `C` in `[0.25, 3]`, both over the full range and within ten
contiguous blocks.

**Results.** The proportion of gaps exceeding `C log p_n` is `0.602` at `C = 0.5`,
`0.338` at `C = 1` and `0.110` at `C = 2`. Across the ten blocks each of these is
constant to within about `3` per cent of its mean; for `C = 1` the block values are
`0.380, 0.394, 0.408, 0.382, 0.355, 0.360, 0.364, 0.367, 0.370, 0.373`, which shows the
stabilisation that a single count over the whole range could not. The largest gap in
range is `248`, at `d_n / log p_n = 13.0`. The heuristic Poisson reference `e^{-C}` lies
above the measured curve for `C > 1.6`. The Erdős--Prachar density is `0.4199`.

## 025. Short Egyptian fractions

**Theorem.** For integers `1 <= a < b` let `N(a, b)` be the least `k` for which
`a/b = 1/n_1 + ... + 1/n_k` with `2 <= n_1 < ... < n_k` distinct and no bound on the
denominators, and set `N(b) = max_{1 <= a < b} N(a, b)`. Then

```
c_1 log log b <= N(b) <= c_2 log log b,
```

which proves Erdős's conjecture and improves the earlier bounds `log b / log log b`
and `sqrt(log b)`.

**Computation.** `N(b)` exactly for every `b <= 220`. For each numerator the expansion
is sought by iterative deepening on the number of terms, exhaustive at each depth, in
exact rational arithmetic, using two bounds: at a node with remainder `r`, `t` terms
remaining and least usable denominator `n`, one has `t/n >= r` and
`1/(n + t - 1) <= r`. A search that fails within the depth cap therefore certifies a
larger minimum length, so the tabulated values are exact. The greedy expansion is
computed alongside for comparison.

**Results.** `N(b)` attains `6` on this range while the worst-case greedy length attains
`9`; the greedy algorithm is not minimal. The least example in the range is `a/b = 5/121`,
where greedy requires five terms and three suffice. The ratio `N(b)/log log b` lies in
`[2.43, 4.80]` and drifts slowly downward over the range, which is the behaviour a
bounded `c_2` predicts. The range `b <= 220` does not separate `log log b` from
`sqrt(log b)` or `log b / log log b` by growth alone, and the figure plots all three
envelopes for that reason.

## 017. The irrationality exponent of `pi`

**Theorem.** The irrationality exponent of an irrational real `x` is

```
mu(x) = sup { nu > 0 : 0 < |x - p/q| < q^{-nu} for infinitely many p, q in Z }.
```

The manuscript proves `mu(pi) = 2`: for every `nu > 2` there is `Q(nu)` such that
`|pi - p/q| >= q^{-nu}` for all integers `p` and `q >= Q(nu)`, the fraction need not be
in lowest terms. Consequently the Flint--Hills series `sum 1/(n^3 sin^2 n)` converges.

**Computation.** The continued fraction of `pi` to forty terms at 150 digits, with
convergents held as exact integers and the errors as arbitrary-precision values. The
representation matters: the convergents reach `q > 6*10^22`, where a float64
representation of `q` destroys the quantity `q^2 |pi - p/q|` entirely.

**Results.** The quantities `q^2 |pi - p/q|` at the convergents have maximum `0.95` and
minimum `0.0034`, so the sequence is bounded above by `1` as the theorem requires,
and the empirical exponent `-log|pi - p/q| / log q` settles at `2.0` over the final ten
convergents. The exponent is defined by a supremum over infinitely many approximants, so
this is a lower estimate and individual convergents take values below `2`; what the
figure shows is the running median of the exponents converging to `2` from above and
remaining there. The control `sqrt(2)`, whose partial quotients are also bounded,
behaves identically. Catalan's constant does not, which is the purpose of the
comparison. The largest partial quotient of `pi` in range is `292`, at index `4`.

## 003. The quasi-Riemann hypothesis

**Theorem.** Every finite-order Hecke `L`-function over `Q(sqrt(-3))` and every
Dirichlet `L`-function, `zeta` included, has no zero in `Re s > 7/8`, the pole at
`s = 1` for a principal character excepted. By the functional equation the nontrivial
zeros of `zeta` then satisfy `1/8 <= Re s <= 7/8`, so an exceptional real zero in
`(7/8, 1)` is excluded. The theorem does not place zeros on `Re s = 1/2`.

**Computation.** The first `600` nontrivial zeros of `zeta`, to height `939.02`, with
real parts recorded as returned rather than assumed, and the count compared with the
Riemann--von Mangoldt smooth count `theta(T)/pi + 1`.

**Results.** All `600` zeros have real part `0.5`; the number lying in `Re s > 7/8` is
`0`, and the margin is `0.375`. The count below `T = 939.02` is `600` against a smooth
prediction of `599.72`, a deviation of `0.28`.

These numbers do not test the theorem. The zeros are located on the critical line, so
their real parts are `1/2` by construction and a zero at, say, `Re s = 0.7` would not
appear in the output at all. The computation verifies only that the zero finder is
consistent with the count formula, and the figure should be read as a drawing of the
claimed region. This limitation is intrinsic to the method rather than to the
implementation, and it is stated here because a reader could otherwise take the empty
intersection for confirmation.

---

## Reproducing

```
python -m pip install -e '.[dev]'
python scripts/run_experiments.py                 # all seven probes
python scripts/run_experiments.py --only 011,026  # a subset
python -m pytest                                  # 99 tests
```

Numerical limits are in `config.json`, read through `openmath.config`. The defaults
complete in a few minutes on one machine. Each probe writes `results/<id>.json` and a
figure under `results/figures/` in both PNG and PDF.

The test suite checks the inputs rather than the conclusions. It covers the Dickman
function against its published values, the `PD(1)` sampler against the exact law of
`L_1`, the Jacobsthal values against the published primorial sequence, the Egyptian
search against exhaustive minimality, the `zeta` zeros against the published ordinates,
and the degenerate cases of each implemented identity.

## The animated site

`docs/` is a static site deployed to GitHub Pages. It loads one file of numerical
series and renders seven animations, each expressed as a function of a single progress
parameter so that the same code drives the autoplay, the scrubber, and the static first
frame. Each card states its theorem, gives a live numerical readout, and captions the
two panels.

The page suppresses autoplay under `prefers-reduced-motion`, provides a play/pause
control and a scrubber for every animation, and begins each one only when it enters the
viewport. To rebuild the series and refresh the previews:

```
python scripts/build_site_data.py                 # writes docs/data.js
python scripts/build_site_data.py --only pi,quasi # replaces the named series only
node scripts/check_site.js                        # renders every card headlessly
node scripts/render_previews.js                   # PNGs into results/site-previews/
```

`scripts/check_site.js` executes the animation code without a browser and catches
runtime errors; `scripts/render_previews.js` requires the `canvas` package from npm and
writes the frames used to inspect layout.

## Layout

```
src/openmath/
  sieves.py               prime, smallest-prime-factor and largest-prime-factor sieves
  poisson_dirichlet.py    PD(1) sampling, the Dickman function, predecessor factoring
  joint_dickman.py        P+(n) coordinates, joint and marginal distribution functions
  jacobsthal.py           periods, longest covered runs, the two bound envelopes
  prime_gaps.py           exceedance density, block stability, the Prachar comparison
  egyptian.py             greedy and shortest unit-fraction expansions
  pi_exponent.py          continued fractions, exact convergents, empirical exponent
  quasi_riemann.py        zeta zeros, zero-free regions, zero counts
  plotting.py             shared figure style
  config.py               numerical limits
scripts/
  build_site_data.py      compact animation series for the site
  check_site.js           headless smoke test of the animation code
  render_previews.js      render animation frames to PNG for review
  run_experiments.py      the seven probes, figures and JSON summaries
docs/                     the static GitHub Pages site
tests/
results/                  JSON summaries, figures, and site previews
```

## References

The probed manuscripts, with their directories in `openai/math`:

- *The Poisson--Dirichlet Law for Prime Predecessors*,
  `preprints/The-Poisson-Dirichlet-Law-for-Prime-Predecessors-September-24-2026`
- *The joint Dickman law for consecutive integers*,
  `preprints/The-joint-Dickman-law-for-consecutive-integers-September-24-2026`
- *A quadratic bound for Jacobsthal's function*,
  `preprints/A-quadratic-bound-for-Jacobsthals-function-September-25-2026`
- *Positive lower density of large prime gaps*,
  `preprints/Positive-lower-density-of-large-prime-gaps-September-25-2026`
- *Short Egyptian fractions*, `preprints/Short-Egyptian-fractions-September-25-2026`
- *The irrationality exponent of pi is 2*,
  `preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026`
- *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane* `Re s > 7/8`,
  `preprints/The-Quasi-Riemann-Hypothesis-September-30-2026`

Reference values and prior bounds are taken from Dickman (1930) for `rho`; Ford,
Konyagin and Luca (2010) for the prime-predecessor conjecture; Erdős and Pomerance
(1978) for the consecutive-integer problem; Iwaniec (1978) and Vaughan (1977) for
Jacobsthal's function; Erdős (1950) and Vose (1985) for Egyptian fractions.
