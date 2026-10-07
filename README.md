# Numerical probes of seven results in the `openai/math` collection

The repository [`openai/math`](https://github.com/openai/math) contains 722 mathematical
manuscripts organised into 372 result families, produced by an internal OpenAI model and
released at various stages of verification. Most of the collection is out of reach of a
laptop: the statements are about motives, geometric Langlands, and von Neumann algebras,
and the quantities involved are not computable in any direct sense.

This repository takes seven of the families whose central objects *are* computable and
probes them numerically. For each one it restates the theorem as it appears in the
manuscript, implements the quantity the theorem is about, and compares the computation
against the law the theorem predicts. The figures live in `results/figures/` and the
summary statistics of each run in `results/`.

## What this is and is not

It is a check that the stated asymptotic laws are consistent with the arithmetic at the
scale a single machine can reach, and a set of pictures of those laws. It is not a
verification of any proof. Where a computation agrees with a manuscript, that is evidence
about the statement, not about the argument. Where it disagrees, the disagreement is more
likely to be a finite-range effect than an error in the manuscript, and the two places
where the numbers do not settle are called out below.

The figures are also not decoration on a claim: for several of these results the plotted
quantity is the theorem's own limit, so a figure that looked wrong would be informative.

## The seven results

| # | Family | Statement probed | Figures |
|---|---|---|---|
| 011 | Prime-factor statistics of `p-1` | The normalised logarithms of the prime factors of `p-1` converge to the Poisson--Dirichlet law `PD(1)` | `011-poisson-dirichlet` |
| 012 | Independent largest prime factors of consecutive integers | `P+(n)` and `P+(n+1)` are asymptotically independent, and `P+(n) < P+(n+1)` has density `1/2` | `012-joint-dickman` |
| 021 | A quadratic bound for Jacobsthal's function | `h(k) <= C k^2 / (log log 3k)^2` | `021-jacobsthal` |
| 026 | Positive lower density of large prime gaps | For each fixed `C > 0` a positive proportion of gaps satisfy `d_n > C log p_n` | `026-prime-gaps` |
| 025 | Short Egyptian fractions | `N(b)`, the worst-case minimum length, is `Theta(log log b)` | `025-egyptian-fractions` |
| 017 | The irrationality exponent of `pi` | `mu(pi) = 2` | `017-pi-exponent` |
| 003 | The quasi-Riemann hypothesis | Every Dirichlet `L`-function is zero-free in `Re s > 7/8` | `003-quasi-riemann` |

## 011. The Poisson--Dirichlet law for prime predecessors

**The theorem.** For a prime `p`, list the prime factors of `p - 1` with multiplicity in
decreasing order as `q_1 >= q_2 >= ...`, and set

```
V_j(p) = log q_j(p) / log(p - 1),      q_j(p) = 1 once the list is exhausted.
```

The manuscript proves that for every fixed `k` and every bounded continuous
`F : [0,1]^k -> R`,

```
(1/(pi(x) - 1)) * sum_{3 <= p <= x} F(V_1(p), ..., V_k(p))  ->  E F(L_1, ..., L_k),
```

where `(L_1, L_2, ...)` is the decreasing rearrangement of the stick-breaking sequence
`B_1 = 1 - U_1`, `B_j = (U_1 ... U_{j-1})(1 - U_j)`, whose law is `PD(1)`. This resolves a
conjecture of Ford, Konyagin and Luca.

**What is computed.** Every prime to `3 * 10^7` is factored by sieving the largest prime
factor of each integer up to the limit, giving `1,857,858` primes. The reference law is
handled two ways: `L_1` has the exact Dickman distribution, `Pr[L_1 <= x] = rho(1/x)`, with
`rho` obtained by fourth-order integration of the delay equation `u rho'(u) = -rho(u-1)`;
and the remaining components are drawn by stick breaking with the leading `12` fragments
kept, which drops a tail of mean `2^-12`.

**Figures and numbers.** The left panel of `011-poisson-dirichlet` overlays the empirical
distribution function of `V_j` on the simulated law for `j = 1, 2, 3, 4`; the middle panel
is the density of the first four components; the right panel tracks the mean of `V_1` over
successive ranges of `x`.

The empirical mean of `V_1` rises from `0.5932` in the smallest window to `0.5996` in the
largest, against a `PD(1)` value of `0.6242`. This shortfall is the one place in this
repository where the computation does not obviously reach the predicted law. It is not
noise: with `1.86 * 10^6` samples the standard error is about `0.001`. The regime of the
limit is `x -> infinity`, and at `x = 3 * 10^7` the largest prime factor of `p - 1` is
still systematically smaller than the limiting stick fragment, because the largest part
is governed by the shifted Dickman law `Pr[P+(p-1) <= x^{1/u}] -> rho(u)`, whose approach
is logarithmic in `x`. The mean is still moving upward across the windows in a way
consistent with slow convergence rather than with a different limit, but this repository's
range is not enough to demonstrate that, and the Kolmogorov--Smirnov distances to the
simulated law (reported in `results/011-poisson-dirichlet.json`) are correspondingly
large for `j >= 3`. The first two components agree well.

## 012. Independence of the largest prime factors of consecutive integers

**The theorem.** For fixed `a, b` in `(0, 1)`,

```
(1/X) # {2 <= n <= X : P+(n) <= n^a, P+(n+1) <= n^b}  ->  rho(1/a) rho(1/b),
```

so that `x_n = log P+(n)/log n` and `y_n = log P+(n+1)/log n` have independent limiting
distributions with common distribution function `a -> rho(1/a)`. A corollary is that
`P+(n) < P+(n+1)` has natural density `1/2`.

**What is computed.** `P+(n)` for every `n <= 6 * 10^7` by the same largest-prime-factor
sieve. The comparison density is measured directly, and the joint distribution of
`(x_n, y_n)` is binned on a grid to compare the joint distribution function against the
product of the two marginals.

**Figures and numbers.** In `012-joint-dickman` the left panel is the joint density, the
middle panel contrasts the joint distribution function with the product of marginals along
the diagonal, and the right panel puts the empirical marginal against `rho(1/a)`.

Over `59,999,998` consecutive pairs the density of `P+(n) < P+(n+1)` is `0.49998`. In six
equal windows it is `0.50000, 0.50005, 0.49995, 0.50000, 0.50001, 0.49989`, so there is no
drift in `x`. The supremum distance between the binned joint distribution function and the
product of marginals is `0.0051`, which is the order of the binning bias rather than
evidence of dependence.

## 021. Jacobsthal's function and the quadratic bound

**The theorem.** Let `j(n)` be the least `m` such that every interval of `m` consecutive
integers contains an integer coprime to `n`, and

```
h(k) = sup { j(n) : omega(n) <= k }.
```

Then `h(k) <= C k^2 / (log log 3k)^2` for every `k >= 1`, with `C` absolute. Since
replacing `n` by its radical does not change `j(n)`, `h(k) - 1` is the greatest length of a
run of consecutive integers that can be covered by the divisibility classes of `k` primes.
This answers Jacobsthal's question whether `h(k) << k^2` and improves Iwaniec's
`h(k) << k^2 log^2 k`.

**What is computed.** `h` is a supremum over all sets of at most `k` primes, so only lower
bounds are computable. For a prime set the multiples are periodic with period equal to the
product of the primes, so `j` is found by sieving one period and taking the longest covered
run, plus one. The product of the first `k` primes grows like `e^{k log k}`, so the
computation is capped at `k = 10`, where the period is `6,469,693,230`. The values found,
`2, 4, 6, 10, 14, 22, 26, 34, 40, 46`, reproduce the published sequence for the primorials
`j(P_k)`.

**Figures and numbers.** `021-jacobsthal` plots the computed values against both the new
bound shape and Iwaniec's, and then the ratio of the values to each envelope. Against the
new envelope the implied constant rises to about `0.71` and then flattens; against
Iwaniec's it falls from `0.55` to `0.09`. Both are consistent with the theorems, since a
lower bound for `h` can only ever sit below an upper bound, and the figure should be read
as showing that the quadratic scale is not vacuous at these `k`, not as estimating `C`.

The primorial case is not the extremal case in general: for `k = 24` one has
`j(P_24) = 234 < 236 = h(24)`, so the values here are lower bounds for `h(k)`.

## 026. Positive lower density of large prime gaps

**The theorem.** With `d_n = p_{n+1} - p_n`, for every fixed real `C > 0` there are
`c(C) > 0` and `N_0(C)` such that

```
# {1 <= n <= N : d_n > C log p_n} >= c(C) N      for all N >= N_0(C).
```

The density is measured by counting prime indices, not integers. A corollary is that the
Erdős--Prachar set `{n : p_n/n < p_{n+1}/(n+1)}` has positive lower density.

**What is computed.** All primes to `2 * 10^8`, giving `11,078,937` primes and as many
gaps, and the exceedance proportion for thresholds `C` from `0.25` to `3`, both over the
whole range and within ten contiguous blocks.

**Figures and numbers.** `026-prime-gaps` shows the exceedance curve with the Cramér
heuristic `e^{-C}` for reference, the distribution of `d_n / log p_n`, and the block-by-block
stability of the density.

The proportion of gaps above `C log p_n` is `0.66` at `C = 0.5`, `0.37` at `C = 1`, and
`0.11` at `C = 2`. Across ten blocks these are flat, staying within about `3%` of their
means, which is the content of the positive-density statement that a single count could not
show. The largest gap in the range is `248`, at a ratio `d_n / log p_n` of `13.0`; the
Cramér heuristic lies above the measured curve at large `C`, as expected. The
Erdős--Prachar density is `0.4199`, and the equivalence with `d_n > p_n/n` used to compute
it is exact.

## 025. Short Egyptian fractions

**The theorem.** For integers `1 <= a < b` let `N(a, b)` be the least `k` with
`a/b = 1/n_1 + ... + 1/n_k` for distinct `2 <= n_1 < ... < n_k`, with no bound on the
denominators, and put `N(b) = max_{1 <= a < b} N(a, b)`. Then

```
c_1 log log b <= N(b) <= c_2 log log b,
```

which proves Erdős's conjecture and improves the earlier `log b / log log b` and
`sqrt(log b)` bounds.

**What is computed.** `N(b)` exactly for every `b <= 220`. For each numerator the search is
by iterative deepening on the number of terms with exhaustive search at each depth, using
exact rational arithmetic and two bounds: the remaining terms are at most `1/n` each and at
least `1/(n+t-1)` each. Failure within the cap therefore certifies a larger minimum length,
so these are true minima and `N(b)` is a true maximum. The greedy expansion is computed for
comparison.

**Figures and numbers.** `025-egyptian-fractions` gives the minimum length against the
worst-case greedy length, the growth against the three bound shapes, and the ratio `N(b) /
log log b`.

`N(b)` reaches `6` in this range and is non-decreasing in steps, while the worst greedy
length reaches `9`; the greedy algorithm is not minimal, and `5/121` is a small example
where it needs five terms against the optimal three. The ratio `N(b)/log log b` lies
between `2.4` and `4.8` and drifts slowly downward as `b` grows, which is the behaviour
expected of a bounded `c_2` in the conjectured order. The range `b <= 220` is far too short
to distinguish `log log b` from `sqrt(log b)` or `log b / log log b` by growth alone, which
is why that panel plots all three.

## 017. The irrationality exponent of `pi`

**The theorem.** The irrationality exponent

```
mu(x) = sup { nu > 0 : 0 < |x - p/q| < q^{-nu} for infinitely many p, q in Z }
```

of `pi` is exactly `2`: for every `nu > 2` there is `Q(nu)` with `|pi - p/q| >= q^{-nu}` for
all integers `p` and `q >= Q(nu)`, whether or not the fraction is reduced. Consequently the
Flint--Hills series `sum 1/(n^3 sin^2 n)` converges.

**What is computed.** The continued fraction of `pi` to forty terms at 150 digits, the
convergents as exact integers, and the quantity `q^2 |pi - p/q|` for each. Using Python
integers and arbitrary-precision errors matters here: the convergents reach denominators
above `10^22`, where a float64 representation of `q` destroys the squared quality.

**Figures and numbers.** `017-pi-exponent` shows the squared quality for `pi` alongside
`sqrt(2)` and Catalan's constant, the empirical exponent `-log|pi - p/q| / log q`, and the
partial quotients.

For `pi` the squared quality oscillates in `[0.0034, 1]`, always below `1`, and the
empirical exponent settles at `2.0` over the last ten convergents. The exponent is a
statement about the best approximations, so it is a lower estimate for `mu` and can dip
below `2` at individual convergents; what the definition asks for is the limsup, and the
figure shows the values approaching `2` from above and staying there. `sqrt(2)` behaves
identically, as it must, since its partial quotients are bounded; Catalan's constant does
not, which is the point of the comparison. The largest partial quotient of `pi` in this
range is `292`, at index 4.

## 003. The quasi-Riemann hypothesis

**The theorem.** Every finite-order Hecke `L`-function over `Q(sqrt(-3))` and every
Dirichlet `L`-function, `zeta` included, has no zero in `Re s > 7/8`, the pole at `s = 1`
for a principal character allowed. By the functional equation the nontrivial zeros of
`zeta` therefore satisfy `1/8 <= Re s <= 7/8`, which excludes an exceptional real zero in
`(7/8, 1)`. The theorem does not place the zeros on `Re s = 1/2`.

**What is computed.** The first `600` nontrivial zeros of `zeta`, to height `939`, with
their real parts as actually returned by the computation rather than assumed, and the
count against the Riemann--von Mangoldt smooth count `theta(T)/pi + 1`.

**Figures and numbers.** `003-quasi-riemann` draws the classical de la Vallée Poussin
boundary `1 - c/log t` against the fixed line `7/8`, plots the counting residual
`k - (theta(t_k)/pi + 1)`, and shows the fraction of zeros with real part above `sigma` for
a range of `sigma`.

The two boundaries are qualitatively different objects. The classical region approaches
`Re s = 1` as the height grows and never bounds a zero-free half-plane; the theorem's
boundary is a vertical line cutting the strip at every height, which is what makes it
stronger. All `600` computed zeros have real part `0.5`, so none lies in `Re s > 7/8`: the
number of violations is `0` and the margin is `0.375`. The count of zeros below `T = 939.02`
is `600` against a smooth prediction of `599.72`, a deviation of `0.28`, which confirms the
zero finder is neither missing nor inventing zeros.

This says nothing about zeros off the critical line. `mpmath.zetazero` solves on the line,
so the real parts are `1/2` by construction, and a zero at, say, `Re s = 0.7` would not
appear in the computation at all. The figure is a picture of the claimed region, not a test
of it.

## The animated site

`docs/` is a static GitHub Pages site. It loads one JavaScript file of numerical
series and renders seven animations, each a function of a single progress parameter
so that the same code drives the autoplay, the scrubber, and the first frame.

- **011** the density of `V_1` fills in as primes accumulate, against the exact PD(1)
  density `x · rho(1/x − 1)`; a second panel tracks the six component means.
- **012** the comparison density converges on `1/2` while the independence gap falls.
- **021** `j(P_k)` grows against both bound envelopes, with the implied constant flat
  or falling.
- **026** the exceedance curves for several `C` build up and hold; a second panel
  shows the same densities in successive blocks.
- **025** the greedy expansion of `5/121` walks down a log-scaled remainder staircase
  while the worst-case ratio `N(b)/log log b` drifts.
- **017** `q²|π − p/q|` accumulates for the convergents and the running median of the
  exponent settles on `2`.
- **003** 600 zeros of `zeta` accumulate on `Re s = 1/2`, with the theorem's fixed
  `7/8` boundary and the classical `1 − c/log t` region for contrast.

The page honours `prefers-reduced-motion` by not autoplaying, gives every animation a
play/pause control and a scrubber, and starts each one only when it scrolls into view.

To rebuild the series and refresh the site previews:

```
python scripts/build_site_data.py                 # writes docs/data.js
python scripts/build_site_data.py --only pi,quasi # replaces named series only
node scripts/check_site.js                        # renders every card headlessly
node scripts/render_previews.js                   # PNGs into results/site-previews/
```

`scripts/check_site.js` catches runtime errors in the animation code without a
browser; `scripts/render_previews.js` needs `canvas` from npm and writes the frames
used to check the layout visually.

## Reproducing

```
python -m pip install -e '.[dev]'
python scripts/run_experiments.py                 # all seven probes
python scripts/run_experiments.py --only 011,026  # a subset
python -m pytest                                  # 99 tests
```

The numerical limits are in `config.json` and read through `openmath.config`; the defaults
finish in a few minutes. Each probe writes `results/<id>.json` and a figure under
`results/figures/` as both PNG and PDF.

The tests check the ingredients rather than the conclusions: the Dickman function against
its published values, the `PD(1)` sampler against the exact `L_1` distribution, the
Jacobsthal values against the published primorial sequence, the Egyptian-fraction search
against exhaustive minimality, the `zeta` zeros against the published ordinates, and every
theorem's degenerate cases.

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

The manuscripts are in `openai/math`. The seven probed here are

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
- *The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane `Re s > 7/8`*,
  `preprints/The-Quasi-Riemann-Hypothesis-September-30-2026`

Classical inputs used for the reference values: Dickman (1930) for `rho`; Ford, Konyagin
and Luca (2010) for the prime-predecessor conjecture; Erdős and Pomerance (1978) for the
consecutive-integer problem; Iwaniec (1978) and Vaughan (1977) for Jacobsthal's function;
Erdős (1950) and Vose (1985) for Egyptian fractions.
