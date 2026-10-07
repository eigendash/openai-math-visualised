/* Animations for the seven probes.
 *
 * Each animation is a function of normalised progress p in [0, 1], so the same code
 * serves the auto-playing animation, the scrubber, and the static first frame.  The
 * numerical series come from `docs/data.js`, produced by scripts/build_site_data.py.
 */
'use strict';

const D = window.OM_DATA;
const C = {
  ink: '#1c1c1a', muted: '#6b6b66', grid: '#ececE6', rule: '#e2e2dc',
  accent: '#2b4c7e', accent2: '#b4451f', accent3: '#3c7a5d',
  accent4: '#8a6d1f', accent5: '#6b3f74', dim: '#c9c9c2',
};

/* ---------------------------------------------------------------- utilities */

function ctx2d(canvas) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const rect = canvas.getBoundingClientRect();
  const w = Math.max(rect.width, 10), h = Math.max(rect.height, 10);
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const c = canvas.getContext('2d');
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { c, w, h };
}

function maxAbs(a) { let m = 0; for (const v of a) { const x = Math.abs(v); if (x > m) m = x; } return m; }
function maxOf(a) { let m = -Infinity; for (const v of a) if (v > m) m = v; return m; }
function minOf(a) { let m = Infinity; for (const v of a) if (v < m) m = v; return m; }
function clamp(v, a, b) { return v < a ? a : v > b ? b : v; }

/* Linear interpolation between successive rows of a series, so that animations are
 * smooth across the discrete stages the data was measured at. */
function atStage(rows, p) {
  const m = rows.length - 1;
  const x = clamp(p, 0, 1) * m;
  const i = Math.floor(x), f = x - i;
  const a = rows[i], b = rows[Math.min(i + 1, m)];
  if (Array.isArray(a)) return a.map((v, k) => v + f * (b[k] - v));
  return a + f * (b - a);
}

function stageIndex(n, p) { return Math.min(n - 1, Math.floor(clamp(p, 0, 1) * n)); }

/* Cartesian plotting frame: margins, scales, grid, ticks, labels. */
function frame(c, w, h, o = {}) {
  const m = Object.assign({ l: 54, r: 14, t: 22, b: 38 }, o.margin || {});
  const xlo = o.xlo ?? 0, xhi = o.xhi ?? 1, ylo = o.ylo ?? 0, yhi = o.yhi ?? 1;
  const X = v => m.l + (v - xlo) / (xhi - xlo || 1) * (w - m.l - m.r);
  const Y = v => h - m.b - (v - ylo) / (yhi - ylo || 1) * (h - m.b - m.t);
  c.clearRect(0, 0, w, h);
  if (o.grid !== false) {
    c.strokeStyle = C.grid; c.lineWidth = 1;
    c.beginPath();
    for (const gx of o.xTicks || []) { const px = X(gx); c.moveTo(px, m.t); c.lineTo(px, h - m.b); }
    for (const gy of o.yTicks || []) { const py = Y(gy); c.moveTo(m.l, py); c.lineTo(w - m.r, py); }
    c.stroke();
  }
  c.strokeStyle = C.rule; c.lineWidth = 1;
  c.beginPath(); c.moveTo(m.l, m.t); c.lineTo(m.l, h - m.b); c.lineTo(w - m.r, h - m.b); c.stroke();
  c.fillStyle = C.muted; c.font = '10px ui-monospace, Menlo, monospace';
  c.textAlign = 'center'; c.textBaseline = 'top';
  for (const gx of o.xTicks || []) c.fillText(fmtTick(gx), X(gx), h - m.b + 5);
  c.textAlign = 'right'; c.textBaseline = 'middle';
  for (const gy of o.yTicks || []) c.fillText(fmtTick(gy), m.l - 6, Y(gy));
  c.textAlign = 'left'; c.textBaseline = 'top';
  if (o.xlabel) c.fillText(o.xlabel, m.l, h - 13);
  if (o.ylabel) {
    c.save(); c.translate(11, m.t); c.rotate(-Math.PI / 2); c.textBaseline = 'bottom';
    c.fillText(o.ylabel, 0, 0); c.restore();
  }
  return { X, Y, m, xlo, xhi, ylo, yhi };
}

function fmtTick(v) {
  const a = Math.abs(v);
  if (a === 0) return '0';
  if (a >= 1e6) return (v / 1e6).toFixed(a >= 1e7 ? 0 : 1) + 'M';
  if (a >= 1e3) return (v / 1e3).toFixed(a >= 1e4 ? 0 : 1) + 'k';
  if (a < 0.01) return v.toExponential(0);
  if (a < 1) return v.toFixed(2).replace(/0$/, '');
  if (Number.isInteger(v)) return String(v);
  return v.toFixed(a < 10 ? 1 : 0);
}

function line(c, X, Y, xs, ys, colour, width = 1.8, dash = null) {
  c.save();
  if (dash) c.setLineDash(dash);
  c.strokeStyle = colour; c.lineWidth = width; c.lineJoin = 'round';
  c.beginPath();
  for (let i = 0; i < xs.length; i++) {
    const x = X(xs[i]), y = Y(ys[i]);
    if (i === 0) c.moveTo(x, y); else c.lineTo(x, y);
  }
  c.stroke(); c.restore();
}

/* Only draw the part of a curve up to progress p, so animations reveal as they go. */
function lineUpTo(c, X, Y, xs, ys, p, colour, width = 1.8, dash = null) {
  const n = Math.max(2, Math.round(clamp(p, 0, 1) * (xs.length - 1)) + 1);
  line(c, X, Y, xs.slice(0, n), ys.slice(0, n), colour, width, dash);
}

function dot(c, X, Y, x, y, colour, r = 3) {
  c.fillStyle = colour; c.beginPath(); c.arc(X(x), Y(y), r, 0, 6.2832); c.fill();
}

function label(c, X, Y, x, y, text, colour, dx = 6, dy = -6) {
  c.fillStyle = colour; c.font = '11px ui-monospace, Menlo, monospace';
  c.textAlign = 'left'; c.textBaseline = 'bottom';
  c.fillText(text, X(x) + dx, Y(y) + dy);
}

function legend(c, items, x, y, lh = 13) {
  c.font = '10px ui-monospace, Menlo, monospace';
  c.textAlign = 'left'; c.textBaseline = 'middle';
  items.forEach((it, i) => {
    const yy = y + i * lh;
    c.strokeStyle = it.colour; c.lineWidth = 2;
    if (it.dash) c.setLineDash(it.dash);
    c.beginPath(); c.moveTo(x, yy); c.lineTo(x + 16, yy); c.stroke(); c.setLineDash([]);
    c.fillStyle = C.muted; c.fillText(it.label, x + 21, yy + 0.5);
  });
}

function axisTicks(lo, hi, n) {
  const out = [];
  for (let i = 0; i <= n; i++) out.push(lo + (hi - lo) * i / n);
  return out;
}

/* Band of bars, used for the density histograms. */
function bars(c, X, Y, centres, values, width, colour, alpha = 1) {
  c.save(); c.fillStyle = colour; c.globalAlpha = alpha;
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (!(v > 0)) continue;
    const x0 = X(centres[i] - width / 2), x1 = X(centres[i] + width / 2);
    const y0 = Y(0), y1 = Y(v);
    c.fillRect(x0, Math.min(y0, y1), Math.max(x1 - x0 - 0.6, 0.8), Math.abs(y1 - y0));
  }
  c.restore();
}

/* ------------------------------------------------------- animation driver */

class Animation {
  constructor(canvas, draw, opts = {}) {
    this.canvas = canvas;
    this.draw = draw;
    this.period = opts.period || 12000;
    this.readout = opts.readout || null;
    this.p = 0;
    this.playing = false;
    this.last = 0;
    this._raf = null;
    this._loop = this._loop.bind(this);
    this.render();
  }
  set(p) { this.p = clamp(p, 0, 1); this.render(); }
  play() {
    if (this.playing) return;
    this.playing = true;
    if (this.p >= 1) this.p = 0;
    this.last = performance.now();
    this._raf = requestAnimationFrame(this._loop);
  }
  pause() { this.playing = false; if (this._raf) cancelAnimationFrame(this._raf); this._raf = null; }
  _loop(now) {
    if (!this.playing) return;
    const dt = now - this.last; this.last = now;
    this.p += dt / this.period;
    if (this.p >= 1) { this.p = 1; this.render(); this.playing = false; this.onEnd && this.onEnd(); return; }
    this.render();
    this._raf = requestAnimationFrame(this._loop);
  }
  render() {
    const { c, w, h } = ctx2d(this.canvas);
    this.draw(c, w, h, this.p);
    if (this.readout) this.readout(this.p);
  }
}

/* ------------------------------------------------------------- the cards */

const CARDS = [
  {
    id: 'pd', tag: '011', title: 'Prime factors of p−1 follow the Poisson–Dirichlet law',
    theorem: `List the prime factors of p−1 with multiplicity in decreasing order, q₁ ≥ q₂ ≥ …, and set Vⱼ(p) = log qⱼ(p) / log(p−1). For every fixed k, the joint law of (V₁,…,V_k) over primes p ≤ x converges to the first k components of PD(1), the decreasing rearrangement of the stick-breaking sequence B₁ = 1−U₁, Bⱼ = (U₁⋯U_{j−1})(1−Uⱼ).`,
    eq: 'Vⱼ(p) = log qⱼ(p) / log(p−1)   ⟶   PD(1)',
    captions: [
      ['Density of V₁', 'The histogram fills in as primes accumulate. The curve is the exact law of the largest PD(1) component, Pr[L₁ ≤ x] = ρ(1/x) with ρ the Dickman function.'],
      ['Mean of each component', 'Solid lines: measured over the primes so far. Dashed: the PD(1) value from a 400,000-draw simulation.'],
    ],
    panels: [['density', 'Density of V₁ to 3×10⁷'], ['means', 'Component means']],
    readout(p) {
      const i = stageIndex(D.pd.x.length, p);
      const near = D.pd.means[i], lim = D.pd.simMean;
      const gap = (lim[0] - near[0]) / lim[0] * 100;
      return `primes ≤ <b>${fmtInt(D.pd.x[i])}</b> · mean V₁ = <b>${near[0].toFixed(4)}</b> · ` +
             `PD(1) value <span class="pred">${lim[0].toFixed(4)}</span> · still <b>${gap.toFixed(1)}%</b> below`;
    },
    draw(c, w, h, p, targets) {
      const i = stageIndex(D.pd.x.length, p);
      const dens = D.pd.histSeries[i];
      const cy = atStage(D.pd.means, p);
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const ymax = Math.max(maxOf(D.pd.theoryDensity), maxOf(dens)) * 1.12;
        const F = frame(c1, w1, h1, {
          xlo: 0, xhi: 1, ylo: 0, yhi: ymax,
          xTicks: [0, 0.2, 0.4, 0.6, 0.8, 1], yTicks: axisTicks(0, ymax, 3),
          xlabel: 'V₁', ylabel: 'density',
        });
        bars(c1, F.X, F.Y, D.pd.grid, dens, D.pd.grid[1] - D.pd.grid[0], C.accent, 0.42);
        line(c1, F.X, F.Y, D.pd.grid, D.pd.theoryDensity, C.accent2, 2.2);
        legend(c1, [
          { label: 'PD(1) density, exact', colour: C.accent2 },
          { label: 'primes so far', colour: C.accent },
        ], F.m.l + 10, F.m.t + 10);
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const rows = [D.pd.simMean, ...D.pd.means];
        const ymax = Math.max(maxOf(D.pd.simMean), maxOf(D.pd.means.flat())) * 1.2;
        const F = frame(c2, w2, h2, {
          xlo: 1, xhi: 6, ylo: 0, yhi: ymax,
          xTicks: [1, 2, 3, 4, 5, 6], yTicks: axisTicks(0, ymax, 3),
          xlabel: 'component j', ylabel: 'mean',
        });
        for (let j = 0; j < 6; j++) {
          const col = j % 2 ? C.accent3 : C.accent;
          line(c2, F.X, F.Y, [j + 1, j + 1], [0, D.pd.simMean[j]], C.dim, 6);
          dot(c2, F.X, F.Y, j + 1, cy[j], col, 3.4);
          dot(c2, F.X, F.Y, j + 1, D.pd.simMean[j], C.accent2, 2.4);
        }
        legend(c2, [
          { label: 'PD(1), simulated', colour: C.accent2 },
          { label: 'primes measured', colour: C.accent },
        ], F.m.l + 10, F.m.t + 10);
      }
    },
  },

  {
    id: 'dickman', tag: '012', title: 'Largest prime factors of n and n+1 are independent',
    theorem: `With P⁺(n) the largest prime factor, log P⁺(n)/log n and log P⁺(n+1)/log n are asymptotically independent in natural density, each with distribution function a ↦ ρ(1/a). Hence P⁺(n) < P⁺(n+1) has density exactly 1/2.`,
    eq: 'density{ P⁺(n) < P⁺(n+1) }  ⟶  1/2',
    captions: [
      ['Density of the comparison', 'Measured over all n ≤ 6×10⁷ so far. The dashed line is the predicted 1/2.'],
      ['Joint law against the product of marginals', 'Supremum gap between the binned joint CDF and the product of the two marginals. Zero means independent.'],
    ],
    panels: [['density', 'Density of P⁺(n) < P⁺(n+1)'], ['gap', 'Independence gap']],
    readout(p) {
      const i = stageIndex(D.dickman.x.length, p);
      return `n ≤ <b>${fmtInt(D.dickman.x[i])}</b> · density = <b>${D.dickman.density[i].toFixed(5)}</b> · ` +
             `predicted <span class="pred">0.50000</span> · gap <b>${D.dickman.gap[i].toFixed(4)}</b>`;
    },
    draw(c, w, h, p, targets) {
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const F = frame(c1, w1, h1, {
          xlo: 0, xhi: D.dickman.x.length - 1, ylo: 0.495, yhi: 0.505,
          xTicks: axisTicks(0, D.dickman.x.length - 1, 4).map(v => Math.round(v)),
          yTicks: [0.496, 0.498, 0.5, 0.502, 0.504],
          xlabel: 'window', ylabel: 'density',
        });
        line(c1, F.X, F.Y, [0, D.dickman.x.length - 1], [0.5, 0.5], C.accent2, 1.6, [5, 4]);
        lineUpTo(c1, F.X, F.Y, D.dickman.density.map((_, k) => k), D.dickman.density, p, C.accent, 2.1);
        const i = stageIndex(D.dickman.density.length, p);
        dot(c1, F.X, F.Y, i, D.dickman.density[i], C.accent, 3.2);
        legend(c1, [
          { label: 'measured', colour: C.accent },
          { label: 'predicted 1/2', colour: C.accent2, dash: [5, 4] },
        ], F.m.l + 10, F.m.t + 10);
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const ymax = Math.max(0.05, maxOf(D.dickman.gap) * 1.25);
        const F = frame(c2, w2, h2, {
          xlo: 0, xhi: D.dickman.x.length - 1, ylo: 0, yhi: ymax,
          xTicks: axisTicks(0, D.dickman.x.length - 1, 4).map(v => Math.round(v)),
          yTicks: axisTicks(0, ymax, 3), xlabel: 'window', ylabel: 'sup gap',
        });
        lineUpTo(c2, F.X, F.Y, D.dickman.gap.map((_, k) => k), D.dickman.gap, p, C.accent4, 2.1);
        legend(c2, [{ label: '|joint − product|', colour: C.accent4 }], F.m.l + 10, F.m.t + 10);
      }
    },
  },

  {
    id: 'jacobsthal', tag: '021', title: 'Jacobsthal’s function is at most quadratic in k',
    theorem: `Let j(n) be the least m such that every interval of m consecutive integers contains an integer coprime to n, and h(k) = sup{j(n) : ω(n) ≤ k}. Then h(k) ≤ C·k²/(log log 3k)² with C absolute, answering Jacobsthal's question whether h(k) ≪ k² and improving Iwaniec's k² log²k.`,
    eq: 'h(k)  ≤  C · k² / (log log 3k)²',
    captions: [
      ['Primorial values j(P_k)', 'Computed exactly by sieving one period of the prime product. The curves are the two bound shapes. Only lower bounds for h(k) are computable, so the points must lie under any valid upper bound.'],
      ['Implied constant', 'j(P_k) divided by each envelope. Falling or flat means the bound shape is not contradicted as k grows.'],
    ],
    panels: [['values', 'j(P_k) against the bounds'], ['ratio', 'value ÷ envelope']],
    readout(p) {
      const k = D.jacobsthal.k[stageIndex(D.jacobsthal.k.length, p)];
      const i = Math.max(0, k - 3);
      return `k = <b>${k}</b> · j(P_k) = <b>${D.jacobsthal.j[k - 1]}</b> · ` +
             `new envelope <b>${i < D.jacobsthal.envelope.length ? D.jacobsthal.envelope[i].toFixed(1) : '—'}</b> · ` +
             `ratio <b>${i < D.jacobsthal.ratio.length ? D.jacobsthal.ratio[i].toFixed(3) : '—'}</b>`;
    },
    draw(c, w, h, p, targets) {
      const k = D.jacobsthal.k[stageIndex(D.jacobsthal.k.length, p)];
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const xs = Array.from({ length: 240 }, (_, i) => 1 + i * (D.jacobsthal.kMax - 1) / 239);
        const envAt = (arr, kk) => {
          // the envelope series starts at k = 3
          const idx = kk - 3;
          if (idx < 0) return null;
          const i0 = Math.floor(idx), f = idx - i0;
          if (i0 + 1 >= arr.length) return arr[arr.length - 1];
          return arr[i0] + f * (arr[i0 + 1] - arr[i0]);
        };
        const F = frame(c1, w1, h1, {
          xlo: 1, xhi: D.jacobsthal.kMax, ylo: 0,
          yhi: Math.max(maxOf(D.jacobsthal.iwaniec), maxOf(D.jacobsthal.j)) * 1.25,
          xTicks: axisTicks(1, D.jacobsthal.kMax, 5).map(Math.round),
          yTicks: axisTicks(0, maxOf(D.jacobsthal.iwaniec) * 1.25, 4),
          xlabel: 'k', ylabel: 'value',
        });
        const envLine = (arr, lo) => {
          const px = [], py = [];
          for (const kk of xs) { const v = envAt(arr, kk); if (v != null) { px.push(kk); py.push(v); } }
          line(c1, F.X, F.Y, px, py, lo, 1.7);
        };
        envLine(D.jacobsthal.iwaniec, C.accent3);
        envLine(D.jacobsthal.envelope, C.accent2);
        const upTo = D.jacobsthal.k.slice(0, k);
        line(c1, F.X, F.Y, upTo, D.jacobsthal.j.slice(0, k), C.accent, 2.1);
        for (let i = 0; i < upTo.length; i++) dot(c1, F.X, F.Y, upTo[i], D.jacobsthal.j[i], C.accent, 3);
        dot(c1, F.X, F.Y, k, D.jacobsthal.j[k - 1], C.accent2, 4.2);
        legend(c1, [
          { label: 'j(P_k) computed', colour: C.accent },
          { label: 'k²/(log log 3k)²', colour: C.accent2 },
          { label: 'k² log² k (Iwaniec)', colour: C.accent3 },
        ], F.m.l + 10, F.m.t + 10);
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const ymax = Math.max(maxOf(D.jacobsthal.ratio), maxOf(D.jacobsthal.ratioIwaniec)) * 1.15;
        const F = frame(c2, w2, h2, {
          xlo: 3, xhi: D.jacobsthal.kMax, ylo: 0, yhi: ymax,
          xTicks: axisTicks(3, D.jacobsthal.kMax, 4).map(Math.round),
          yTicks: axisTicks(0, ymax, 3), xlabel: 'k', ylabel: 'value ÷ envelope',
        });
        line(c2, F.X, F.Y, D.jacobsthal.envelopeK, D.jacobsthal.ratio, C.accent2, 2.1);
        line(c2, F.X, F.Y, D.jacobsthal.envelopeK, D.jacobsthal.ratioIwaniec, C.accent3, 1.8, [5, 4]);
        const i = Math.max(0, Math.min(k - 3, D.jacobsthal.ratio.length - 1));
        if (k >= 3) {
          dot(c2, F.X, F.Y, D.jacobsthal.envelopeK[i], D.jacobsthal.ratio[i], C.accent2, 3.4);
          dot(c2, F.X, F.Y, D.jacobsthal.envelopeK[i], D.jacobsthal.ratioIwaniec[i], C.accent3, 3);
        }
        legend(c2, [
          { label: 'vs new bound', colour: C.accent2 },
          { label: 'vs Iwaniec', colour: C.accent3, dash: [5, 4] },
        ], F.m.l + 10, F.m.t + 10);
      }
    },
  },

  {
    id: 'gaps', tag: '026', title: 'A positive proportion of prime gaps exceed C log p',
    theorem: `For every fixed C > 0 there is c(C) > 0 such that, for all sufficiently large N, at least c(C)·N of the indices n ≤ N satisfy p_{n+1} − p_n > C log p_n. The density is over prime indices, not integers.`,
    eq: '#{ n ≤ N : p_{n+1} − p_n > C log p_n }  ≥  c(C)·N',
    captions: [
      ['Exceedance density', 'Proportion of gaps so far exceeding C log p_n, for thresholds C. The dotted curve is the heuristic Poisson model e^{−C}. Each curve has settled to a positive value and stays there.'],
      ['Stability across the range', 'The same densities measured in successive tens of millions of primes. Flatness is the content of a positive density; a single count could not show it.'],
    ],
    panels: [['curve', 'Proportion with d > C log p'], ['stable', 'Density in successive blocks']],
    readout(p) {
      const i = stageIndex(D.gaps.x.length, p);
      const cs = D.gaps.thresholds, cur = D.gaps.curves[i];
      const at = t => cur[cs.findIndex(v => v >= t)];
      return `primes ≤ <b>${fmtInt(D.gaps.limit * D.gaps.x[i] / D.gaps.x[D.gaps.x.length - 1])}</b> · ` +
             `C = 0.5 → <b>${at(0.5).toFixed(3)}</b> · C = 1 → <b>${at(1).toFixed(3)}</b> · ` +
             `C = 2 → <b>${at(2).toFixed(3)}</b> · largest gap <b>${D.gaps.recordGap[i]}</b>`;
    },
    draw(c, w, h, p, targets) {
      const i = stageIndex(D.gaps.x.length, p);
      const cs = D.gaps.thresholds;
      const cur = atStage(D.gaps.curves, p);
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const F = frame(c1, w1, h1, {
          xlo: 0.25, xhi: 3, ylo: 0, yhi: 1,
          xTicks: [0.5, 1, 1.5, 2, 2.5, 3], yTicks: [0, 0.25, 0.5, 0.75, 1],
          xlabel: 'C', ylabel: 'proportion',
        });
        line(c1, F.X, F.Y, cs, D.gaps.cramer, C.dim, 1.6, [3, 3]);
        line(c1, F.X, F.Y, cs, cur, C.accent, 2.3);
        for (const t of [0.5, 1, 2]) {
          const idx = cs.findIndex(v => v >= t);
          if (idx >= 0 && cur[idx] > 0) dot(c1, F.X, F.Y, cs[idx], cur[idx], C.accent2, 3);
        }
        legend(c1, [
          { label: 'measured so far', colour: C.accent },
          { label: 'Cramér e^−C', colour: C.dim, dash: [3, 3] },
        ], F.m.l + 10, F.m.t + 10);
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const ymax = Math.max(0.6, maxOf(D.gaps.curves.flat()) * 1.2);
        const F = frame(c2, w2, h2, {
          xlo: 0, xhi: D.gaps.x.length - 1, ylo: 0, yhi: ymax,
          xTicks: axisTicks(0, D.gaps.x.length - 1, 4).map(v => Math.round(v)),
          yTicks: axisTicks(0, ymax, 3), xlabel: 'block', ylabel: 'density within block',
        });
        const series = [0.5, 1, 2].map(t => cs.findIndex(v => v >= t));
        const colours = [C.accent, C.accent2, C.accent3];
        series.forEach((idx, s) => {
          const ys = D.gaps.curves.map(row => row[idx]);
          const m = D.gaps.x.length;
          const xs = ys.map((_, k) => Math.floor(k * (m - 1) / (ys.length - 1)));
          lineUpTo(c2, F.X, F.Y, xs, ys, p, colours[s], 1.8);
        });
        legend(c2, [
          { label: 'C = 0.5', colour: C.accent },
          { label: 'C = 1', colour: C.accent2 },
          { label: 'C = 2', colour: C.accent3 },
        ], F.m.l + 10, F.m.t + 10);
      }
    },
  },

  {
    id: 'egyptian', tag: '025', title: 'Egyptian fraction expansions need Θ(log log b) terms',
    theorem: `For 1 ≤ a < b let N(a,b) be the least number of distinct unit fractions summing to a/b, and N(b) = max_a N(a,b). Then c₁ log log b ≤ N(b) ≤ c₂ log log b, proving Erdős's conjecture and improving the earlier log b / log log b and √log b bounds.`,
    eq: 'N(b)  =  Θ( log log b )',
    captions: [
      ['Greedy expansion of 5/121', 'The greedy algorithm takes the largest unit fraction that fits. Here it needs five terms; the optimal expansion needs three.'],
      ['N(b)/log log b', 'The worst-case minimum length divided by log log b. A bounded, slowly drifting ratio is what the Θ statement predicts.'],
    ],
    panels: [['worked', 'Greedy vs optimal, 5/121'], ['ratio', 'N(b) ÷ log log b']],
    readout(p) {
      const w = D.egyptian.worked;
      const step = stageIndex(w.greedy.length + 1, p);
      if (step === 0) return `5/121 = ? · greedy is about to start · optimal is <b>${w.optimal.length}</b> terms`;
      const s = w.greedy[step - 1];
      return `5/121 · greedy term <b>${step}</b> of ${w.greedy.length}: 1/${fmtDen(s.den)} · ` +
             `remaining <b>${s.rest === 0 ? '0' : s.rest.toExponential(1)}</b> · ` +
             `optimal <span class="pred">${w.optimal.length} terms</span>`;
    },
    draw(c, w, h, p, targets) {
      const work = D.egyptian.worked;
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const g = work.greedy;
        // Each greedy step subtracts a unit fraction.  The remainder is the story:
        // it falls from 5/121 to exact zero in five steps, whereas three suffice.
        const sup = Math.log10(1.0);            // log10 of the starting value, 5/121
        const rems = [Math.log10(5 / 121), ...g.map(s => s.rest > 0 ? Math.log10(s.rest) : -26)];
        const lo = Math.min(-24, minOf(rems) - 2);
        const frames = g.length + 1;
        const step = stageIndex(frames, p);
        const F = frame(c1, w1, h1, {
          xlo: 0, xhi: frames, ylo: lo, yhi: 0.4,
          xTicks: Array.from({ length: frames }, (_, i) => i + 0.5),
          yTicks: [-24, -18, -12, -6, 0],
          xlabel: 'greedy step', ylabel: 'log₁₀ remainder',
        });
        // the step values, as a descending staircase
        const px = [], py = [];
        for (let i = 0; i <= step; i++) {
          px.push(i + 0.5); py.push(rems[i]);
          if (i < step) { px.push(i + 1.0); py.push(rems[i]); }
        }
        if (px.length > 1) line(c1, F.X, F.Y, px, py, C.accent, 2.2);
        for (let i = 0; i <= step && i < g.length; i++) {
          dot(c1, F.X, F.Y, i + 0.5, rems[i], C.accent, 3.2);
          c1.fillStyle = C.ink; c1.font = '10px ui-monospace, Menlo, monospace';
          c1.textAlign = 'center'; c1.textBaseline = 'bottom';
          c1.fillText('1/' + fmtDen(g[i].den), F.X(i + 0.5), F.Y(rems[i]) - 8);
        }
        if (step >= g.length) {
          dot(c1, F.X, F.Y, g.length + 0.5, rems[g.length], C.accent3, 3.6);
          c1.fillStyle = C.accent3; c1.textAlign = 'right'; c1.textBaseline = 'top';
          c1.fillText('exact', F.X(g.length + 0.5) + 2, F.Y(rems[g.length]) + 6);
          c1.fillStyle = C.accent2; c1.textAlign = 'right'; c1.textBaseline = 'top';
          c1.fillText(`optimal: ${work.optimal.length} terms`, F.X(frames) - 6, F.Y(0.0) + 8);
        }
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const ratio = D.egyptian.ratio;
        const ymax = Math.max(maxOf(ratio) * 1.15, 1);
        const F = frame(c2, w2, h2, {
          xlo: D.egyptian.ratioB[0], xhi: D.egyptian.bMax, ylo: 0, yhi: ymax,
          xTicks: axisTicks(D.egyptian.ratioB[0], D.egyptian.bMax, 4).map(Math.round),
          yTicks: axisTicks(0, ymax, 3), xlabel: 'b', ylabel: 'N(b) ÷ log log b',
        });
        lineUpTo(c2, F.X, F.Y, D.egyptian.ratioB, ratio, p, C.accent, 2.1);
        const i = stageIndex(ratio.length, p);
        dot(c2, F.X, F.Y, D.egyptian.ratioB[i], ratio[i], C.accent2, 3.2);
        legend(c2, [{ label: 'ratio stays bounded', colour: C.accent }], F.m.l + 10, F.m.t + 10);
      }
    },
  },

  {
    id: 'pi', tag: '017', title: 'The irrationality exponent of π is exactly 2',
    theorem: `The irrationality exponent μ(x) = sup{ν : 0 < |x − p/q| < q^{−ν} for infinitely many p/q} satisfies μ(π) = 2: for every ν > 2 there is Q with |π − p/q| ≥ q^{−ν} for all q ≥ Q. Consequently the Flint–Hills series Σ 1/(n³ sin²n) converges.`,
    eq: 'μ(π) = 2',
    captions: [
      ['Squared quality of convergents', 'q²·|π − p/q| for each convergent, computed with exact integer numerators and denominators. The theorem needs this to stay bounded; it never exceeds 1.'],
      ['Running estimate of the exponent', 'The running maximum of −log|π − p/q| / log q. The theorem says the limsup is exactly 2, so this must flatten at 2 rather than keep climbing. √2 behaves the same way; its partial quotients are also bounded.'],
    ],
    panels: [['quality', 'q²·|π − p/q|'], ['exponent', 'running exponent estimate']],
    readout(p) {
      const i = stageIndex(D.pi.exponent.length, p);
      return `convergent <b>${i + 1}</b> of ${D.pi.exponent.length} · q ≈ 10^<b>${D.pi.q[i].toFixed(1)}</b> · ` +
             `q²|π − p/q| = <b>${D.pi.quality[i].toFixed(4)}</b> · ` +
             `running median <b>${D.pi.smoothedExponent[i].toFixed(3)}</b> → <span class="pred">2</span>`;
    },
    draw(c, w, h, p, targets) {
      const i = stageIndex(D.pi.exponent.length, p);
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const F = frame(c1, w1, h1, {
          xlo: 0, xhi: maxOf(D.pi.q) * 1.05, ylo: 0, yhi: 1.05,
          xTicks: [0, 5, 10, 15, 20, 25], yTicks: [0, 0.25, 0.5, 0.75, 1],
          xlabel: 'log₁₀ q', ylabel: 'q²·|π − p/q|',
        });
        line(c1, F.X, F.Y, [0, maxOf(D.pi.q) * 1.05], [1, 1], C.accent2, 1.4, [5, 4]);
        lineUpTo(c1, F.X, F.Y, D.pi.q, D.pi.quality, p, C.accent, 2.0);
        for (let k = 0; k <= i; k++) dot(c1, F.X, F.Y, D.pi.q[k], D.pi.quality[k], C.accent, 2.2);
        dot(c1, F.X, F.Y, D.pi.q[i], D.pi.quality[i], C.accent2, 3.6);
        legend(c1, [
          { label: 'π convergents', colour: C.accent },
          { label: 'bound q²·|π−p/q| = 1', colour: C.accent2, dash: [5, 4] },
        ], F.m.l + 10, F.m.t + 26);
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const qmax = maxOf(D.pi.q) * 1.05;
        const F = frame(c2, w2, h2, {
          xlo: 0, xhi: qmax, ylo: 0, yhi: 3.7,
          xTicks: [0, 5, 10, 15, 20, 25], yTicks: [0, 1, 2, 3],
          xlabel: 'log₁₀ q', ylabel: 'running median',
        });
        line(c2, F.X, F.Y, [0, qmax], [2, 2], C.accent2, 1.6, [5, 4]);
        lineUpTo(c2, F.X, F.Y, D.pi.q, D.pi.smoothedExponent, p, C.accent, 2.2);
        dot(c2, F.X, F.Y, D.pi.q[i], D.pi.smoothedExponent[i], C.accent2, 3.6);
        // the curve starts high and falls, so the legend sits in the empty upper right
        legend(c2, [
          { label: 'running median of the exponent', colour: C.accent },
          { label: 'predicted μ(π) = 2', colour: C.accent2, dash: [5, 4] },
        ], F.X(qmax * 0.42), F.m.t + 12);
      }
    },
  },

  {
    id: 'quasi', tag: '003', title: 'No Dirichlet L-function has a zero in Re s > 7/8',
    theorem: `Every finite-order Hecke L-function over ℚ(√−3) and every Dirichlet L-function, ζ included, is zero-free in the half-plane Re s > 7/8, the pole at s = 1 for a principal character allowed. So the nontrivial zeros of ζ lie in 1/8 ≤ Re s ≤ 7/8, which excludes an exceptional real zero in (7/8, 1).`,
    eq: 'Re s > 7/8  ⟹  L(s, χ) ≠ 0',
    captions: [
      ['Zeros of ζ in the critical strip', 'The first 600 nontrivial zeros. The classical de la Vallée Poussin region 1 − c/log t creeps towards Re s = 1 and bounds nothing at fixed height; the theorem gives the vertical line 7/8, which does.'],
      ['Counting residual', 'k − (θ(t_k)/π + 1), the classical S(T). It stays within ±1 here, confirming the zero finder is not missing or inventing zeros.'],
    ],
    panels: [['strip', 'Zeros and the zero-free boundary'], ['residual', 'Riemann–von Mangoldt residual']],
    readout(p) {
      const i = stageIndex(D.quasi.t.length, p);
      return `zeros plotted <b>${i + 1}</b> of ${D.quasi.count} · height t ≤ <b>${D.quasi.t[i].toFixed(1)}</b> · ` +
             `in Re s > 7/8: <b>${D.quasi.violations}</b> · <span class="pred">theorem predicts 0</span>`;
    },
    draw(c, w, h, p, targets) {
      const n = stageIndex(D.quasi.t.length, p) + 1;
      {
        const { c: c1, w: w1, h: h1 } = ctx2d(targets[0]);
        const F = frame(c1, w1, h1, {
          xlo: 0.45, xhi: 1.0, ylo: 0, yhi: D.quasi.tMax * 1.02,
          xTicks: [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
          yTicks: axisTicks(0, D.quasi.tMax, 4), xlabel: 'Re s', ylabel: 'Im s',
        });
        // the classical zero-free boundary
        line(c1, F.X, F.Y, D.quasi.classical, D.quasi.tGrid, C.accent3, 1.8);
        // the theorem's boundary
        line(c1, F.X, F.Y, [D.quasi.boundary, D.quasi.boundary], [0, D.quasi.tMax], C.accent2, 2.2);
        // zeros
        c1.fillStyle = C.accent;
        for (let i = 0; i < n; i++) {
          c1.beginPath(); c1.arc(F.X(0.5), F.Y(D.quasi.t[i]), 1.5, 0, 6.2832); c1.fill();
        }
        legend(c1, [
          { label: 'zeros, all at Re s = 1/2', colour: C.accent },
          { label: 'Re s = 7/8, zero-free', colour: C.accent2 },
          { label: 'classical 1 − c/log t', colour: C.accent3 },
        ], F.m.l + Math.min(w1 * 0.42, 210), F.m.t + 10);
      }
      {
        const { c: c2, w: w2, h: h2 } = ctx2d(targets[1]);
        const F = frame(c2, w2, h2, {
          xlo: 0, xhi: D.quasi.tMax, ylo: -3, yhi: 3,
          xTicks: axisTicks(0, D.quasi.tMax, 4).map(Math.round),
          yTicks: [-2, -1, 0, 1, 2], xlabel: 'height t', ylabel: 'S(T)',
        });
        line(c2, F.X, F.Y, [0, D.quasi.tMax], [0, 0], C.dim, 1.2, [3, 3]);
        const xs = D.quasi.t.slice(0, n), ys = D.quasi.residual.slice(0, n);
        line(c2, F.X, F.Y, xs, ys, C.accent4, 1.6);
        legend(c2, [{ label: 'k − (θ(t_k)/π + 1)', colour: C.accent4 }], F.m.l + 10, F.m.t + 10);
      }
    },
  },
];

/* ------------------------------------------------------------------ build */

/* Denominators in the Egyptian expansion reach 1.5e24; show them compactly. */
function fmtDen(n) {
  if (n < 1e6) return String(n);
  const e = Math.floor(Math.log10(n));
  const m = n / Math.pow(10, e);
  return `${m.toFixed(m < 10 ? 2 : 1)}e${e}`;
}

function fmtInt(v) {
  if (v >= 1e9) return (v / 1e9).toFixed(2) + '×10⁹';
  if (v >= 1e6) return (v / 1e6).toFixed(v >= 1e7 ? 0 : 1) + '×10⁶';
  if (v >= 1e4) return (v / 1e3).toFixed(0) + 'k';
  return String(Math.round(v));
}

function build() {
  const main = document.getElementById('main');
  const toc = document.getElementById('toc');
  const anims = [];

  for (const card of CARDS) {
    const sec = document.createElement('section');
    sec.className = 'card';
    sec.id = card.id;

    const h2 = document.createElement('h2');
    h2.innerHTML = `<span class="tag">${card.tag}</span><span>${card.title}</span>`;
    sec.appendChild(h2);

    const th = document.createElement('div');
    th.className = 'theorem';
    th.innerHTML = card.theorem + `<span class="eq">${escapeEq(card.eq)}</span>`;
    sec.appendChild(th);

    const ro = document.createElement('p');
    ro.className = 'readout';
    sec.appendChild(ro);

    const panels = document.createElement('div');
    panels.className = 'panels two';
    const canvases = [];
    card.panels.forEach(([key, cap], idx) => {
      const wrap = document.createElement('div');
      wrap.className = 'panel';
      const cv = document.createElement('canvas');
      cv.setAttribute('aria-label', cap);
      wrap.appendChild(cv);
      const p = document.createElement('p');
      p.className = 'cap';
      p.innerHTML = `<b>${cap}.</b> ${card.captions[idx][1]}`;
      wrap.appendChild(p);
      panels.appendChild(wrap);
      canvases.push(cv);
    });
    sec.appendChild(panels);

    const controls = document.createElement('div');
    controls.className = 'controls';
    const play = document.createElement('button');
    play.type = 'button';
    play.textContent = 'Play';
    const reset = document.createElement('button');
    reset.type = 'button';
    reset.textContent = 'Reset';
    const scrubWrap = document.createElement('div');
    scrubWrap.className = 'scrub-wrap';
    const scrub = document.createElement('input');
    scrub.type = 'range'; scrub.min = '0'; scrub.max = '1000'; scrub.value = '0';
    scrub.className = 'scrub';
    scrub.setAttribute('aria-label', `Scrub the ${card.title} animation`);
    const stage = document.createElement('span');
    stage.className = 'stage';
    scrubWrap.appendChild(scrub); scrubWrap.appendChild(stage);
    controls.appendChild(play); controls.appendChild(reset); controls.appendChild(scrubWrap);
    sec.appendChild(controls);

    const updateStage = p => {
      stage.textContent = `progress ${(p * 100).toFixed(0)}%`;
      scrub.value = String(Math.round(p * 1000));
    };

    const anim = new Animation(
      canvases[0],
      (c, w, h, p) => card.draw(c, w, h, p, canvases),
      {
        period: 14000,
        readout: p => { ro.innerHTML = card.readout(p); updateStage(p); },
      }
    );
    anim.onEnd = () => { play.textContent = 'Replay'; play.setAttribute('aria-pressed', 'false'); };
    anims.push(anim);

    play.addEventListener('click', () => {
      if (anim.playing) { anim.pause(); play.textContent = 'Play'; play.setAttribute('aria-pressed', 'false'); }
      else { anim.play(); play.textContent = 'Pause'; play.setAttribute('aria-pressed', 'true'); }
    });
    reset.addEventListener('click', () => {
      anim.pause(); anim.set(0); play.textContent = 'Play'; play.setAttribute('aria-pressed', 'false');
    });
    scrub.addEventListener('input', () => {
      anim.pause(); play.textContent = 'Play'; play.setAttribute('aria-pressed', 'false');
      anim.set(Number(scrub.value) / 1000);
    });

    main.appendChild(sec);
    const link = document.createElement('a');
    link.href = '#' + card.id;
    link.textContent = `${card.tag} · ${card.short || card.title.split(' ').slice(0, 4).join(' ')}`;
    toc.appendChild(link);
  }

  // Redraw on resize so canvases stay crisp.
  let rt = null;
  window.addEventListener('resize', () => {
    clearTimeout(rt);
    rt = setTimeout(() => anims.forEach(a => a.render()), 120);
  });

  // Autoplay only when the reader has not asked for reduced motion, and only when
  // the card is on screen, so the page is not all motion at once.
  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  if (!reduced && 'IntersectionObserver' in window) {
    const seen = new Set();
    const io = new IntersectionObserver(entries => {
      for (const e of entries) {
        if (e.isIntersecting && !seen.has(e.target)) {
          seen.add(e.target);
          const idx = CARDS.findIndex(cd => cd.id === e.target.id);
          if (idx >= 0) {
            const a = anims[idx];
            a.play();
            const btn = e.target.querySelector('button');
            btn.textContent = 'Pause'; btn.setAttribute('aria-pressed', 'true');
          }
        }
      }
    }, { threshold: 0.35 });
    CARDS.forEach(cd => io.observe(document.getElementById(cd.id)));
  }

  // Keep the animation independent of scroll for readers who want to inspect.
  document.addEventListener('visibilitychange', () => {
    if (document.hidden) anims.forEach(a => a.pause());
  });
}

function escapeEq(s) {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', build);
else build();
