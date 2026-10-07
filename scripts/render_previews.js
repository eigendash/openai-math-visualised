/* Render each card's canvases to PNG at several progress values, for visual review. */
const fs = require('fs');
const path = require('path');
const { createCanvas } = require('/tmp/canvacheck/node_modules/canvas');

function realCtx(canvas) {
  const ctx = canvas.getContext('2d');
  const noop = () => {};
  // node-canvas has everything except these two, which we add for the stub contract
  if (!ctx.setTransform) ctx.setTransform = noop;
  return ctx;
}

const elements = {};
function makeEl(tag) {
  const el = {
    tagName: tag, children: [], style: {}, _id: '', className: '', type: '',
    value: '0', min: '0', max: '1000', textContent: '', innerHTML: '',
    get id() { return this._id; }, set id(v) { this._id = v; elements[v] = this; },
    setAttribute(k, v) { this[k] = v; }, getAttribute(k) { return this[k]; },
    appendChild(c) { this.children.push(c); return c; },
    addEventListener() {},
    getBoundingClientRect() { return { width: el._w || 520, height: el._h || 250 }; },
    querySelector() { return makeEl('button'); },
    width: 0, height: 0,
  };
  if (tag === 'canvas') {
    const cv = createCanvas(el._w || 520, el._h || 250);
    el.getContext = t => cv.getContext(t);
    Object.defineProperty(el, 'width', { get: () => cv.width, set: v => { cv.width = v; } });
    Object.defineProperty(el, 'height', { get: () => cv.height, set: v => { cv.height = v; } });
    el.toBuffer = (t) => cv.toBuffer(t);
    el._canvas = cv;
  }
  return el;
}

global.window = { devicePixelRatio: 2, matchMedia: () => ({ matches: true }), addEventListener: () => {} };
global.document = {
  readyState: 'complete', hidden: false,
  getElementById: id => elements[id] || (elements[id] = Object.assign(makeEl('div'), { _id: id })),
  createElement: makeEl, addEventListener: () => {},
};
global.performance = { now: () => Date.now() };
global.requestAnimationFrame = () => 0;
global.cancelAnimationFrame = () => {};

const ROOT = path.resolve(__dirname, '..');
const docs = path.join(ROOT, 'docs');
eval(fs.readFileSync(path.join(docs, 'data.js'), 'utf8'));
eval(fs.readFileSync(path.join(docs, 'app.js'), 'utf8') + '\n;globalThis.__CARDS = CARDS;');

console.log('HARNESS data keys:', Object.keys(window.OM_DATA || {}).join(', '));
const outDir = path.join(ROOT, 'results', 'site-previews');
fs.mkdirSync(outDir, { recursive: true });

// Rebuild each card's canvases at a chosen size, then draw.
function renderCard(card, p, w, h) {
  const canvases = card._canvases;
  for (const cv of canvases) { cv._w = w; cv._h = h; }
  const mocks = canvases.map(cv => {
    // reset the backing store to the requested size
    cv._canvas.width = w * 2; cv._canvas.height = h * 2;
    return { getContext: () => cv._canvas.getContext('2d'), width: cv._canvas.width, height: cv._canvas.height,
             getBoundingClientRect: () => ({ width: w, height: h }) };
  });
  try {
    card.draw(Object.create(null), w, h, p, mocks);
  } catch (e) { console.error('DRAW FAIL', card.id, p, e.stack.split('\n').slice(0,4).join('\n')); throw e; }
  return canvases;
}

const PROGRESS = [0.12, 0.55, 1.0];
for (const card of globalThis.__CARDS) {
  // collect the canvases created during build for this card
  const el = elements[card.id];
  const canvases = [];
  (function walk(n) { for (const c of n.children) { if (c.tagName === 'canvas') canvases.push(c); walk(c); } })(el);
  card._canvases = canvases;
  for (const p of PROGRESS) {
    const mocks = renderCard(card, p, 520, 250);
    mocks.forEach((m, i) => {
      const buf = canvases[i]._canvas.toBuffer('image/png');
      fs.writeFileSync(path.join(outDir, `${card.id}-p${String(Math.round(p*100)).padStart(3,'0')}-${i}.png`), buf);
    });
    // readout text for the same stage
    console.log(`${card.id.padEnd(11)} p=${p.toFixed(2)}  ${card.readout(p).replace(/<[^>]+>/g, '')}`);
  }
}
console.log('\npreviews written to results/site-previews/');
