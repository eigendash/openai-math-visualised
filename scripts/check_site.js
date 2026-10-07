// Minimal DOM + canvas stubs so app.js can be executed outside a browser.
const fs = require('fs');
const path = require('path');

function makeCtx() {
  const noop = () => {};
  const ctx = {
    canvas: null,
    setTransform: noop, clearRect: noop, beginPath: noop, moveTo: noop, lineTo: noop,
    stroke: noop, fill: noop, arc: noop, fillRect: noop, strokeRect: noop,
    save: noop, restore: noop, translate: noop, rotate: noop, setLineDash: noop,
    fillText: noop, measureText: () => ({ width: 10 }),
    globalAlpha: 1, fillStyle: '', strokeStyle: '', lineWidth: 1, font: '',
    textAlign: '', textBaseline: '', lineJoin: '',
  };
  return ctx;
}

const listeners = {};
const elements = {};
function makeEl(tag) {
  const el = {
    tagName: tag, children: [], style: {}, dataset: {}, _id: '',
    className: '', type: '', value: '0', min: '0', max: '1000',
    get id() { return this._id; }, set id(v) { this._id = v; elements[v] = this; },
    textContent: '', innerHTML: '',
    setAttribute(k, v) { this[k] = v; },
    getAttribute(k) { return this[k]; },
    appendChild(c) { this.children.push(c); return c; },
    addEventListener(ev, fn) { (listeners[ev] = listeners[ev] || []).push(fn); },
    getBoundingClientRect() { return { width: 480, height: 240 }; },
    getContext() { return makeCtx(); },
    querySelector() { return makeEl('button'); },
    width: 0, height: 0,
  };
  return el;
}

global.window = {
  devicePixelRatio: 1,
  matchMedia: () => ({ matches: true }),   // reduced motion: no autoplay, no rAF
  addEventListener: () => {},
  IntersectionObserver: undefined,
};
global.document = {
  readyState: 'complete',
  hidden: false,
  getElementById: id => elements[id] || (elements[id] = Object.assign(makeEl('div'), { _id: id })),
  createElement: makeEl,
  addEventListener: (ev, fn) => { (listeners[ev] = listeners[ev] || []).push(fn); },
  querySelector: () => null,
};
global.performance = { now: () => Date.now() };
global.requestAnimationFrame = () => 0;
global.cancelAnimationFrame = () => {};
global.clearTimeout = clearTimeout;
global.setTimeout = setTimeout;

// load the data and the app, as the browser would
const ROOT = path.resolve(__dirname, '..');
const docs = path.join(ROOT, 'docs');
eval(fs.readFileSync(path.join(docs, 'data.js'), 'utf8'));
const src = fs.readFileSync(path.join(docs, 'app.js'), 'utf8');
eval(src + '\n;globalThis.__CARDS = CARDS;');

// exercise every card at several progress values
const results = [];
const main = document.getElementById('main');
console.log('cards built:', main.children.length);
for (const card of globalThis.__CARDS) {
  const canvases = [];
  const el = document.getElementById(card.id);
  // find canvases by walking the element tree
  (function walk(n) { for (const c of n.children) { if (c.tagName === 'canvas') canvases.push(c); walk(c); } })(el);
  if (canvases.length !== 2) { results.push(`${card.id}: expected 2 canvases, got ${canvases.length}`); continue; }
  for (const p of [0, 0.13, 0.5, 0.87, 1]) {
    try {
      card.draw(makeCtx(), 480, 240, p, canvases);
      card.readout(p);
    } catch (e) {
      results.push(`${card.id} @ p=${p}: ${e.message}`);
    }
  }
}
const newData = fs.readdirSync(docs).length;
console.log(newData ? '' : '');
if (results.length === 0) console.log('ALL CARDS OK at p = 0, 0.13, 0.5, 0.87, 1');
else { console.log('FAILURES:'); results.forEach(r => console.log('  ' + r)); process.exitCode = 1; }
