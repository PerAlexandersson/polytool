const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const { test } = require('node:test');

const html = fs.readFileSync(process.env.POLYTOOL_HTML || path.join(__dirname, 'index.html'), 'utf8');

test('comparison selection changes exported formula and preserves checked domain', () => {
  const start = html.indexOf('function renderRecurrenceMenu(');
  const end = html.indexOf('function createWorker(', start);
  assert.ok(start >= 0 && end > start);
  let displayed = '';
  let copied = '';
  const slot = {};
  const candidate = (lag, derivative, first) => ({
    recurrence: `candidate-${lag}-${derivative}`, latex: 'P(n)',
    recurrence_json: {schema: 'polytool.recurrence.v1', initial_polynomials: Array(first - 1).fill(['1'])},
    complexity: {rec_len: lag, diff_deg: derivative, var_deg: 2, idx_deg: 1,
      denom_var_deg: 0, denom_idx_deg: 1},
    first_checked_index: first, last_checked_index: 20,
    fit_polynomials: 10, verification_polynomials: 10, singular_checked_indices: [],
  });
  const report = {
    representations: [candidate(3, 0, 4), candidate(1, 2, 2)],
    frontier: [0, 1], status: 'budget_exhausted', source_rows: 20,
    configurations_considered: 10, search: {max_var_deg: 3, max_idx_deg: 1},
    coverage: [{rec_len: 3, diff_deg: 0, status: 'found', configurations_considered: 5},
      {rec_len: 1, diff_deg: 1, status: 'incomplete', configurations_considered: 5}],
  };
  const sandbox = vm.createContext({
    window: {_recurrenceMenu: report},
    document: {getElementById: () => slot},
    navigator: {clipboard: {writeText: async text => { copied = text; }}},
    esc: value => String(value).replaceAll('<', '&lt;'),
    setResult: value => { displayed = value; },
    renderRecurrenceExtension: n => `<extension rows="${n}">`,
  });
  vm.runInContext(html.slice(start, end), sandbox);
  sandbox.renderRecurrenceMenu(0);
  assert.match(displayed, /partial comparison/);
  assert.match(displayed, /not proofs for all n/);
  assert.match(displayed, /incomplete/);
  assert.match(displayed, /extension rows="20"/);
  assert.equal(JSON.parse(sandbox.window._lastRecJson).initial_polynomials.length, 3);
  sandbox.renderRecurrenceMenu(1);
  assert.equal(sandbox.window._lastRecPlain, 'candidate-1-2');
  assert.equal(JSON.parse(sandbox.window._lastRecJson).initial_polynomials.length, 1);
  sandbox.copyRecurrenceMenu({});
  assert.equal(JSON.parse(copied).status, 'budget_exhausted');
  report.representations[1].singular_checked_indices = [5];
  sandbox.renderRecurrenceMenu(1);
  assert.match(displayed, /undefined at checked indices 5/);
  assert.doesNotMatch(displayed, /<extension/);
  report.representations = [];
  report.frontier = [];
  sandbox.renderRecurrenceMenu(0);
  assert.match(displayed, /No candidate returned/);
});

test('actual WASM returns a menu and its selection generates exact rows', {
  skip: !process.env.POLYTOOL_WASM_DIR,
}, async () => {
  const dir = process.env.POLYTOOL_WASM_DIR;
  const sandbox = vm.createContext({TextDecoder, TextEncoder, WebAssembly});
  vm.runInContext(fs.readFileSync(path.join(dir, 'polytool_web.js'), 'utf8'), sandbox);
  const wasm = vm.runInContext('wasm_bindgen', sandbox);
  await wasm(fs.readFileSync(path.join(dir, 'polytool_web_bg.wasm')));
  const huge = 9007199254740993n;
  const rows = Array.from({length: 8}, (_, n) => String(huge * 2n ** BigInt(n))).join('\n');
  const menu = JSON.parse(wasm.find_recurrence_representations(rows, 2, 0, 0, 0, false, false, false, 100, 0, 0));
  assert.equal(menu.schema, 'polytool.recurrence-representations.v1');
  assert.equal(menu.status, 'complete');
  assert.ok(menu.representations.length > 0);
  const chosen = menu.representations[menu.frontier[0]];
  assert.equal(chosen.evidence, 'exact_fit');
  assert.equal(chosen.recurrence_json.initial_polynomials[0][0], String(huge));
  const extended = JSON.parse(wasm.generate_recurrence_rows(JSON.stringify(chosen.recurrence_json), 12));
  assert.equal(extended.polynomials[11][0], String(huge * 2n ** 11n));
  const partial = JSON.parse(wasm.find_recurrence_representations(rows, 2, 0, 0, 0, false, false, false, 0, 0, 0));
  assert.equal(partial.status, 'budget_exhausted');
  const bad = JSON.parse(wasm.find_recurrence_representations('1\n2\n???\n4', 2, 0, 0, 0, false, false, false, 10, 0, 0));
  assert.match(bad.error, /invalid polynomial input/);
  // Exercise the actual browser example and the denominator-zero-constant chart.
  const sampleStart = html.indexOf('  A171608: {');
  const sampleEnd = html.indexOf('  A001263: {', sampleStart);
  const sample = vm.runInNewContext('({' + html.slice(sampleStart, sampleEnd) + '})').A171608;
  assert.equal(sample.lines.length, 41);
  const tradeoffs = JSON.parse(wasm.find_recurrence_representations(
    sample.lines.slice(1).join('\n'), 3, 3, 1, 2, false, true, false, 2000, 0, 1));
  assert.equal(tradeoffs.status, 'complete');
  for (const [r, d] of [[3, 0], [2, 1], [1, 2]]) {
    assert.ok(tradeoffs.frontier.some(i =>
      tradeoffs.representations[i].complexity.rec_len === r &&
      tradeoffs.representations[i].complexity.diff_deg === d));
  }
});
