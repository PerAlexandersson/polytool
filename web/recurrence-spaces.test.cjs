const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const {test} = require('node:test');
const html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8');

test('space display separates evidence and exports the complete certificate', async () => {
  let displayed = '', copied = '';
  const sandbox = vm.createContext({
    window: {}, setResult: text => {displayed = text;},
    navigator: {clipboard: {writeText: async text => {copied = text;}}},
    esc: value => String(value).replaceAll('<', '&lt;'),
  });
  const start = html.indexOf('function renderRecurrenceSpaces(');
  vm.runInContext(html.slice(start, html.indexOf('function createWorker(', start)), sandbox);
  const report = {request: {domain_start: 3}, frontier: [0], results: [
    {status: 'successor_exists', evidence: 'conditional_on_source', certificate: {basis: [['1']]},
      successor_projection_rank: 1, residual_count: 4,
      representatives: [{formula: '<unsafe>', singular_ranks: ['51'], requires_x_division: true}]},
    {status: 'successor_exists', evidence: 'exact_fit', representatives: []},
    {status: 'absent_in_shape', evidence: 'rows_only_obstruction', representatives: []},
    {status: 'incomplete_budget', evidence: 'none', reason: 'budget exhausted', representatives: []},
  ]};
  sandbox.renderRecurrenceSpaces(report);
  assert.match(displayed, /domain starts at n = 3/);
  assert.match(displayed, /conditional on the supplied source/);
  assert.match(displayed, /no all-rank conclusion/);
  assert.match(displayed, /no successor recurrence in this shape/);
  assert.match(displayed, /No mathematical conclusion/);
  assert.match(displayed, /Singular base ranks: 51/);
  assert.match(displayed, /&lt;unsafe>/);
  assert.doesNotMatch(displayed, /<unsafe>/);
  await sandbox.copyRecurrenceSpaces({});
  assert.deepEqual(JSON.parse(copied), report);
});

test('actual WASM finds and replays a complete space', {skip: !process.env.POLYTOOL_WASM_DIR}, async () => {
  const dir = process.env.POLYTOOL_WASM_DIR;
  const sandbox = vm.createContext({TextDecoder, TextEncoder, WebAssembly});
  vm.runInContext(fs.readFileSync(path.join(dir, 'polytool_web.js'), 'utf8'), sandbox);
  const wasm = vm.runInContext('wasm_bindgen', sandbox);
  await wasm(fs.readFileSync(path.join(dir, 'polytool_web_bg.wasm')));
  const request = {source: {coefficients: [['0', '1']], seeds: [['1']]},
    shapes: [{kind: 'box', length: 1, variable_degree: 1}]};
  const report = JSON.parse(wasm.find_recurrence_spaces(JSON.stringify(request)));
  assert.equal(report.schema, 'polytool.recurrence-spaces.v1');
  assert.equal(report.results[0].evidence, 'conditional_on_source');
  assert.equal(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(report))).verified_completed_claims, true);
  const example = html.match(/<textarea id="spaceRequest"[^>]*>([^<]+)<\/textarea>/)[1];
  const menu = JSON.parse(wasm.find_recurrence_spaces(example));
  assert.equal(menu.results.length, 3);
  assert.equal(menu.frontier.length, 3);
  assert.ok(menu.results.every(r => r.evidence === 'conditional_on_source'));
  assert.equal(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(menu))).verified_completed_claims, true);
  report.results[0].certificate.basis[0][0] = '100';
  assert.ok(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(report))).error);
});
