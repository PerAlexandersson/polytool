const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const {test} = require('node:test');
const pageJs = require('./page_source.cjs').loadPage().js;

test('space display leads with formulas, preserves scope, and omits certificates', async () => {
  let displayed = '', copied = '', status = '';
  const sandbox = vm.createContext({
    window: {}, setResult: text => {displayed = text;},
    navigator: {clipboard: {writeText: async text => {copied = text;}}},
    setTimeout: () => {}, setStatus: text => {status = text;},
  });
  const start = pageJs.indexOf('function recurrenceSpaceSummary(');
  vm.runInContext(pageJs.slice(pageJs.indexOf('function esc('), pageJs.indexOf('async function run(')), sandbox);
  vm.runInContext(pageJs.slice(start, pageJs.indexOf('function createWorker(', start)), sandbox);
  const report = {request: {domain_start: 3, shapes: [{kind: 'box', length: 2, derivative_order: 1}]}, results: [
    {status: 'successor_exists', evidence: 'exact_fit', certificate: {basis: [['1']]},
      representatives: [{formula: '<unsafe>', singular_ranks: ['51'], requires_x_division: true},
        {formula: 'second formula', singular_ranks: [], requires_x_division: false}]},
    {status: 'absent_in_shape', evidence: 'rows_only_obstruction', certificate: {}, representatives: []},
    {status: 'incomplete_budget', evidence: 'none', reason: 'budget exhausted', representatives: [],
      automatic_degree_bound: {variable_degree: 1}},
  ]};
  assert.match(sandbox.spaceFormulaHtml('(1) P[n+1] = (2X) P[n+0] + (-1) D^1(P[n+0])'), /P_\{n\+1\} = \(2x\) P_\{n\} - P_\{n\}'/);
  sandbox.renderRecurrenceSpaces(report);
  assert.match(displayed, /2 of 3 searches finished; 1 unfinished/);
  assert.match(displayed, /2 previous terms · derivatives up to order 1/);
  assert.match(displayed, /Base index starts at n = 3/);
  assert.match(displayed, /fit the supplied rows exactly/);
  assert.match(displayed, /No recurrence within these bounds/);
  assert.match(displayed, /does not tell us whether a recurrence exists/);
  assert.match(displayed, /automatic x-degree search is unfinished/);
  assert.match(displayed, /Cannot determine the next term at base n = 51/);
  assert.match(displayed, /&lt;unsafe&gt;/);
  assert.doesNotMatch(displayed, /<unsafe>|kernel|projection|certificate|verify-space|Copy.*JSON/i);
  assert.match(displayed, /<details><summary>More formulas/);
  await sandbox.copySpaceFormula({}, 0, 0);
  assert.equal(copied, '<unsafe>');
  sandbox.navigator.clipboard.writeText = async () => {throw new Error('permission denied');};
  await sandbox.copySpaceFormula({}, 0, 0);
  assert.match(status, /Could not copy the formula/);
});

test('actual WASM finds and replays a complete space', {skip: !process.env.POLYTOOL_WASM_DIR}, async () => {
  const dir = process.env.POLYTOOL_WASM_DIR;
  const sandbox = vm.createContext({TextDecoder, TextEncoder, WebAssembly});
  vm.runInContext(fs.readFileSync(path.join(dir, 'polytool_web.js'), 'utf8'), sandbox);
  const wasm = vm.runInContext('wasm_bindgen', sandbox);
  await wasm(fs.readFileSync(path.join(dir, 'polytool_web_bg.wasm')));
  const parsedRows = JSON.parse(wasm.parse_and_format('1 + x\n9007199254740993, 2'));
  assert.deepEqual(parsedRows.map(row => row.coefficients), [['1', '1'], ['9007199254740993', '2']]);
  assert.ok(JSON.parse(wasm.parse_and_format('1\ninvalid row'))[1].error);
  const request = {source: {coefficients: [['0', '1']], seeds: [['1']]},
    shapes: [{kind: 'box', length: 1, variable_degree: 1}]};
  const report = JSON.parse(wasm.find_recurrence_spaces(JSON.stringify(request)));
  assert.equal(report.schema, 'polytool.recurrence-spaces.v1');
  assert.equal(report.results[0].evidence, 'conditional_on_source');
  assert.equal(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(report))).verified_completed_claims, true);
  const exampleStart = pageJs.indexOf('function recurrenceSpaceExample(');
  const exampleContext = vm.createContext({});
  vm.runInContext(pageJs.slice(exampleStart, pageJs.indexOf('function renderRecurrenceSpaces(', exampleStart)), exampleContext);
  const menu = JSON.parse(wasm.find_recurrence_spaces(JSON.stringify(exampleContext.recurrenceSpaceExample('default'))));
  assert.equal(menu.results.length, 3);
  assert.equal(menu.frontier.length, 3);
  assert.ok(menu.results.every(r => r.evidence === 'conditional_on_source'));
  assert.equal(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(menu))).verified_completed_claims, true);
  for (const kind of ['hermite', 'automatic']) {
    const extended = JSON.parse(wasm.find_recurrence_spaces(JSON.stringify(exampleContext.recurrenceSpaceExample(kind))));
    assert.equal(extended.results[0].evidence, 'conditional_on_source');
    assert.ok(extended.results[0].propagation);
    if (kind === 'automatic') assert.ok(extended.results[0].automatic_degree_bound);
    assert.equal(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(extended))).verified_completed_claims, true);
  }
  report.results[0].certificate.basis[0][0] = '100';
  assert.ok(JSON.parse(wasm.verify_recurrence_spaces(JSON.stringify(report))).error);
});
