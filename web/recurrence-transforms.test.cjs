// Run: node --test web/recurrence-transforms.test.cjs
// POLYTOOL_WASM_DIR=/path/to/no-modules-bundle also runs the actual WASM checks;
// POLYTOOL_CLI=/path/to/polytool additionally replays browser reports with the CLI.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const {spawnSync} = require('node:child_process');
const {test} = require('node:test');

const html = fs.readFileSync(process.env.POLYTOOL_HTML || path.join(__dirname, 'index.html'), 'utf8');

function slice(start, end) {
  const from = html.indexOf(start);
  const to = html.indexOf(end, from);
  assert.ok(from >= 0 && to > from, `missing ${start} .. ${end}`);
  return html.slice(from, to);
}

// A minimal DOM: elements are created on first use and keep value/checked/hidden.
function fakeDocument() {
  const elements = new Map();
  const element = id => {
    if (!elements.has(id)) {
      elements.set(id, {id, value: '', checked: false, hidden: false, disabled: false, textContent: '',
        style: {}, dataset: {}, appendChild() {}});
    }
    return elements.get(id);
  };
  const params = ['schedule', 'schedule', 'schedule', 'ab', 'ab', 'ce', 'ce', 'power', 'weight', 'weight', 'antecedent']
    .map(param => ({dataset: {transformParam: param}, hidden: false}));
  return {
    elements, params,
    getElementById: element,
    querySelectorAll: selector => (selector === '[data-transform-param]' ? params : []),
    querySelector: selector => ({textContent: (/label\[for="([^"]+)"\]/.exec(selector) || [, 'field'])[1]}),
    createElement: () => ({}),
  };
}

function pageContext(extra) {
  const state = {displayed: '', status: '', copied: ''};
  const document = fakeDocument();
  const context = vm.createContext(Object.assign({
    window: {}, document, state, BigInt, URL, Blob: class {}, setTimeout: () => {},
    performance: {now: () => 0},
    navigator: {clipboard: {writeText: async text => { state.copied = text; }}},
    setResult: value => { state.displayed = value; },
    setStatus: value => { state.status = value; },
    katexInline: tex => '<k>' + String(tex).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;') + '</k>',
    getInput: () => document.getElementById('input').value,
    wasmReady: true, ASSET_VERSION: 'test', location: {href: 'http://localhost/'},
  }, extra || {}));
  vm.runInContext(slice('function span(cls, text)', 'function katexInline(tex)'), context);
  vm.runInContext(slice('// --- Transform tab', '// --- OEIS lookup ---'), context);
  return {context, state, document};
}

const unsafe = '<img src=x onerror=alert(1)>';

test('evidence labels are visually distinct and unknown labels stay visible and escaped', () => {
  const {context} = pageContext();
  const tones = new Set();
  for (const code of ['certified_all_indices', 'transported_identity', 'conditional_theorem_match',
    'checked_supplied_rows', 'finite_fit', 'proposed_interval', 'failed_on_rows']) {
    const badge = context.evidenceBadge(code);
    const tone = /evidence-(\w+)/.exec(badge)[1];
    tones.add(tone);
    assert.match(badge, new RegExp('data-evidence="' + code + '"'));
  }
  assert.equal(tones.size, 7, 'each evidence class has its own style');
  assert.match(context.evidenceBadge('finite_fit'), />Finite fit</);
  assert.match(context.evidenceBadge('certified_by_reconstruction'), /evidence-certified/);
  assert.match(context.evidenceBadge('exact_fit'), /evidence-fit/);
  const unknown = context.evidenceBadge(unsafe);
  assert.match(unknown, /evidence-neutral/);
  assert.doesNotMatch(unknown, /<img/);
  assert.match(context.evidenceBadge('future_label_kind'), />Future label kind</);
  assert.match(context.evidenceBadge('__proto__'), /evidence-neutral/);
  assert.match(context.evidenceLegendHtml(), /What the labels mean/);
});

test('rule texts become TeX without changing coefficients', () => {
  const {context} = pageContext();
  const tex = (text, options) => context.ruleTextTex(text, options);
  assert.equal(tex('D Q[n+1] - D Q[n] - Q[n] = 0'), "Q'_{n+1} - Q'_{n} - Q_{n} = 0");
  assert.equal(tex('x D^2 Q[n+2] + (-3 x - 2 n x^2) D Q[n+1]'), "x Q''_{n+2} + (-3 x - 2 n x^{2}) Q'_{n+1}");
  assert.equal(tex('D^5 Q[n]'), 'Q^{(5)}_{n}');
  assert.equal(tex('[x^0] Q[n+1] - [x^0] Q[n] = 0'), '[x^{0}] Q_{n+1} - [x^{0}] Q_{n} = 0');
  assert.equal(tex('(1) P[n+1] = (2X) P[n+0] + (-1) D^1(P[n+0])'), "P_{n+1} = (2x) P_{n} + (-1) P'_{n}");
  assert.equal(tex("P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)", {name: 'Q', variable: 'x'}),
    "Q_{n} = (1 + 2x) Q_{n-1} + (x + x^{2}) Q'_{n-1}");
  assert.equal(tex('\\left(1 + nt^{2}\\right) P(n) = 3/2 P^{(3)}(n-2)', {name: 'Q', variable: 'x'}),
    '\\left(1 + nx^{2}\\right) Q_{n} = \\frac{3}{2} Q^{(3)}_{n-2}');
  assert.equal(tex("(1 - t) a(n+1) = (t - t^2) a'(n) + (t) b(n)", {variable: 'x'}),
    "(1 - x) a_{n+1} = (x - x^{2}) a'_{n} + (x) b_{n}");
  assert.equal(context.rowPolyTex(['1', '-3/2', '0', '9007199254740993']), '1 - \\frac{3}{2}x + 9007199254740993x^{3}');
  assert.equal(context.rowPolyTex(['0', '-1']), '-x');
  assert.equal(context.rowPolyTex([]), '0');
  assert.equal(context.rowPolyTex(['\\href{x}']), '\\text{?}');
  assert.equal(context.intervalText({lower: {value: '-1', closed: true}, upper: {value: '0'}}), '[-1, 0]');
  assert.equal(context.intervalText({upper: {value: '0', closed: false}}), '(−∞, 0)');
  assert.equal(context.intervalText({lower: {value: '0', per_index: '-1'}}), '[-1·n, ∞)');
  assert.equal(context.indexRangeText(1, null), 'for every n ≥ 1');
  assert.equal(context.indexRangeText(1, 8), 'for n = 1…8');
});

function syntheticTransformReport() {
  return {
    schema: 'polytool.recurrence-transform.v1', request: {transformation: {kind: 'f_polynomial'}},
    transformation: {kind: 'f_polynomial', row_map: 'Q_n(x) = ' + unsafe, inverse: 'P_n(y) = ...',
      index_domain: 'n >= 1', degree_schedule: 'd_n = n - 1', components: [{name: 'Q'}]},
    source_rule: {schema: 'polytool.recurrence.v1'}, row_origin: 1,
    source_rows: [['1'], ['1', '1']], transformed_states: [[['1']], [['1', '2']]],
    reconstruction_verified: true,
    conditions: [{condition: 'deg P_n <= d_n', status: 'certified_all_indices', detail: unsafe}],
    transported: {evidence: 'transported_identity', valid_from_index: 2, text: "P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)",
      holds_on_transformed_rows: true, checked_targets: 1, scalar: {schema: 'polytool.recurrence.v1'}},
    independent_search: {evidence: 'finite_fit', status: 'found', text: 'P(n) = (1 + 2t) P(n-1)', holds_on_transformed_rows: true},
    root_domain: {claims: [
      {evidence: 'finite_exact_observation', subject: 'source', component: 0, property: 'real_rooted_in',
        interval: {lower: {value: '-7'}, upper: {value: '0'}}, first_index: 1, last_index: 8, derivation: 'exact', conditional_on: []},
      {evidence: 'proposed_interval', subject: 'transformed', component: 0, property: 'real_roots_in',
        interval: {lower: {value: '-1'}, upper: {value: '0'}}, first_index: 1, derivation: 'window', conditional_on: [unsafe]},
      {evidence: 'certified_family_invariant', subject: 'source', component: 0, property: 'real_rooted_in',
        interval: {upper: {value: '0'}}, first_index: 1, derivation: 'Möbius', conditional_on: []},
    ], implications: [{theorem: 'real_mobius_equivalence', direction: 'equivalent', statement: 's',
      hypotheses: [{condition: 'c', status: 'checked_supplied_rows', detail: 'd'}]}],
    recurrence_theorem: {matched: true, theorem: 'rolle_derivative_of_product', form: 'Q_m = ...',
      attempts: [{interval: {lower: {value: '-1'}, upper: {value: '0'}}, origin: 'default_window', success: true, detail: 'ok'}]},
    observations: [{subject: 'transformed', component: 0, index: 2, real_rooted: true,
      roots: [{lower: '-1/2', upper: '-1/2', multiplicity: 1}, {lower: '-1', upper: '-3/4', multiplicity: 2}]}]},
  };
}

test('transform reports separate transported, fitted, checked and certified results', () => {
  const {context} = pageContext();
  const card = context.renderReportCard({report: syntheticTransformReport(), origin: 'computed', replay: null}, 0);
  assert.doesNotMatch(card, /<img/);
  assert.match(card, /&lt;img src=x onerror=alert\(1\)&gt;/);
  const transported = card.slice(card.indexOf('Transported recurrence'), card.indexOf('Independent fit'));
  assert.match(transported, /evidence-transported[^>]*>Transported identity/);
  assert.doesNotMatch(transported, /evidence-fit/);
  assert.match(transported, /Q_\{n\} = \(1 \+ 2x\) Q_\{n-1\}/);
  const fit = card.slice(card.indexOf('Independent fit'), card.indexOf('Root location'));
  assert.match(fit, /evidence-fit[^>]*>Finite fit/);
  assert.doesNotMatch(fit, /evidence-certified|evidence-transported/);
  const roots = card.slice(card.indexOf('Root location'));
  // Certified claims are listed first, finite observations later.
  assert.ok(roots.indexOf('Certified for the family') < roots.indexOf('Exact on listed rows'));
  assert.match(roots, /Proposed window/);
  assert.match(roots, /Depends on: &lt;img/);
  assert.match(roots, /\[-1, 0\]/);
  assert.match(roots, /\(-1, -3\/4\) ×2/);
  assert.match(card, /Certified for all indices/);
  assert.match(card, /The inverse reproduces all 2 source rows exactly/);
  assert.match(card, /checkReportTheorems\(0,&#39;transported&#39;\)/);
  assert.match(card, /Copy report JSON/);
  assert.match(card, /Replay \(verify\)/);
  assert.doesNotMatch(card, /undefined|NaN/);
});

test('replay states, imports and unknown future report types', async () => {
  const {context, state} = pageContext();
  const report = syntheticTransformReport();
  assert.match(context.renderReportCard({report, origin: 'imported', replay: null}, 0), /not yet replayed/);
  assert.match(context.renderReportCard({report, origin: 'imported', replay: {ok: true}}, 0), /Replay passed/);
  const failed = context.renderReportCard({report, origin: 'imported', replay: {error: unsafe}}, 0);
  assert.match(failed, /Replay failed: &lt;img/);
  assert.match(failed, /<details><summary>Show the unverified content/);
  // A later report type renders through the generic evidence renderer.
  const future = {schema: 'polytool.recurrence-representation-search.v1', request: {evidence: 'ignored'},
    candidates: [{evidence: 'certified_by_reconstruction', text: 'Q[n+1] = x Q[n]'},
      {nested: {evidence: 'some_new_label', statement: unsafe}}]};
  const generic = context.renderReportCard({report: future, origin: 'imported', replay: null}, 0);
  assert.match(generic, /No dedicated view for polytool.recurrence-representation-search.v1/);
  assert.match(generic, /Proved by reconstruction/);
  assert.match(generic, /Some new label/);
  assert.doesNotMatch(generic, /ignored|<img/);
  // Import validation and clipboard export.
  context.document.getElementById('transformImportText').value = '{"no": "schema"}';
  context.replayImportedReport();
  assert.match(state.displayed, /Cannot import: a polytool report is a JSON object with a schema field/);
  vm.runInContext('transformEntries = [{report: ' + JSON.stringify(report) + ', origin: "computed", replay: null}]', context);
  await context.copyTransformReport({}, 0);
  assert.deepEqual(JSON.parse(state.copied), report);
  context.navigator.clipboard.writeText = async () => { throw new Error('denied'); };
  await context.copyTransformReport({}, 0);
  assert.match(state.status, /Could not copy the report/);
});

test('theorem and root reports state what remains to be shown', () => {
  const {context} = pageContext();
  const theorem = {schema: 'polytool.recurrence-theorems.v1', recurrence_text: 'P(n) = (1 + t) P(n-1)',
    recurrence_json: {}, finite_evidence: {source: 'supplied_recurrence', status: 'finite_evidence_only', note: 'n'},
    theorem_matches: [{theorem: 'fixed_delay_predecessor_window', status: 'conditional_theorem_match',
      proof_status: 'paper_proof_not_lean_checked', family: 'f', conclusion: 'c',
      window: {inequality: 'B - A <= 1', statement: 's'}, hypotheses: [{statement: 'h', check: 'exact_structural'}],
      boundary: {status: 'canonical', rows_match_canonical: true, contract: 'k'},
      remaining_obligations: ['Show that the intended sequence satisfies this recurrence for every n'],
      reference: 'javascript:alert(1)'}], rejected_theorems: [{theorem: 'weighted_renewal_window', reason: 'r'}]};
  const card = context.renderReportCard({report: theorem, origin: 'computed', replay: null}, 0);
  assert.match(card, /evidence-fit[^>]*>Finite evidence only/);
  assert.match(card, /evidence-theorem[^>]*>Conditional theorem match/);
  assert.match(card, /Paper proof, not Lean-checked/);
  assert.match(card, /Still to show:.*intended sequence satisfies/s);
  assert.match(card, /Fixed-delay predecessor window/);
  assert.doesNotMatch(card, /href="javascript/);
  const roots = {schema: 'polytool.root-enclosures.v1', request: {interval: {lower: {value: '-1'}, upper: {value: '0'}}},
    items: [{index: 0, degree: 2, real_rooted: true, real_root_count: 2, real_rooted_in_interval: false,
      roots: [{lower: '-3', upper: '-5/2', multiplicity: 1}, {lower: '0', upper: '0', multiplicity: 1}]}]};
  const rootCard = context.renderReportCard({report: roots, origin: 'computed', replay: null}, 0);
  assert.match(rootCard, /Exact for these polynomials only/);
  assert.match(rootCard, /\(-3, -5\/2\), 0/);
});

test('guided requests follow the form and reject bad parameters', () => {
  const {context, document} = pageContext();
  vm.runInContext('transformCatalog = [{id: "A008292", name: "Eulerian"}, {id: "A007318", name: "Pascal"}]', context);
  context.loadTransformExample('eulerian-f');
  let job = context.buildTransformRequest();
  assert.equal(job.endpoint, 'transform');
  assert.deepEqual(JSON.parse(JSON.stringify(job.request)), {oeis_id: 'A008292', row_count: 8,
    transformation: {kind: 'f_polynomial', schedule: {first_index: 1, first_degree: 0, step: 1}},
    search: {}, root_domain: {}});
  document.getElementById('transformIntervalLower').value = '-1';
  job = context.buildTransformRequest();
  assert.deepEqual(JSON.parse(JSON.stringify(job.request.root_domain)), {candidate_intervals: [{lower: {value: '-1'}}]});
  context.loadTransformExample('binomial-factorial');
  job = context.buildTransformRequest();
  assert.equal(job.endpoint, 'implicit');
  assert.equal(job.request.root_domain, undefined);
  assert.equal(document.getElementById('transformRootDomain').disabled, true);
  assert.deepEqual(JSON.parse(JSON.stringify(job.request.transformation)),
    {kind: 'diagonal', weight: {kind: 'factorial_power', power: 1}});
  document.getElementById('transformOeis').value = 'A999999';
  assert.throws(() => context.buildTransformRequest(), /not in the bundled catalog/);
  document.getElementById('transformOeis').value = 'A007318';
  document.getElementById('transformKind').value = 'affine';
  document.getElementById('transformA').value = '1/0';
  assert.throws(() => context.buildTransformRequest(), /rational number/);
  document.getElementById('transformSource').value = 'recurrence';
  document.getElementById('transformRecurrence').value = '{';
  assert.throws(() => context.buildTransformRequest(), /not valid JSON/);
  // Editing the JSON switches to a custom request until a form setting changes.
  document.getElementById('transformRequest').value = '{"custom": true}';
  document.getElementById('transformEndpoint').value = 'implicit';
  context.editTransformRequest();
  assert.deepEqual(JSON.parse(JSON.stringify(context.prepareTransformRequest())), {endpoint: 'implicit', text: '{"custom": true}'});
  assert.match(document.getElementById('transformSummary').textContent, /Custom request/);
});

async function loadWasm() {
  const dir = process.env.POLYTOOL_WASM_DIR;
  const sandbox = vm.createContext({TextDecoder, TextEncoder, WebAssembly});
  vm.runInContext(fs.readFileSync(path.join(dir, 'polytool_web.js'), 'utf8'), sandbox);
  const wasm = vm.runInContext('wasm_bindgen', sandbox);
  await wasm(fs.readFileSync(path.join(dir, 'polytool_web_bg.wasm')));
  return wasm;
}

function cliVerify(subcommand, report) {
  const result = spawnSync(process.env.POLYTOOL_CLI, [subcommand, '--verify'], {input: JSON.stringify(report), encoding: 'utf8'});
  assert.equal(result.status, 0, result.stderr || result.stdout);
}

test('actual WASM: every example transforms, renders, replays and rejects tampering', {
  skip: !process.env.POLYTOOL_WASM_DIR,
}, async () => {
  const wasm = await loadWasm();
  const {context, document} = pageContext({wasm_bindgen: wasm});
  for (const fn of vm.runInContext('TRANSFORM_WORKER_FUNCTIONS', context)) assert.equal(typeof wasm[fn], 'function', fn);
  context.populateTransformCatalog();
  const catalog = vm.runInContext('transformCatalog', context);
  assert.ok(catalog.length > 100 && catalog.some(entry => entry.id === 'A065826'));
  const expected = {
    'eulerian-f': report => {
      assert.equal(report.transported.evidence, 'transported_identity');
      assert.equal(report.independent_search.evidence, 'finite_fit');
      assert.ok(report.root_domain.claims.some(c => c.evidence === 'certified_family_invariant'));
      assert.ok(report.root_domain.claims.some(c => c.evidence === 'completed_recurrence_induction'));
    },
    'eulerian-symmetric': report => {
      assert.ok(report.transported.vector);
      assert.deepEqual(report.transformation.components.map(c => c.name), ['a', 'b']);
    },
    'narayana-reversal': report => assert.ok(report.transported.scalar),
    'binomial-factorial': report => assert.equal(report.assessment.status, 'equivalent_description'),
    'weighted-eulerian': report => assert.equal(report.independent_search.evidence, 'certified_by_reconstruction'),
  };
  for (const [name, check] of Object.entries(expected)) {
    context.loadTransformExample(name);
    const job = context.prepareTransformRequest();
    assert.equal(JSON.parse(document.getElementById('transformRequest').value).transformation.kind,
      JSON.parse(job.text).transformation.kind);
    const raw = (job.endpoint === 'implicit' ? wasm.transform_recurrence_implicit : wasm.transform_recurrence)(job.text);
    const report = JSON.parse(raw);
    assert.ok(!report.error, name + ': ' + raw);
    check(report);
    const card = context.renderReportCard({report, origin: 'computed', replay: null}, 0);
    assert.match(card, /class="evidence evidence-/);
    assert.doesNotMatch(card, /undefined|NaN|\[object Object\]/, name);
    // Compact, pretty-printed and reparsed forms all replay exactly.
    for (const text of [raw, JSON.stringify(report, null, 2), JSON.stringify(JSON.parse(JSON.stringify(report)))]) {
      assert.equal(JSON.parse(wasm.verify_report(text)).verified, true, name);
    }
    assert.equal(JSON.parse((job.endpoint === 'implicit' ? wasm.verify_recurrence_implicit
      : wasm.verify_recurrence_transform)(raw)).verified, true);
    if (process.env.POLYTOOL_CLI) {
      const subcommand = job.endpoint === 'implicit' ? 'recurrence-implicit' : 'recurrence-transform';
      cliVerify(subcommand, report);
      // The CLI computes the same report from the same request, and the browser replays it.
      const cli = spawnSync(process.env.POLYTOOL_CLI, [subcommand], {input: job.text, encoding: 'utf8'});
      assert.equal(cli.status, 0, cli.stderr);
      assert.deepEqual(JSON.parse(cli.stdout), report, name + ': CLI and browser reports agree');
      assert.equal(JSON.parse(wasm.verify_report(cli.stdout)).verified, true);
    }
    const tampered = JSON.parse(raw);
    const rows = job.endpoint === 'implicit' ? tampered.transformed_rows : tampered.transformed_states.map(s => s[0]);
    rows[rows.length - 1][0] = String(BigInt(rows[rows.length - 1][0]) + 1n);
    assert.ok(JSON.parse(wasm.verify_report(JSON.stringify(tampered))).error, name + ' tampered rows');
    const relabelled = JSON.parse(raw.replace('"finite_fit"', '"certified_by_reconstruction"')
      .replace('"transported_identity"', '"certified_all_indices"'));
    if (JSON.stringify(relabelled) !== JSON.stringify(report)) {
      assert.ok(JSON.parse(wasm.verify_report(JSON.stringify(relabelled))).error, name + ' relabelled');
    }
    // Theorem and root sub-reports for the same result.
    if (report.transported && report.transported.scalar) {
      const theorems = JSON.parse(wasm.recurrence_theorems(JSON.stringify(report.transported.scalar)));
      assert.equal(theorems.schema, 'polytool.recurrence-theorems.v1');
      assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(theorems))).verified, true);
    }
    const states = report.transformed_rows || report.transformed_states.map(s => s[0]);
    const roots = JSON.parse(wasm.root_enclosures(JSON.stringify({polynomials: states, first_index: report.row_origin})));
    assert.equal(roots.items.length, states.length);
    assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(roots))).verified, true);
  }
});

test('actual WASM: theorem matches, enclosures, spaces and browser limits', {
  skip: !process.env.POLYTOOL_WASM_DIR,
}, async () => {
  const wasm = await loadWasm();
  const {context} = pageContext({wasm_bindgen: wasm});
  const record = {schema: 'polytool.recurrence.v1', first_index: 1, initial_polynomials: [['1', '1']],
    recurrence: {terms: [{offset: 1, deriv_order: 0, sign: 'none', coeff: {coeffs: [['1', '1']]}}]}};
  const theorems = JSON.parse(wasm.recurrence_theorems(JSON.stringify(record)));
  assert.equal(theorems.theorem_matches[0].status, 'conditional_theorem_match');
  assert.equal(theorems.finite_evidence.status, 'finite_evidence_only');
  const card = context.renderReportCard({report: theorems, origin: 'computed', replay: null}, 0);
  assert.match(card, /Conditional theorem match/);
  assert.match(card, /Finite evidence only/);
  assert.match(card, /P_\{n\} = \\bigl\(1 \+ x\\bigr\) P_\{n-1\}/);
  assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(theorems))).verified, true);
  theorems.theorem_matches[0].window.inequality = 'B - A <= 9';
  assert.match(JSON.parse(wasm.verify_report(JSON.stringify(theorems))).error, /differs/);
  assert.match(JSON.parse(wasm.recurrence_theorems('{"schema": 1}')).error, /invalid recurrence JSON/);

  const huge = '9007199254740993';
  const enclosures = JSON.parse(wasm.root_enclosures(JSON.stringify({polynomials: [[huge, '1']],
    interval: {upper: {value: '0'}}})));
  assert.deepEqual(enclosures.items[0].roots, [{lower: '-' + huge, upper: '-' + huge, multiplicity: 1}]);
  assert.equal(enclosures.items[0].real_rooted_in_interval, true);
  enclosures.items[0].real_rooted_in_interval = false;
  assert.match(JSON.parse(wasm.verify_report(JSON.stringify(enclosures))).error, /differ/);
  assert.match(JSON.parse(wasm.root_enclosures(JSON.stringify({polynomials: [['1', '1']],
    max_width: '1/' + (2n ** 300n).toString()}))).error, /2\^-256/);
  assert.match(JSON.parse(wasm.transform_recurrence(JSON.stringify({oeis_id: 'A008292',
    transformation: {kind: 'f_polynomial', schedule: {first_index: 1, first_degree: 0, step: 1}},
    limits: {max_rows: 100000}}))).error, /service ceiling/);
  assert.match(JSON.parse(wasm.verify_report('{"schema": "polytool.something-else.v1"}')).error, /no browser replay/);

  // Factorial rescaling in the full recurrence search (Q rows, uniform bound).
  const spaceContext = vm.createContext({});
  vm.runInContext(slice('function recurrenceSpaceExample(', 'function invalidateRecurrenceSpaces('), spaceContext);
  const spaces = JSON.parse(wasm.find_recurrence_spaces(JSON.stringify(spaceContext.recurrenceSpaceExample('binomial'))));
  assert.equal(spaces.results[0].evidence, 'conditional_on_source');
  assert.ok(spaces.results[0].rescaling_bound);
  assert.match(spaces.results[0].representatives[0].formula, /Q\[n\+1\]/);
  assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(spaces))).verified, true);
  const formulaContext = vm.createContext({katexInline: tex => tex});
  vm.runInContext(slice('function spaceFormulaHtml(', 'function renderRecurrenceSpaces('), formulaContext);
  assert.match(formulaContext.spaceFormulaHtml(spaces.results[0].representatives[0].formula), /Q_\{n\+1\}.*Q'_\{n\}|Q_\{n\}'/);
});
