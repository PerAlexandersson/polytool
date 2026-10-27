// Run: node --test web/recurrence-transforms.test.cjs
// POLYTOOL_WASM_DIR=/path/to/no-modules-bundle also runs the actual WASM checks;
// POLYTOOL_CLI=/path/to/polytool additionally replays browser reports with the CLI.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const {spawnSync} = require('node:child_process');
const {test} = require('node:test');

const pageJs = require('./tests/page_source.cjs').loadPage().js;

function slice(start, end) {
  const from = pageJs.indexOf(start);
  const to = pageJs.indexOf(end, from);
  assert.ok(from >= 0 && to > from, `missing ${start} .. ${end}`);
  return pageJs.slice(from, to);
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
  assert.equal(context.intervalText({lower: {value: '0', per_index: '-1'}}), '[-n, ∞)');
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
  const future = {schema: 'polytool.future-report.v1', request: {evidence: 'ignored'},
    candidates: [{evidence: 'certified_by_reconstruction', text: 'Q[n+1] = x Q[n]'},
      {nested: {evidence: 'some_new_label', statement: unsafe}}]};
  const generic = context.renderReportCard({report: future, origin: 'imported', replay: null}, 0);
  assert.match(generic, /No dedicated view for polytool.future-report.v1/);
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

test('index-dependent interval ends read as formulas in n', () => {
  const {context} = pageContext();
  const end = (value, per_index, closed) => ({value, per_index, closed});
  assert.equal(context.intervalText({lower: end('-1', '-1', true), upper: {value: '0'}}), '[-1-n, 0]');
  assert.equal(context.intervalText({lower: {value: '-1'}, upper: end('0', '1', false)}), '[-1, n)');
  assert.equal(context.intervalText({lower: end('0', '-1/2', true), upper: end('3/2', '3/2', true)}), '[-n/2, 3/2+3n/2]');
  assert.equal(context.intervalText({lower: end('-5', '2', true)}), '[-5+2n, ∞)');
  assert.equal(context.intervalText({lower: end('-5', '0', true), upper: {value: '1'}}), '[-5, 1]');
  assert.equal(context.intervalText({lower: end('<img>', '1', true)}), '[?, ∞)');
  assert.equal(context.intervalText({lower: end('0', 'x', true)}), '[?, ∞)');
});

function syntheticImplicitReport() {
  const certifiedInduction = {evidence: 'completed_recurrence_induction', subject: 'source', component: 0, property: 'real_rooted_in',
    interval: {lower: {value: '-1', per_index: '-1'}, upper: {value: '0'}}, first_index: 0, derivation: 'Rolle on the defining rule', conditional_on: []};
  return {
    schema: 'polytool.recurrence-implicit.v1', request: {transformation: {kind: 'diagonal'}},
    transformation: {kind: 'diagonal', row_map: 'Q_n = T_n P_n', inverse: 'P_n = ...'},
    map_conditions: [], row_origin: 0, source_rows: [['1'], ['1', '1']], transformed_rows: [['1'], ['1', '1']],
    reconstruction_verified: true, transports: [],
    description: {first_index: 0, seeds: [['1']], relations: [], exceptional_index_conditions: [], boundary_conditions: []},
    assessment: {status: 'equivalent_description', seeds_consistent: true},
    root_domain: {
      request: {}, conventions: ['c'], assumption_checks: [], notes: ['note ' + unsafe],
      implications: [{theorem: 'polya_schur_multiplier_sequence', direction: 'source_to_target', source_property: 'P_n real-rooted',
        target_property: 'Q_n real-rooted', statement: 's', hypotheses: [
          {condition: 'lambda defined', status: 'certified_all_indices', detail: 'd'},
          {condition: 'f_j real roots', status: 'exact_fixed_polynomial', detail: 'd'},
          {condition: 'bad', status: 'violated', detail: unsafe}]}],
      claims: [
        {evidence: 'finite_exact_observation', subject: 'transformed', component: 0, property: 'real_rooted',
          first_index: 0, last_index: 7, derivation: 'exact', conditional_on: []},
        {evidence: 'certified_family_invariant', subject: 'transformed', component: 0, property: 'real_rooted',
          first_index: 0, derivation: 'induction and multiplier sequence', conditional_on: []},
        certifiedInduction],
      source_rule_theorem: {theorem: 'rolle_derivative_of_product', matched: true, statement: 's', form: 'P_m = (1 + x) P_(m-1)',
        attempts: [{interval: {lower: {value: '-1'}, upper: {value: '-1'}}, origin: 'root_hull', success: true, detail: 'ok'},
          {interval: {lower: {value: '0'}}, origin: 'supplied', success: false, detail: 'no'}]},
      fitted_rule_theorem: {theorem: 'rolle_derivative_of_product', matched: false, statement: 's', reason: 'not of Rolle form', attempts: []},
      jensen_fixed_cap: {evidence: 'fixed_cap', cap: 7, multipliers: ['1', '1'], jensen_polynomial: ['1', '7', '21/2'],
        roots: [{lower: '-2', upper: '-3/2', multiplicity: 1}], real_rooted: true, status: 'preserves_up_to_cap',
        statement: 'every real-rooted polynomial of degree at most 7 has a real-rooted image; ' + unsafe},
      observations: [{subject: 'transformed', component: 0, index: 1, real_rooted: true,
        roots: [{lower: '-1', upper: '-1', multiplicity: 1}]}],
    },
  };
}

test('implicit reports render the root-domain section with weaker fixed-cap evidence', () => {
  const {context} = pageContext();
  const card = context.renderReportCard({report: syntheticImplicitReport(), origin: 'computed', replay: null}, 0);
  assert.doesNotMatch(card, /<img|undefined|NaN|\[object Object\]/);
  const roots = card.slice(card.indexOf('<h4>Root location</h4>'));
  assert.ok(roots.length > 0, 'the implicit report has a root-location section');
  // Certified claims first, finite observation last; intervals read as formulas in n.
  assert.ok(roots.indexOf('Certified for the family') < roots.indexOf('Exact on listed rows'));
  assert.ok(roots.indexOf('Proved by induction') < roots.indexOf('Exact on listed rows'));
  assert.match(roots, /every root in \[-1-n, 0\]/);
  assert.match(roots, /Defining rule of Pₙ: matches rolle derivative of product/);
  assert.match(roots, /Certified successor rule of Qₙ: no match/);
  assert.match(roots, /\[-1, -1\] \(root hull\)/);
  assert.match(roots, /Map theorems \(1\)/);
  assert.match(roots, /evidence-fail[^>]*data-evidence="violated"/);
  assert.match(roots, /evidence-certified[^>]*data-evidence="certified_all_indices"/);
  // The Jensen test is in its own muted block with the fixed-cap label and an explicit warning.
  const cap = roots.slice(roots.indexOf('class="cap-note"'), roots.indexOf('Jensen polynomial and its roots'));
  assert.match(cap, /evidence-checked[^>]*data-evidence="fixed_cap"[^>]*>Fixed degree cap only/);
  assert.match(cap, /evidence-checked[^>]*data-evidence="preserves_up_to_cap"/);
  assert.match(cap, /Jensen test at degree 7/);
  assert.match(cap, /Finite check only: not a proof for degree 8 or more/);
  assert.doesNotMatch(cap, /evidence-certified/);
  assert.match(roots, /J_\{7\}\(x\) = 1 \+ 7x \+ \\frac\{21\}\{2\}x\^\{2\}/);
  assert.match(roots, /&lt;img/);
  const tone = code => /evidence-(\w+)/.exec(context.evidenceBadge(code))[1];
  assert.equal(tone('fixed_cap'), 'checked');
  assert.notEqual(tone('fixed_cap'), tone('certified_family_invariant'));
  assert.notEqual(tone('fixed_cap'), tone('completed_recurrence_induction'));
  assert.equal(tone('fails_at_cap'), 'fail');
  // A report without the section renders as before; an unknown Jensen status stays visible.
  const plain = syntheticImplicitReport();
  delete plain.root_domain;
  assert.doesNotMatch(context.renderReportCard({report: plain, origin: 'computed', replay: null}, 0), /Root location/);
  const odd = syntheticImplicitReport();
  odd.root_domain.jensen_fixed_cap.status = 'future_status';
  assert.match(context.renderReportCard({report: odd, origin: 'computed', replay: null}, 0), /evidence-neutral[^>]*>Future status/);
});

test('gamma reports show the transported rule and rows in z', () => {
  const {context} = pageContext();
  const report = syntheticTransformReport();
  report.transformation = {kind: 'gamma', row_map: 'P_n(x) = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)', inverse: 'inverse',
    index_domain: 'n >= 1', components: [{name: 'Q'}]};
  report.transported.text = 'P(n) = (1 - 4t + 2nt) P(n-1) + (t - 4t^2) P\'(n-1)';
  const card = context.renderReportCard({report, origin: 'computed', replay: null}, 0);
  assert.match(card, /Q_\{n\} = \(1 - 4z \+ 2nz\) Q_\{n-1\} \+ \(z - 4z\^\{2\}\) Q'_\{n-1\}/);
  assert.match(card, /gamma polynomial in z = x\/\(1 \+ x\)²/);
  assert.match(card, /<td><k>1 \+ 2z<\/k><\/td>/);
  assert.match(card, /<td><k>1 \+ x<\/k><\/td>/);
  assert.doesNotMatch(card, /undefined|NaN/);
  const refused = syntheticTransformReport();
  refused.transformation = {kind: 'gamma', components: [{name: 'Q'}]};
  refused.transported = null;
  refused.transport_error = 'the even and odd parts are independent: ' + unsafe;
  const refusedCard = context.renderReportCard({report: refused, origin: 'computed', replay: null}, 0);
  assert.match(refusedCard, /Not transported: the even and odd parts are independent: &lt;img/);
  assert.doesNotMatch(refusedCard, /<img/);
  assert.match(vm.runInContext('TRANSFORM_KIND_HELP.gamma', context), /transported when its symmetric and antisymmetric parts are proportional/);
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
  assert.equal(document.getElementById('transformRootDomain').disabled, false);
  assert.deepEqual(JSON.parse(JSON.stringify(job.request.transformation)),
    {kind: 'diagonal', weight: {kind: 'factorial_power', power: 1}});
  // Root location applies to implicit (diagonal) transforms too, with the same candidate intervals.
  context.loadTransformExample('binomial-factorial-roots');
  job = context.buildTransformRequest();
  assert.equal(job.endpoint, 'implicit');
  assert.deepEqual(JSON.parse(JSON.stringify(job.request.root_domain)), {});
  assert.equal(document.getElementById('transformIntervalFields').hidden, false);
  document.getElementById('transformIntervalLower').value = '-1-n';
  document.getElementById('transformIntervalUpper').value = '0';
  job = context.buildTransformRequest();
  assert.deepEqual(JSON.parse(JSON.stringify(job.request.root_domain)),
    {candidate_intervals: [{lower: {value: '-1', per_index: '-1'}, upper: {value: '0'}}]});
  context.loadTransformExample('binomial-factorial');
  assert.equal(context.buildTransformRequest().request.root_domain, undefined);
  assert.equal(document.getElementById('transformIntervalFields').hidden, true);
  // Interval ends: rationals and linear functions of the index.
  const endpoint = text => {
    document.getElementById('transformIntervalLower').value = text;
    document.getElementById('transformIntervalUpper').value = '';
    return JSON.parse(JSON.stringify(context.transformCandidateInterval())).lower;
  };
  assert.deepEqual(endpoint('-3/2'), {value: '-3/2'});
  assert.deepEqual(endpoint('n'), {value: '0', per_index: '1'});
  assert.deepEqual(endpoint('-n'), {value: '0', per_index: '-1'});
  assert.deepEqual(endpoint('\u22121 \u2212 n'), {value: '-1', per_index: '-1'});
  assert.deepEqual(endpoint('1/2+3n/2'), {value: '1/2', per_index: '3/2'});
  assert.deepEqual(endpoint('2n'), {value: '0', per_index: '2'});
  assert.deepEqual(endpoint('-1/2-n/2'), {value: '-1/2', per_index: '-1/2'});
  for (const bad of ['x', '1/0', '1+', 'n/0', '2^n', '--1']) {
    document.getElementById('transformIntervalLower').value = bad;
    assert.throws(() => context.transformCandidateInterval(), /rational number or a linear function of n/, bad);
  }
  document.getElementById('transformIntervalLower').value = '';
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
    'eulerian-gamma': report => {
      assert.equal(report.transported.evidence, 'transported_identity');
      assert.ok(report.transported.scalar, 'the gamma rule is transported');
      assert.equal(report.conditions[0].status, 'certified_all_indices');
      assert.equal(report.root_domain.implications[0].theorem, 'gamma_quarter_bound');
    },
    'eulerian-growing': report => {
      const attempt = report.root_domain.recurrence_theorem.attempts[0];
      assert.equal(attempt.success, true);
      assert.deepEqual(attempt.interval.lower, {value: '-1', per_index: '-1', closed: true});
      assert.ok(report.root_domain.claims.some(c => c.evidence === 'completed_recurrence_induction'));
    },
    'binomial-factorial': report => assert.equal(report.assessment.status, 'equivalent_description'),
    'binomial-factorial-roots': report => {
      assert.equal(report.assessment.status, 'equivalent_description');
      assert.equal(report.root_domain.jensen_fixed_cap.evidence, 'fixed_cap');
      assert.ok(report.root_domain.claims.some(c => c.evidence === 'certified_family_invariant' && c.subject === 'transformed'));
    },
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

test('actual WASM: implicit root-domain export, rendering, replay and tampering', {
  skip: !process.env.POLYTOOL_WASM_DIR,
}, async () => {
  const wasm = await loadWasm();
  const {context, document} = pageContext({wasm_bindgen: wasm});
  context.populateTransformCatalog();
  context.loadTransformExample('binomial-factorial-roots');
  const job = context.prepareTransformRequest();
  assert.equal(job.endpoint, 'implicit');
  assert.deepEqual(JSON.parse(job.text).root_domain, {});
  const raw = wasm.transform_recurrence_implicit(job.text);
  const report = JSON.parse(raw);
  assert.ok(!report.error, raw);
  const root = report.root_domain;
  assert.equal(root.implications[0].theorem, 'polya_schur_multiplier_sequence');
  assert.equal(root.jensen_fixed_cap.evidence, 'fixed_cap');
  assert.equal(root.jensen_fixed_cap.status, 'preserves_up_to_cap');

  // The page shows certified claims, the fixed-cap block and the Rolle matches of both rules.
  const card = context.renderReportCard({report, origin: 'computed', replay: null}, 0);
  assert.doesNotMatch(card, /undefined|NaN|\[object Object\]/);
  assert.match(card, /<h4>Root location<\/h4>/);
  assert.match(card, /evidence-certified[^>]*data-evidence="certified_family_invariant"/);
  assert.match(card, /evidence-checked[^>]*data-evidence="fixed_cap"[^>]*>Fixed degree cap only/);
  assert.match(card, /Finite check only: not a proof for degree 8 or more/);
  assert.match(card, /Defining rule of Pₙ: matches rolle derivative of product/);
  assert.match(card, /Polya schur multiplier sequence/i);
  assert.equal(card.includes('Certified for the family'), true);
  assert.ok(card.indexOf('Certified for the family') < card.indexOf('Exact on listed rows'));

  // Replay: compact, pretty and reparsed text, through both exports; the section is recognised.
  for (const text of [raw, JSON.stringify(report, null, 2)]) {
    const checked = JSON.parse(wasm.verify_report(text));
    assert.equal(checked.verified, true);
    assert.equal(checked.root_domain, true);
    assert.equal(JSON.parse(wasm.verify_recurrence_implicit(text)).root_domain, true);
  }
  if (process.env.POLYTOOL_CLI) {
    cliVerify('recurrence-implicit', report);
    const cli = spawnSync(process.env.POLYTOOL_CLI, ['recurrence-implicit'], {input: job.text, encoding: 'utf8'});
    assert.equal(cli.status, 0, cli.stderr);
    assert.deepEqual(JSON.parse(cli.stdout), report, 'CLI and browser reports agree');
  }

  // Tampering inside the root-domain section is rejected, whatever it touches.
  const edits = {
    'Jensen evidence promoted': r => { r.root_domain.jensen_fixed_cap.evidence = 'certified_all_indices'; },
    'Jensen verdict changed': r => { r.root_domain.jensen_fixed_cap.status = 'fails_at_cap'; },
    'Jensen polynomial changed': r => { r.root_domain.jensen_fixed_cap.jensen_polynomial[1] = '8'; },
    'finite claim promoted': r => {
      r.root_domain.claims.find(c => c.evidence === 'finite_exact_observation').evidence = 'certified_family_invariant';
    },
    'claim dropped': r => { r.root_domain.claims.pop(); },
    'hypothesis status changed': r => { r.root_domain.implications[0].hypotheses[0].status = 'assumed'; },
    'request changed': r => { r.root_domain.request.max_observed_rows = 3; },
    'observation changed': r => { r.root_domain.observations[2].real_rooted = false; },
    'source rows changed': r => { r.source_rows[3][1] = '99'; },
  };
  for (const [name, edit] of Object.entries(edits)) {
    const tampered = JSON.parse(raw);
    edit(tampered);
    assert.ok(JSON.parse(wasm.verify_report(JSON.stringify(tampered))).error, name);
  }
  // A report without the section is a different, valid report: it replays and says it has no root-domain section.
  const stripped = JSON.parse(raw);
  delete stripped.root_domain;
  const strippedCheck = JSON.parse(wasm.verify_report(JSON.stringify(stripped)));
  assert.equal(strippedCheck.verified, true);
  assert.equal(strippedCheck.root_domain, false);

  // A bad root-domain section is reported by the export, and the page shows it as rejected.
  const bad = JSON.parse(job.text);
  bad.root_domain = {candidate_intervals: [{lower: {value: 'abc'}}]};
  assert.ok(JSON.parse(wasm.transform_recurrence_implicit(JSON.stringify(bad))).error);
  bad.root_domain = {unknown_field: 1};
  assert.match(JSON.parse(wasm.transform_recurrence_implicit(JSON.stringify(bad))).error, /root_domain/);

  // An index-dependent candidate interval typed into the form reaches the library and is displayed as a formula.
  context.loadTransformExample('eulerian-growing');
  assert.equal(document.getElementById('transformIntervalLower').value, '-1-n');
  const growing = context.prepareTransformRequest();
  const growingReport = JSON.parse(wasm.transform_recurrence(growing.text));
  assert.ok(!growingReport.error, JSON.stringify(growingReport));
  const growingCard = context.renderReportCard({report: growingReport, origin: 'computed', replay: null}, 0);
  assert.match(growingCard, /\[-1-n, 0\]/);
  assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(growingReport))).verified, true);

  // Gamma: the Eulerian rule is transported and shown in z; the report replays.
  context.loadTransformExample('eulerian-gamma');
  const gamma = context.prepareTransformRequest();
  const gammaReport = JSON.parse(wasm.transform_recurrence(gamma.text));
  const gammaCard = context.renderReportCard({report: gammaReport, origin: 'computed', replay: null}, 0);
  assert.match(gammaCard, /evidence-transported[^>]*data-evidence="transported_identity"/);
  assert.doesNotMatch(gammaCard, /not transported/i);
  assert.match(gammaCard, /gamma polynomial in z/);
  assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(gammaReport))).verified, true);
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

// --- Representation search view ---

function syntheticSearchReport() {
  const cost = {shift_length: 1, derivative_order: 1, coefficient_degree: 1, index_degree: 0, component_count: 1,
    reconstruction_cost: 2};
  const description = (label, evidence, kind, text, extra) => Object.assign({label, kind, evidence, text, cost,
    provenance: 'steps[0]'}, extra || {});
  return {
    schema: 'polytool.recurrence-representation-search.v1', status: 'budget_exhausted',
    request: {transformations: [{kind: 'diagonal', weight: {kind: 'factorial_power', power: 1}},
      {kind: 'affine', a: '-1', b: '0'}, {kind: 'gamma', schedule: {first_index: 0, first_degree: 0, step: 1}}],
    limits: {max_depth: 1}},
    source: {kind: 'oeis', row_origin: 0, row_count: 3},
    representations: [
      {id: 0, depth: 0, chain: [], map: {identity: true, segments: [], index_domain: 'n >= 0'}, components: ['P'],
        row_origin: 0, states: [[['1']], [['1', '1']], [['1', '2', '1']]], reconstruction_verified: true, expanded: true,
        descriptions: [description('equivalent_generating_description', 'defining', 'defining_rule', 'P(n) = (1 + t) P(n-1)')]},
      {id: 1, depth: 1, parent: 0, transformation: 0, chain: [0], map: {identity: false, index_domain: 'n >= 0',
        segments: ['Q_n(x) = ' + unsafe]}, components: ['Q'], row_origin: 0,
      states: [[['1']], [['1', '1']], [['1', '2', '1/2']]], reconstruction_verified: true, expanded: false,
      notes: ['note ' + unsafe],
      descriptions: [
        description('equivalent_generating_description', 'certified_complete_space', 'certified_space',
          '(n) P(n) = (n + t) P(n-1) + t P\'(n-1)', {assessment: {status: 'equivalent_description',
            uniqueness: {status: 'certified', statement: 'unique', reasons: [unsafe]}}, also_found_by: ['fits']}),
        description('certified_annihilating_relation', 'transported_identity', 'implicit_transport', 'D Q[n+1] - D Q[n] - Q[n] = 0'),
        description('finite_fit', 'finite_fit', 'finite_fit', 'P(n) = (1 + t) P(n-1)'),
      ],
      fits: {status: 'complete', configurations_considered: 42, fits_found: 3},
      steps: [{source_description: 0, error: 'step ' + unsafe}], spaces: []},
      {id: 2, depth: 1, parent: 0, transformation: 1, chain: [1], map: {identity: false, index_domain: 'n >= 0', segments: ['Q_n(x) = P_n(-x)']},
        components: ['Q'], row_origin: 0, states: [[['1']]], reconstruction_verified: false, expanded: true,
        equal_sequence_as: {representation: 0, other_description: 0, description: 0, certificate: 'same relations ' + unsafe},
        descriptions: [description('finite_fit', 'finite_fit', 'finite_fit', 'P(n) = (1 - t) P(n-1)')]},
    ],
    edges: [
      {parent: 0, transformation: 0, outcome: 'new_representation', representation: 1},
      {parent: 0, transformation: 1, outcome: 'new_representation', representation: 2},
      {parent: 0, transformation: 2, outcome: 'rejected', detail: 'gamma: ' + unsafe},
      {parent: 2, transformation: 1, outcome: 'duplicate_map', representation: 0, detail: 'certified equal map normal forms'},
      {parent: 2, transformation: 0, outcome: 'skipped_budget', detail: 'transformations[0] was not applied'},
    ],
    ranked: [{representation: 0, description: 0, label: 'equivalent_generating_description', text: 'P(n) = (1 + t) P(n-1)', cost},
      {representation: 1, description: 0, label: 'equivalent_generating_description', text: 'x', cost}],
    depth_frontier: [1],
    budget: [{limit: 'max_representations', representation: 2, detail: 'transformations[0] was not applied ' + unsafe}],
    steps_used: 12, conventions: ['convention ' + unsafe],
  };
}

test('representation search view marks maps, duplicates, budget events and evidence', () => {
  const {context} = pageContext();
  const card = context.renderReportCard({report: syntheticSearchReport(), origin: 'computed', replay: null}, 3);
  assert.doesNotMatch(card, /No dedicated view|<img|undefined|NaN|\[object Object\]/);
  assert.match(card, /<h3>Representation search<\/h3>/);
  assert.match(card, /&lt;img src=x onerror=alert\(1\)&gt;/);
  // The tree: source, new representations, duplicate, rejected and skipped edges, equal sequence, frontier.
  const tree = card.slice(card.indexOf('class="rep-tree"'), card.indexOf('Simplest descriptions overall'));
  assert.match(tree, /rep-link" onclick="showRepresentation\(3,0\)">#0<\/button> <strong>Source/);
  assert.match(tree, /<strong>Diagonal: coefficient of xᵏ ÷ k!<\/strong>/);
  assert.match(tree, /rep-tag-equal">same sequence as #0/);
  assert.match(tree, /rep-tag-frontier">depth limit: not expanded/);
  assert.match(tree, /rep-edge-duplicate"><span class="rep-outcome">Duplicate map<\/span>[^<]*Substitute x ↦ -1x \+ 0 has the same normal form as <button[^>]*>#0/);
  assert.match(tree, /rep-edge-rejected"><span class="rep-outcome">Rejected<\/span> Gamma/);
  assert.match(tree, /rep-edge-skipped"><span class="rep-outcome">Skipped: budget<\/span>/);
  // Budget events are visible and state that the search is incomplete.
  assert.match(card, /class="rep-status rep-status-partial">Budget exhausted: incomplete/);
  assert.match(card, /class="rep-budget"[^>]*>.*Max representations.*at representation 2/s);
  assert.match(card, /Depth limit \(1\) reached.*#1/);
  // Descriptions keep the six-class scheme: label and evidence badges with their own classes.
  const one = card.slice(card.indexOf('id="rep-3-1"'), card.indexOf('id="rep-3-2"'));
  assert.match(one, /evidence-certified[^>]*data-evidence="equivalent_generating_description"[^>]*>Equivalent generating description/);
  assert.match(one, /evidence-certified[^>]*data-evidence="certified_complete_space"[^>]*>Certified complete space/);
  assert.match(one, /evidence-transported[^>]*data-evidence="certified_annihilating_relation"[^>]*>Certified annihilating relation/);
  assert.match(one, /evidence-transported[^>]*data-evidence="transported_identity"/);
  assert.match(one, /evidence-fit[^>]*data-evidence="finite_fit"[^>]*>Finite fit/);
  assert.equal((one.match(/data-evidence="finite_fit"/g) || []).length, 1, 'label equal to evidence is shown once');
  assert.match(one, /\(n\) Q_\{n\} = \(n \+ x\) Q_\{n-1\} \+ x Q'_\{n-1\}/);
  assert.match(one, /Q'_\{n\+1\} - Q'_\{n\} - Q_\{n\} = 0/);
  assert.match(one, /shift 1 · derivative order 1 · coefficient degree 1 · index degree 0 · 1 component · reconstruction cost 2/);
  assert.match(one, /Found by: steps\[0\], fits/);
  assert.match(one, /Why the description is unique.*evidence-certified[^>]*>Certified/s);
  assert.match(one, /Fit menu: complete; 42 configurations considered, 3 fit\(s\) found/);
  assert.match(one, /A transport step failed: step &lt;img/);
  assert.match(one, /<th><k>Q_n<\/k><\/th>/);
  const two = card.slice(card.indexOf('id="rep-3-2"'));
  assert.match(two, /Same sequence as representation 0\./);
  assert.match(two, /evidence-fail[^>]*>Fails on rows/);
  // Every evidence class used by the search has a distinct tone and none is neutral.
  for (const code of ['equivalent_generating_description', 'certified_complete_space', 'certified_annihilating_relation', 'defining']) {
    assert.doesNotMatch(context.evidenceBadge(code), /evidence-neutral/, code);
  }
  assert.equal(context.transformReportSummary(syntheticSearchReport()),
    '3 representation(s), 2 certified description(s); budget exhausted, the search is incomplete.');
  // A malformed or minimal report renders without failing.
  const minimal = context.renderReportCard({report: {schema: 'polytool.recurrence-representation-search.v1'},
    origin: 'imported', replay: null}, 0);
  assert.match(minimal, /Representation search/);
  assert.doesNotMatch(minimal, /undefined|NaN/);
});

test('representation search requests follow the form and stay within browser bounds', () => {
  const {context, document} = pageContext();
  const set = (id, value) => { document.getElementById(id).value = value; };
  const request = () => JSON.parse(JSON.stringify(context.buildTransformRequest()));
  vm.runInContext(`transformCatalog = [{id: "A007318", name: "Pascal"}, {id: "A065826", name: "Weighted"}, {id: "A008292", name: "Eulerian"}]`, context);
  context.loadTransformExample('search-pascal');
  assert.equal(document.getElementById('transformMode').value, 'search');
  context.updateTransformControls();
  assert.equal(document.getElementById('transformSearchFields').hidden, false);
  assert.equal(document.getElementById('transformSingleFields').hidden, true);
  assert.equal(document.getElementById('btnTransform').textContent, 'Search representations');
  assert.match(document.getElementById('transformSummary').textContent, /representations/);
  assert.deepEqual(request(), {endpoint: 'representations', request: {oeis_id: 'A007318', row_count: 10,
    transformations: [{kind: 'diagonal', weight: {kind: 'factorial_power', power: 1}}],
    limits: {max_depth: 1, max_representations: 8}}});
  // The checklist, the shared schedule and the optional switches reach the request.
  context.loadTransformExample('search-eulerian');
  const eulerian = request();
  assert.deepEqual(eulerian.request.transformations, [{kind: 'diagonal', weight: {kind: 'inverse_binomial', m: 1}},
    {kind: 'f_polynomial', schedule: {first_index: 1, first_degree: 0, step: 1}}]);
  assert.deepEqual(eulerian.request.limits, {max_depth: 2, max_representations: 12});
  for (const id of ['reversal', 'symmetric_decomposition', 'gamma', 'shift', 'negate', 'multiply_x', 'divide_x',
    'derivative_antecedent', 'diag_binomial']) {
    document.getElementById('repMap_' + id).checked = true;
  }
  document.getElementById('repFits').checked = false;
  document.getElementById('repSpaces').checked = false;
  const all = request().request;
  assert.equal(all.transformations.length, 11);
  assert.deepEqual(all.transformations.filter(t => t.kind === 'affine'), [{kind: 'affine', a: '1', b: '1'}, {kind: 'affine', a: '-1', b: '0'}]);
  assert.deepEqual(all.transformations.find(t => t.kind === 'derivative_antecedent'), {kind: 'derivative_antecedent', g: ['0', '1']});
  assert.equal(all.fits, false);
  assert.equal(all.certified_spaces, false);
  // Bounds and validation.
  set('repMaxDepth', '4');
  assert.throws(() => context.buildTransformRequest(), /from 1 to 3/);
  set('repMaxDepth', '3');
  set('repMaxRepresentations', '25');
  assert.throws(() => context.buildTransformRequest(), /from 2 to 24/);
  set('repMaxRepresentations', '24');
  set('transformRowCount', '25');
  assert.throws(() => context.buildTransformRequest(), /from 1 to 24/);
  set('transformRowCount', '10');
  vm.runInContext('REPRESENTATION_MAPS.forEach(map => { document.getElementById("repMap_" + map.id).checked = false; })', context);
  assert.throws(() => context.buildTransformRequest(), /Select at least one transformation/);
  // Rows only and recurrence sources share the single-transformation parsing.
  document.getElementById('repMap_reversal').checked = true;
  set('transformSource', 'recurrence');
  set('transformRecurrence', '{');
  assert.throws(() => context.buildTransformRequest(), /not valid JSON/);
  set('transformRecurrence', '{"schema": "polytool.recurrence.v1"}');
  assert.equal(request().request.source_recurrence.schema, 'polytool.recurrence.v1');
  // A custom request keeps its endpoint and is sent unchanged.
  document.getElementById('transformRequest').value = '{"oeis_id": "A007318"}';
  document.getElementById('transformEndpoint').value = 'representations';
  context.editTransformRequest();
  assert.deepEqual(JSON.parse(JSON.stringify(context.prepareTransformRequest())),
    {endpoint: 'representations', text: '{"oeis_id": "A007318"}'});
  // The worker exposes the search and the single-transformation modes still build as before.
  assert.ok(vm.runInContext('TRANSFORM_WORKER_FUNCTIONS', context).includes('search_recurrence_representations'));
  assert.equal(vm.runInContext('TRANSFORM_ENDPOINT_FUNCTIONS.representations', context), 'search_recurrence_representations');
  context.loadTransformExample('eulerian-f');
  assert.equal(document.getElementById('transformMode').value, 'single');
  assert.equal(context.buildTransformRequest().endpoint, 'transform');
});

test('actual WASM: representation search runs, renders, replays and rejects tampering', {
  skip: !process.env.POLYTOOL_WASM_DIR,
}, async () => {
  const wasm = await loadWasm();
  const {context, document} = pageContext({wasm_bindgen: wasm});
  for (const name of ['search_recurrence_representations', 'verify_recurrence_representation_search']) {
    assert.equal(typeof wasm[name], 'function', name);
  }
  context.populateTransformCatalog();
  const best = report => report.ranked[0];
  const expected = {
    'search-pascal': report => {
      assert.equal(report.status, 'complete');
      const found = report.representations[1].descriptions.filter(d => d.evidence === 'certified_complete_space');
      assert.ok(found.length >= 2 && found.every(d => d.label === 'equivalent_generating_description'));
    },
    'search-eulerian': report => {
      assert.ok(report.representations.some(r => r.descriptions.some(d => d.evidence === 'certified_by_reconstruction')));
      assert.ok(report.representations.length >= 3);
    },
    'search-budget': report => {
      assert.equal(report.status, 'budget_exhausted');
      assert.ok(report.budget.some(event => event.limit === 'max_representations'));
      assert.ok(report.representations.some(r => r.equal_sequence_as));
      assert.ok(report.edges.some(edge => edge.outcome === 'skipped_budget'));
    },
  };
  for (const [name, check] of Object.entries(expected)) {
    context.loadTransformExample(name);
    const job = context.prepareTransformRequest();
    assert.equal(job.endpoint, 'representations');
    const raw = wasm.search_recurrence_representations(job.text);
    const report = JSON.parse(raw);
    assert.ok(!report.error, name + ': ' + raw);
    assert.equal(report.schema, 'polytool.recurrence-representation-search.v1');
    check(report);
    assert.ok(best(report));
    const card = context.renderReportCard({report, origin: 'computed', replay: null}, 0);
    assert.match(card, /class="rep-tree"/);
    assert.match(card, /data-evidence="equivalent_generating_description"/);
    assert.doesNotMatch(card, /No dedicated view|undefined|NaN|\[object Object\]/, name);
    if (report.status === 'budget_exhausted') assert.match(card, /Budget exhausted: this search is incomplete/);
    for (const rep of report.representations) assert.ok(card.includes('id="rep-0-' + rep.id + '"'), name + ' rep ' + rep.id);
    assert.match(context.transformReportSummary(report), /representation\(s\)/);
    // Compact, pretty-printed and reparsed forms replay, through both exports.
    for (const text of [raw, JSON.stringify(report, null, 2)]) {
      const checked = JSON.parse(wasm.verify_report(text));
      assert.equal(checked.verified, true, name);
      assert.equal(checked.representations, report.representations.length);
      assert.equal(JSON.parse(wasm.verify_recurrence_representation_search(text)).status, report.status);
    }
    if (process.env.POLYTOOL_CLI) {
      cliVerify('recurrence-search-representations', report);
      const cli = spawnSync(process.env.POLYTOOL_CLI, ['recurrence-search-representations'], {input: job.text, encoding: 'utf8'});
      assert.equal(cli.status, 0, cli.stderr);
      assert.deepEqual(JSON.parse(cli.stdout), report, name + ': CLI and browser reports agree');
    }
    // Tampered numbers, labels and budget events are rejected.
    const tamperedState = JSON.parse(raw);
    tamperedState.representations[1].states[2][0][0] = '7';
    assert.ok(JSON.parse(wasm.verify_report(JSON.stringify(tamperedState))).error, name + ' state');
    const relabelled = JSON.parse(raw);
    const target = relabelled.representations.flatMap(r => r.descriptions).find(d => d.label === 'finite_fit' || d.label === 'certified_annihilating_relation');
    if (target) {
      target.label = 'equivalent_generating_description';
      assert.ok(JSON.parse(wasm.verify_recurrence_representation_search(JSON.stringify(relabelled))).error, name + ' label');
    }
    const hidden = JSON.parse(raw);
    if (hidden.budget) {
      delete hidden.budget;
      hidden.status = 'complete';
      assert.ok(JSON.parse(wasm.verify_report(JSON.stringify(hidden))).error, name + ' hidden budget');
    }
    const reranked = JSON.parse(raw);
    reranked.ranked.reverse();
    assert.ok(JSON.parse(wasm.verify_report(JSON.stringify(reranked))).error, name + ' ranking');
  }

  // Duplicate maps, rejected maps and the depth frontier from a hand-written request.
  const mixed = {oeis_id: 'A007318', row_count: 8, limits: {max_depth: 2, max_representations: 12},
    transformations: [{kind: 'affine', a: '-1', b: '0'},
      {kind: 'gamma', schedule: {first_index: 0, first_degree: 0, step: 1}}, {kind: 'divide_power', power: 1},
      {kind: 'multiply_power', power: 1}]};
  const mixedReport = JSON.parse(wasm.search_recurrence_representations(JSON.stringify(mixed)));
  assert.ok(mixedReport.edges.some(edge => edge.outcome === 'duplicate_map'));
  assert.ok(mixedReport.edges.some(edge => edge.outcome === 'rejected'));
  // The symmetric Pascal rows have a transported gamma description (Gamma is no longer rejected at the source).
  assert.ok(mixedReport.edges.some(edge => edge.parent === 0 && edge.transformation === 1 && edge.outcome === 'new_representation'));
  const mixedCard = context.renderReportCard({report: mixedReport, origin: 'computed', replay: null}, 0);
  assert.match(mixedCard, /rep-edge-duplicate/);
  assert.match(mixedCard, /rep-edge-rejected/);
  assert.match(mixedCard, /Depth limit \(2\) reached/);
  assert.equal(JSON.parse(wasm.verify_report(JSON.stringify(mixedReport))).verified, true);
  // Rows only: every description is a finite fit, and the page says so.
  const rowsOnly = JSON.parse(wasm.search_recurrence_representations(JSON.stringify({
    rows: [['1'], ['1', '1'], ['1', '2', '1'], ['1', '3', '3', '1'], ['1', '4', '6', '4', '1'], ['1', '5', '10', '10', '5', '1']],
    row_origin: 0, limits: {max_depth: 1}, transformations: [{kind: 'diagonal', weight: {kind: 'factorial_power', power: 1}}]})));
  assert.equal(rowsOnly.source.kind, 'rows_only');
  assert.ok(rowsOnly.ranked.every(item => item.label === 'finite_fit'));
  assert.match(context.renderReportCard({report: rowsOnly, origin: 'computed', replay: null}, 0), /every description is a finite fit/);

  // Limits above the service ceiling, request size and malformed requests are rejected, never run.
  const base = {oeis_id: 'A007318', transformations: [{kind: 'reversal', schedule: {first_index: 0, first_degree: 0, step: 1}}]};
  for (const limits of [{max_depth: 5}, {max_representations: 65}, {max_steps: 100000}, {max_operations: 100000000}]) {
    assert.match(JSON.parse(wasm.search_recurrence_representations(JSON.stringify(Object.assign({limits}, base)))).error,
      /service ceiling/, JSON.stringify(limits));
  }
  assert.match(JSON.parse(wasm.search_recurrence_representations(' '.repeat(1048577))).error, /1 MiB/);
  assert.ok(JSON.parse(wasm.search_recurrence_representations('{"oeis_id": 1}')).error);
  assert.ok(JSON.parse(wasm.search_recurrence_representations('{"unknown": true}')).error);
  const oversized = JSON.parse(wasm.search_recurrence_representations(JSON.stringify(Object.assign({}, base, {limits: {max_depth: 1}}))));
  oversized.request.limits.max_depth = 99;
  assert.match(JSON.parse(wasm.verify_report(JSON.stringify(oversized))).error, /service ceiling/);
  const wrongSchema = JSON.parse(JSON.stringify(oversized));
  wrongSchema.request.limits.max_depth = 1;
  wrongSchema.schema = 'polytool.recurrence-transform.v1';
  assert.ok(JSON.parse(wasm.verify_recurrence_representation_search(JSON.stringify(wrongSchema))).error);
});
