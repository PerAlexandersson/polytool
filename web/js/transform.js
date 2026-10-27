// --- Requests ---

const TRANSFORM_PARAMETERS = {
  multiply_power: ['power'], divide_power: ['power'], reversal: ['schedule'], affine: ['ab'],
  mobius: ['ab', 'ce', 'schedule'], f_polynomial: ['schedule'], gamma: ['schedule'],
  symmetric_decomposition: ['schedule'], diagonal: ['weight'], derivative_antecedent: ['antecedent'],
};
const TRANSFORM_KIND_HELP = {
  f_polynomial: 'Qₙ = (1 + x)^dₙ Pₙ(x/(1 + x)) with the declared degree bound dₙ ≥ deg Pₙ.',
  reversal: 'Qₙ = x^dₙ Pₙ(1/x) with the declared degree bound dₙ ≥ deg Pₙ.',
  affine: 'Qₙ = Pₙ(ax + b) with rational a ≠ 0 and b, such as 3/2.',
  mobius: 'Qₙ = (cx + e)^dₙ Pₙ((ax + b)/(cx + e)) with ae − bc ≠ 0.',
  multiply_power: 'Qₙ = xʳ Pₙ.',
  divide_power: 'Qₙ = Pₙ/xʳ; every row must be divisible by xʳ.',
  symmetric_decomposition: 'Pₙ = aₙ + x bₙ with aₙ, bₙ symmetric; the recurrence becomes a two-component system.',
  gamma: 'Pₙ = Σᵢ Qₙ,ᵢ xⁱ(1 + x)^(dₙ−2i) for symmetric rows, so Qₙ is a polynomial in z = x/(1 + x)². The recurrence is transported when its symmetric and antisymmetric parts are proportional (Eulerian, Narayana, type B Eulerian); otherwise the reason is shown and only the rows and a fit are reported.',
  diagonal: 'Qₙ has the coefficient of xᵏ multiplied by λ(k); the recurrence becomes an implicit relation.',
  derivative_antecedent: 'Qₙ with Pₙ = Dᵐ(g Qₙ)/m!, m = deg g; the recurrence becomes an implicit relation.',
};
const IMPLICIT_KINDS = ['diagonal', 'derivative_antecedent'];

function transformValue(id) {
  return document.getElementById(id).value;
}

function readTransformInteger(id, min, max) {
  const field = document.getElementById(id);
  const value = Number(field.value);
  if (!String(field.value).trim() || !Number.isInteger(value) || value < min || value > max) {
    throw new Error('Choose a whole number from ' + min + ' to ' + max + ' for ' +
      document.querySelector('label[for="' + id + '"]').textContent.trim().toLowerCase() + '.');
  }
  return value;
}

function readTransformRational(id, optional) {
  const text = String(document.getElementById(id).value).trim().replace(/−/g, '-');
  if (!text && optional) return null;
  if (!/^-?\d+(\/[1-9]\d*)?$/.test(text)) {
    throw new Error('Enter a rational number such as 3 or -3/2 for ' +
      document.querySelector('label[for="' + id + '"]').textContent.trim().toLowerCase() + '.');
  }
  return text;
}

// Representation search: browser bounds are inside the library's service ceiling.
const REPRESENTATION_BOUNDS = {depth: [1, 3, 2], representations: [2, 24, 12], rows: [1, 24]};
const REPRESENTATION_MAPS = [
  {id: 'diag_fact', label: 'Diagonal: divide the coefficient of xᵏ by (k!)ᵐ', checked: true,
    build: o => ({kind: 'diagonal', weight: {kind: 'factorial_power', power: o.order}})},
  {id: 'diag_inverse_binomial', label: 'Diagonal: divide the coefficient of xᵏ by binom(k+m, m)',
    build: o => ({kind: 'diagonal', weight: {kind: 'inverse_binomial', m: o.order}})},
  {id: 'diag_binomial', label: 'Diagonal: multiply the coefficient of xᵏ by binom(k+m, m)',
    build: o => ({kind: 'diagonal', weight: {kind: 'binomial', m: o.order}})},
  {id: 'f_polynomial', label: 'f-polynomial · (1+x)ᵈ P(x/(1+x))', checked: true,
    build: o => ({kind: 'f_polynomial', schedule: o.schedule})},
  {id: 'reversal', label: 'Reversal · xᵈ P(1/x)', checked: true,
    build: o => ({kind: 'reversal', schedule: o.schedule})},
  {id: 'symmetric_decomposition', label: 'Symmetric decomposition · P = a + x b',
    build: o => ({kind: 'symmetric_decomposition', schedule: o.schedule})},
  {id: 'gamma', label: 'Gamma coordinates · P = Σ Qᵢ xⁱ(1+x)^(d−2i)', build: o => ({kind: 'gamma', schedule: o.schedule})},
  {id: 'shift', label: 'Substitute x ↦ x + 1', build: () => ({kind: 'affine', a: '1', b: '1'})},
  {id: 'negate', label: 'Substitute x ↦ −x', build: () => ({kind: 'affine', a: '-1', b: '0'})},
  {id: 'multiply_x', label: 'Multiply by x', build: () => ({kind: 'multiply_power', power: 1})},
  {id: 'divide_x', label: 'Divide by x', build: () => ({kind: 'divide_power', power: 1})},
  {id: 'derivative_antecedent', label: 'Derivative antecedent · P = D(x Q)',
    build: () => ({kind: 'derivative_antecedent', g: ['0', '1']})},
];

function transformMode() {
  return transformValue('transformMode') === 'search' ? 'search' : 'single';
}

function populateRepresentationMaps() {
  const list = document.getElementById('repMapList');
  if (!list || list.dataset.ready) return;
  list.dataset.ready = '1';
  list.innerHTML = REPRESENTATION_MAPS.map(map => '<label><input type="checkbox" id="repMap_' + map.id + '"' +
    (map.checked ? ' checked' : '') + ' onchange="transformControlsChanged()"> <span>' + esc(map.label) + '</span></label>').join('');
}

function buildRepresentationFields(request) {
  populateRepresentationMaps();
  const options = {
    schedule: {first_index: readTransformInteger('repScheduleIndex', 0, 1000000),
      first_degree: readTransformInteger('repScheduleDegree', 0, 512),
      step: readTransformInteger('repScheduleStep', 0, 16)},
    order: readTransformInteger('repWeightOrder', 1, 8),
  };
  const maps = REPRESENTATION_MAPS.filter(map => document.getElementById('repMap_' + map.id).checked);
  if (!maps.length) throw new Error('Select at least one transformation to search with.');
  request.transformations = maps.map(map => map.build(options));
  request.limits = {max_depth: readTransformInteger('repMaxDepth', REPRESENTATION_BOUNDS.depth[0], REPRESENTATION_BOUNDS.depth[1]),
    max_representations: readTransformInteger('repMaxRepresentations', REPRESENTATION_BOUNDS.representations[0],
      REPRESENTATION_BOUNDS.representations[1])};
  if (!document.getElementById('repFits').checked) request.fits = false;
  if (!document.getElementById('repSpaces').checked) request.certified_spaces = false;
}

function transformEndpointFor(kind, source) {
  return IMPLICIT_KINDS.includes(kind) && source !== 'rows' ? 'implicit' : 'transform';
}

function buildTransformation(kind) {
  const schedule = () => ({first_index: readTransformInteger('transformScheduleIndex', 0, 1000000),
    first_degree: readTransformInteger('transformScheduleDegree', 0, 512),
    step: readTransformInteger('transformScheduleStep', 0, 16)});
  switch (kind) {
    case 'multiply_power': case 'divide_power':
      return {kind, power: readTransformInteger('transformPower', 0, 64)};
    case 'affine':
      return {kind, a: readTransformRational('transformA'), b: readTransformRational('transformB')};
    case 'mobius':
      return {kind, a: readTransformRational('transformA'), b: readTransformRational('transformB'),
        c: readTransformRational('transformC'), e: readTransformRational('transformE'), schedule: schedule()};
    case 'reversal': case 'f_polynomial': case 'gamma': case 'symmetric_decomposition':
      return {kind, schedule: schedule()};
    case 'diagonal': {
      const weight = transformValue('transformWeight');
      const order = readTransformInteger('transformWeightOrder', 1, 8);
      return {kind, weight: weight === 'factorial_power' ? {kind: weight, power: order} : {kind: weight, m: order}};
    }
    case 'derivative_antecedent': {
      const g = String(transformValue('transformG')).split(/[\s,;]+/).filter(Boolean);
      if (!g.length || !g.every(c => /^-?\d+(\/[1-9]\d*)?$/.test(c))) {
        throw new Error('Enter the coefficients of g as rational numbers, lowest degree first, such as 0, 1.');
      }
      return {kind, g};
    }
  }
  throw new Error('Unknown transformation.');
}

// An interval end: a rational, or an affine function of the index n such as -1-n, n/2 or 1/2+3n/2.
function readTransformEndpoint(id) {
  const text = String(document.getElementById(id).value).trim().replace(/−/g, '-').replace(/\s+/g, '').replace(/\*/g, '');
  if (!text) return null;
  const constant = /^(-?\d+(?:\/[1-9]\d*)?)$/.exec(text);
  if (constant) return {value: constant[1]};
  // [constant] [+|-] [integer] n [/ denominator]: -1-n, n/2, 1/2+3n/2.
  const affine = /^(?:(-?\d+(?:\/[1-9]\d*)?)(?=[+-]))?([+-]?)(\d+)?n(?:\/([1-9]\d*))?$/.exec(text);
  if (affine) {
    return {value: affine[1] || '0',
      per_index: (affine[2] === '-' ? '-' : '') + (affine[3] || '1') + (affine[4] ? '/' + affine[4] : '')};
  }
  throw new Error('Enter an interval end as a rational number or a linear function of n, such as -1, 3/2, -1-n or n/2, for ' +
    document.querySelector('label[for="' + id + '"]').textContent.trim().toLowerCase() + '.');
}

function transformCandidateInterval() {
  const lower = readTransformEndpoint('transformIntervalLower');
  const upper = readTransformEndpoint('transformIntervalUpper');
  if (lower === null && upper === null) return null;
  const interval = {};
  if (lower !== null) interval.lower = lower;
  if (upper !== null) interval.upper = upper;
  return interval;
}

function buildTransformRequest() {
  const source = transformValue('transformSource');
  const searching = transformMode() === 'search';
  const kind = transformValue('transformKind');
  const endpoint = searching ? 'representations' : transformEndpointFor(kind, source);
  const rowMax = searching ? REPRESENTATION_BOUNDS.rows[1] : 64;
  const request = {};
  if (source === 'oeis') {
    const id = String(transformValue('transformOeis')).trim().toUpperCase();
    if (!/^A\d{6}$/.test(id)) throw new Error('Enter an OEIS A-number such as A008292.');
    if (transformCatalog && !transformCatalog.some(entry => entry.id === id)) {
      throw new Error(id + ' is not in the bundled catalog. Paste its recurrence JSON instead.');
    }
    request.oeis_id = id;
    request.row_count = readTransformInteger('transformRowCount', 1, rowMax);
  } else if (source === 'recurrence') {
    const text = String(transformValue('transformRecurrence')).trim();
    if (!text) throw new Error('Paste a recurrence JSON record, or use the recurrence found in the Recurrence tab.');
    if (text.length > 1048576) throw new Error('Recurrence JSON is too large (limit: 1 MiB).');
    let record;
    try { record = JSON.parse(text); } catch (error) { throw new Error('Recurrence JSON is not valid JSON: ' + error.message); }
    if (!record || typeof record !== 'object' || Array.isArray(record)) throw new Error('Recurrence JSON must be an object.');
    request.source_recurrence = record;
    request.row_count = readTransformInteger('transformRowCount', 1, rowMax);
  } else {
    const input = getInput();
    if (input.length > 1048576) throw new Error('Polynomial input is too large (limit: 1 MiB).');
    const parsed = JSON.parse(wasm_bindgen.parse_and_format(input));
    const errors = parsed.filter(row => row.error).map(row => row.error);
    if (errors.length) throw new Error(errors.join('; '));
    if (!parsed.length) throw new Error('Enter polynomial rows above, or choose a bundled family.');
    if (parsed.length > 256) throw new Error('At most 256 rows can be transformed.');
    request.rows = parsed.map(row => row.coefficients);
    request.row_origin = readTransformInteger('transformOrigin', 0, 1000000);
  }
  if (searching) {
    buildRepresentationFields(request);
    return {endpoint, request};
  }
  request.transformation = buildTransformation(kind);
  if (document.getElementById('transformSearch').checked) request.search = {};
  if (document.getElementById('transformRootDomain').checked) {
    const interval = transformCandidateInterval();
    request.root_domain = interval ? {candidate_intervals: [interval]} : {};
  }
  return {endpoint, request};
}

function transformExample(name) {
  const examples = {
    'eulerian-f': {source: 'oeis', oeis: 'A008292', rows: 8, kind: 'f_polynomial', schedule: [1, 0, 1], search: true, roots: true},
    'eulerian-symmetric': {source: 'oeis', oeis: 'A008292', rows: 8, kind: 'symmetric_decomposition', schedule: [1, 0, 1], search: false, roots: true},
    'narayana-reversal': {source: 'oeis', oeis: 'A001263', rows: 8, kind: 'reversal', schedule: [1, 0, 1], search: true, roots: true},
    'eulerian-gamma': {source: 'oeis', oeis: 'A008292', rows: 8, kind: 'gamma', schedule: [1, 0, 1], search: true, roots: true},
    'eulerian-growing': {source: 'oeis', oeis: 'A008292', rows: 8, kind: 'f_polynomial', schedule: [1, 0, 1], search: false, roots: true,
      interval: ['-1-n', '0']},
    'binomial-factorial': {source: 'oeis', oeis: 'A007318', rows: 8, kind: 'diagonal', weight: ['factorial_power', 1], search: true, roots: false},
    'binomial-factorial-roots': {source: 'oeis', oeis: 'A007318', rows: 8, kind: 'diagonal', weight: ['factorial_power', 1], search: false, roots: true},
    'weighted-eulerian': {source: 'oeis', oeis: 'A065826', rows: 9, kind: 'diagonal', weight: ['inverse_binomial', 1], search: true, roots: false},
    'search-pascal': {mode: 'search', source: 'oeis', oeis: 'A007318', rows: 10, maps: ['diag_fact'], schedule: [0, 0, 1],
      depth: 1, representations: 8},
    'search-eulerian': {mode: 'search', source: 'oeis', oeis: 'A065826', rows: 9, maps: ['diag_inverse_binomial', 'f_polynomial'],
      schedule: [1, 0, 1], depth: 2, representations: 12},
    'search-budget': {mode: 'search', source: 'oeis', oeis: 'A007318', rows: 8, maps: ['reversal', 'f_polynomial', 'diag_fact'],
      schedule: [0, 0, 1], depth: 2, representations: 6},
  };
  return examples[name] || null;
}

function loadTransformExample(name) {
  const example = transformExample(name);
  if (!example) return;
  populateRepresentationMaps();
  const set = (id, value) => { document.getElementById(id).value = value; };
  set('transformMode', example.mode || 'single');
  if (example.mode === 'search') {
    REPRESENTATION_MAPS.forEach(map => { document.getElementById('repMap_' + map.id).checked = example.maps.includes(map.id); });
    set('repScheduleIndex', example.schedule[0]);
    set('repScheduleDegree', example.schedule[1]);
    set('repScheduleStep', example.schedule[2]);
    set('repWeightOrder', example.order || 1);
    set('repMaxDepth', example.depth);
    set('repMaxRepresentations', example.representations);
    document.getElementById('repFits').checked = true;
    document.getElementById('repSpaces').checked = true;
  }
  set('transformSource', example.source);
  set('transformOeis', example.oeis);
  set('transformRowCount', example.rows);
  if (example.kind) set('transformKind', example.kind);
  if (example.mode !== 'search' && example.schedule) {
    set('transformScheduleIndex', example.schedule[0]);
    set('transformScheduleDegree', example.schedule[1]);
    set('transformScheduleStep', example.schedule[2]);
  }
  if (example.weight) {
    set('transformWeight', example.weight[0]);
    set('transformWeightOrder', example.weight[1]);
  }
  if (example.mode !== 'search') {
    document.getElementById('transformSearch').checked = example.search;
    document.getElementById('transformRootDomain').checked = example.roots;
  }
  set('transformIntervalLower', example.interval ? example.interval[0] : '');
  set('transformIntervalUpper', example.interval ? example.interval[1] : '');
  transformControlsChanged();
  setStatus('Example loaded. Press Transform.');
}

function populateTransformCatalog() {
  if (transformCatalog || !wasmReady) return;
  try {
    transformCatalog = JSON.parse(wasm_bindgen.oeis_catalog_index());
  } catch (_) {
    transformCatalog = null;
    return;
  }
  const list = document.getElementById('transformOeisList');
  list.textContent = '';
  transformCatalog.forEach(entry => {
    const option = document.createElement('option');
    option.value = entry.id;
    option.label = entry.id + ' — ' + String(entry.name).slice(0, 90);
    list.appendChild(option);
  });
  updateTransformControls();
}

function updateTransformControls() {
  populateRepresentationMaps();
  const source = transformValue('transformSource');
  const kind = transformValue('transformKind');
  const searching = transformMode() === 'search';
  const endpoint = transformEndpointFor(kind, source);
  document.getElementById('transformSingleFields').hidden = searching;
  document.getElementById('transformSearchFields').hidden = !searching;
  document.getElementById('btnTransform').textContent = searching ? 'Search representations' : 'Transform';
  document.getElementById('transformRowCount').max = searching ? REPRESENTATION_BOUNDS.rows[1] : 64;
  const params = TRANSFORM_PARAMETERS[kind] || [];
  document.querySelectorAll('[data-transform-param]').forEach(field => {
    field.hidden = !params.includes(field.dataset.transformParam);
  });
  document.getElementById('transformOeisField').hidden = source !== 'oeis';
  document.getElementById('transformRecurrenceField').hidden = source !== 'recurrence';
  document.getElementById('transformRowCountField').hidden = source === 'rows';
  document.getElementById('transformOriginField').hidden = source !== 'rows';
  const rootBox = document.getElementById('transformRootDomain');
  document.getElementById('transformIntervalFields').hidden = !rootBox.checked;
  document.getElementById('transformKindHelp').textContent = TRANSFORM_KIND_HELP[kind] || '';
  const id = String(transformValue('transformOeis')).trim().toUpperCase();
  const entry = transformCatalog && transformCatalog.find(item => item.id === id);
  document.getElementById('transformOeisName').textContent = entry ? entry.name :
    (transformCatalog ? 'Not in the bundled catalog (' + transformCatalog.length + ' families).' : '');
  const sourceText = {oeis: 'the defining recurrence of the bundled family', recurrence: 'the supplied recurrence',
    rows: 'the rows entered above only (no recurrence, so every result is checked on these rows only)'}[source];
  document.getElementById('transformSummary').textContent = transformCustom
    ? 'Custom request: the JSON below is sent unchanged.'
    : (searching ? 'Bounded search over equivalent representations, starting from '
      : (endpoint === 'implicit' ? 'Implicit relations from ' : 'Coordinate change using ')) + sourceText + '.';
}

function refreshTransformRequestPreview() {
  if (transformCustom || !wasmReady) return;
  try {
    const job = buildTransformRequest();
    document.getElementById('transformEndpoint').value = job.endpoint;
    document.getElementById('transformRequest').value = JSON.stringify(job.request, null, 2);
  } catch (_) {
    document.getElementById('transformRequest').value = '';
  }
}

function invalidateTransformResults() {
  if (transformRunning) cancelTransformJob(true);
  if (transformEntries.length) {
    transformEntries = [];
    setResult('Settings changed. Press Transform to update the results.', 'transform');
  }
}

function transformControlsChanged() {
  transformCustom = false;
  updateTransformControls();
  refreshTransformRequestPreview();
  invalidateTransformResults();
}

function transformMainInputChanged() {
  if (transformValue('transformSource') === 'rows') {
    refreshTransformRequestPreview();
    invalidateTransformResults();
  }
}

function editTransformRequest() {
  transformCustom = true;
  updateTransformControls();
  invalidateTransformResults();
}

function prepareTransformRequest() {
  if (transformCustom) {
    const text = transformValue('transformRequest');
    if (!text.trim()) throw new Error('The request JSON is empty.');
    return {endpoint: transformValue('transformEndpoint'), text};
  }
  const job = buildTransformRequest();
  document.getElementById('transformEndpoint').value = job.endpoint;
  document.getElementById('transformRequest').value = JSON.stringify(job.request, null, 2);
  return {endpoint: job.endpoint, text: JSON.stringify(job.request)};
}

// --- Worker and actions ---

const TRANSFORM_WORKER_FUNCTIONS = ['transform_recurrence', 'transform_recurrence_implicit', 'root_enclosures',
  'recurrence_theorems', 'verify_report', 'search_recurrence_representations'];
const TRANSFORM_ENDPOINT_FUNCTIONS = {transform: 'transform_recurrence', implicit: 'transform_recurrence_implicit',
  representations: 'search_recurrence_representations'};
let transformEntries = [];
let transformCatalog = null;
let transformCustom = false;
let transformWorker = null;
let transformRunning = false;
let transformJob = null;
let transformJobSerial = 0;

function createTransformWorker() {
  const base = new URL('.', location.href).href;
  const code = `
    importScripts('${base}pkg/polytool_web.js?v=${ASSET_VERSION}');
    const ready = wasm_bindgen({ module_or_path: '${base}pkg/polytool_web_bg.wasm?v=${ASSET_VERSION}' });
    const allowed = new Set(${JSON.stringify(TRANSFORM_WORKER_FUNCTIONS)});
    onmessage = async function(e) {
      await ready;
      const {id, fn, arg} = e.data;
      if (!allowed.has(fn)) { postMessage({id, data: JSON.stringify({error: 'unsupported call'})}); return; }
      try {
        postMessage({id, data: wasm_bindgen[fn](arg)});
      } catch (error) {
        postMessage({id, poisoned: true, data: JSON.stringify({error: 'internal error: ' + (error && error.message || error)})});
      }
    };
  `;
  const blob = new Blob([code], {type: 'application/javascript'});
  return new Worker(URL.createObjectURL(blob));
}

function setTransformBusy(busy, label) {
  transformRunning = busy;
  document.getElementById('btnTransform').disabled = busy;
  document.getElementById('btnTransformImport').disabled = busy;
  document.getElementById('btnTransformCancel').style.display = busy ? '' : 'none';
  if (busy) setStatus(label + '…');
}

function callTransformWorker(fn, arg, label) {
  if (transformRunning) return Promise.reject(new Error('Another computation is running; wait for it or cancel it.'));
  if (!transformWorker) transformWorker = createTransformWorker();
  setTransformBusy(true, label);
  const id = ++transformJobSerial;
  return new Promise((resolve, reject) => {
    transformJob = {id, reject};
    transformWorker.onmessage = e => {
      if (e.data.id !== id) return;
      transformJob = null;
      setTransformBusy(false);
      if (e.data.poisoned) { transformWorker.terminate(); transformWorker = null; }
      resolve(e.data.data);
    };
    transformWorker.onerror = e => {
      transformJob = null;
      setTransformBusy(false);
      if (transformWorker) transformWorker.terminate();
      transformWorker = null;
      reject(new Error('Worker error: ' + (e.message || 'unknown')));
    };
    transformWorker.postMessage({id, fn, arg});
  });
}

function cancelTransformJob(silent) {
  if (!transformRunning) return;
  if (transformWorker) transformWorker.terminate();
  transformWorker = null;
  const job = transformJob;
  transformJob = null;
  setTransformBusy(false);
  if (job) job.reject(Object.assign(new Error('cancelled'), {cancelled: true}));
  if (!silent) setStatus('Cancelled.');
}

function showTransformError(prefix, error) {
  if (error && error.cancelled) {
    if (!transformEntries.length) setResult('Computation cancelled.', 'transform');
    return;
  }
  setResult(span('err', prefix + (error && error.message || error)), 'transform');
  setStatus('Error.');
}

function displayTransformResults(status) {
  try {
    renderTransformResults();
    if (status) setStatus(status);
  } catch (error) {
    setResult(span('err', 'Could not display the report: ' + (error && error.message || error)), 'transform');
    setStatus('Error displaying the report.');
  }
}

async function runTransform() {
  if (!wasmReady) { setStatus('WASM not loaded yet.'); return; }
  if (transformRunning) return;
  let job;
  try {
    job = prepareTransformRequest();
  } catch (error) {
    setResult(span('err', 'Error: ' + error.message), 'transform');
    setStatus('Check the transformation settings.');
    return;
  }
  transformEntries = [];
  setResult(job.endpoint === 'representations' ? 'Searching representations…' : 'Transforming…', 'transform');
  const started = performance.now();
  try {
    const text = await callTransformWorker(TRANSFORM_ENDPOINT_FUNCTIONS[job.endpoint] || 'transform_recurrence',
      job.text, job.endpoint === 'representations' ? 'Searching representations' : 'Transforming');
    const data = JSON.parse(text);
    if (data.error) {
      setResult(span('err', 'Request rejected: ' + data.error), 'transform');
      setStatus('Request rejected.');
      return;
    }
    transformEntries = [{report: data, origin: 'computed', replay: null}];
    displayTransformResults(transformReportSummary(data) + ' Returned in ' + (performance.now() - started).toFixed(0) + ' ms.');
  } catch (error) {
    showTransformError('Error: ', error);
  }
}

async function addTransformSubReport(fn, arg, title, label) {
  if (transformRunning) { setStatus('Wait for the current computation or cancel it.'); return; }
  try {
    const data = JSON.parse(await callTransformWorker(fn, arg, label));
    if (data.error) { setStatus(label + ' rejected: ' + data.error); return; }
    transformEntries.push({report: data, origin: 'computed', replay: null, title});
    displayTransformResults(title + ': ' + transformReportSummary(data));
  } catch (error) {
    if (!(error && error.cancelled)) setStatus(label + ' failed: ' + (error && error.message || error));
  }
}

function checkReportTheorems(index, which) {
  const report = transformEntries[index] && transformEntries[index].report;
  if (!report) return;
  const record = which === 'source' ? report.source_rule
    : (which === 'transported' ? report.transported && report.transported.scalar
      : report.independent_search && report.independent_search.scalar);
  if (!record) { setStatus('No scalar recurrence to check.'); return; }
  const titles = {source: 'Interlacing theorems: source recurrence', transported: 'Interlacing theorems: transported recurrence',
    fit: 'Interlacing theorems: fitted recurrence'};
  addTransformSubReport('recurrence_theorems', JSON.stringify(record), titles[which], 'Checking theorems');
}

function checkSourceTheorems() {
  const text = String(transformValue('transformRecurrence')).trim();
  if (!text) { setStatus('Paste a recurrence JSON record first.'); return; }
  addTransformSubReport('recurrence_theorems', text, 'Interlacing theorems: supplied recurrence', 'Checking theorems');
}

function isolateReportRows(index, which) {
  const report = transformEntries[index] && transformEntries[index].report;
  if (!report) return;
  let rows = which === 'source' ? report.source_rows
    : (report.transformed_rows || (report.transformed_states || []).map(state => state[0] || []));
  rows = (rows || []).slice(0, 64);
  if (!rows.length) { setStatus('No rows to isolate.'); return; }
  const request = {polynomials: rows, first_index: report.row_origin || 0};
  try {
    const interval = transformCandidateInterval();
    // Root enclosures take fixed endpoints; an index-dependent candidate interval is not applied to them.
    if (interval && !(interval.lower && interval.lower.per_index) && !(interval.upper && interval.upper.per_index)) {
      request.interval = interval;
    }
  } catch (_) {}
  addTransformSubReport('root_enclosures', JSON.stringify(request),
    'Exact roots of the ' + (which === 'source' ? 'source' : 'transformed') + ' rows', 'Isolating roots');
}

async function copyTransformReport(button, index) {
  const entry = transformEntries[index];
  if (!entry) return;
  try {
    await navigator.clipboard.writeText(JSON.stringify(entry.report, null, 2));
    button.textContent = 'Copied!';
    setTimeout(() => { button.textContent = 'Copy report JSON'; }, 1500);
  } catch (_) {
    setStatus('Could not copy the report. Check clipboard permission, or use Download JSON.');
  }
}

function downloadTransformReport(index) {
  const entry = transformEntries[index];
  if (!entry) return;
  const name = String(entry.report.schema || 'polytool-report').replace(/^polytool\./, '').replace(/[^A-Za-z0-9.-]/g, '-');
  const link = document.createElement('a');
  link.href = URL.createObjectURL(new Blob([JSON.stringify(entry.report, null, 2) + '\n'], {type: 'application/json'}));
  link.download = name + '.json';
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(link.href), 0);
}

async function replayTransformReport(index) {
  const entry = transformEntries[index];
  if (!entry) return;
  if (transformRunning) { setStatus('Wait for the current computation or cancel it.'); return; }
  entry.replay = 'running';
  displayTransformResults();
  try {
    const data = JSON.parse(await callTransformWorker('verify_report', JSON.stringify(entry.report), 'Replaying'));
    entry.replay = data.error ? {error: data.error} : {ok: true, result: data};
    displayTransformResults(data.error ? 'Replay failed.' : 'Replay passed.');
  } catch (error) {
    entry.replay = error && error.cancelled ? null : {error: String(error && error.message || error)};
    displayTransformResults(error && error.cancelled ? 'Replay cancelled.' : 'Replay failed.');
  }
}

function importTransformFile(input) {
  const file = input.files && input.files[0];
  if (!file) return;
  if (file.size > 1048576) { setStatus('The report file exceeds 1 MiB.'); return; }
  const reader = new FileReader();
  reader.onload = () => { document.getElementById('transformImportText').value = String(reader.result); setStatus('Report loaded. Press Verify report.'); };
  reader.onerror = () => setStatus('Could not read the report file.');
  reader.readAsText(file);
}

function replayImportedReport() {
  if (transformRunning) { setStatus('Wait for the current computation or cancel it.'); return; }
  const text = String(document.getElementById('transformImportText').value).trim();
  let report;
  try {
    if (text.length > 1048576) throw new Error('the report exceeds 1 MiB');
    report = JSON.parse(text);
    if (!report || typeof report !== 'object' || Array.isArray(report) || typeof report.schema !== 'string') {
      throw new Error('a polytool report is a JSON object with a schema field');
    }
  } catch (error) {
    setResult(span('err', 'Cannot import: ' + error.message), 'transform');
    setStatus('Import rejected.');
    return;
  }
  transformEntries = [{report, origin: 'imported', replay: null}];
  replayTransformReport(0);
}

function useFoundRecurrence() {
  if (!window._lastRecJson) { setStatus('Find a recurrence in the Recurrence tab first.'); return; }
  document.getElementById('transformRecurrence').value = window._lastRecJson;
  document.getElementById('transformSource').value = 'recurrence';
  transformControlsChanged();
  setStatus('Using the recurrence found in the Recurrence tab.');
}

function openTransformWithFoundRecurrence() {
  switchTab('transform');
  useFoundRecurrence();
}
