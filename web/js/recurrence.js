// --- Recurrence via Web Worker ---

function renderRecurrenceMenu(selected) {
  const report = window._recurrenceMenu;
  if (!report) return;
  const items = report.representations;
  const coverage = report.coverage.map(cell =>
    'Length ' + cell.rec_len + ', derivative order ' + cell.diff_deg + ': ' +
    cell.status.replaceAll('_', ' ') + ' (' + cell.configurations_considered + ' configurations)'
  ).join('\n');
  const incomplete = report.status === 'budget_exhausted';
  let html = '<div class="recurrence-note">' +
    (incomplete ? 'Search limit reached; this is a partial comparison. ' : 'Bounded comparison completed. ') +
    'Exact fits to the supplied rows, not proofs for all n. One preferred fit per length/order bound pair; ' +
    '“Tradeoff” marks candidates not dominated in length, derivative order, coefficient degrees, size, and checked start.</div>' +
    '<div class="recurrence-stats">' + esc(items.length + ' representations · ' +
      report.configurations_considered + ' configurations inspected · t-degree ≤ ' +
      report.search.max_var_deg + ', n-degree ≤ ' + report.search.max_idx_deg) + '</div>';
  if (items.length) {
    html += '<div class="representation-table"><table><thead><tr>' +
      '<th>Choose</th><th>Length</th><th>Derivative</th><th>t / n degree</th><th>Checked n</th>' +
      '</tr></thead><tbody>' + items.map((item, i) => {
        const c = item.complexity;
        return '<tr aria-selected="' + (i === selected) + '"><td>' +
          '<button class="copy-btn" onclick="renderRecurrenceMenu(' + i + ')">' +
          (i === selected ? 'Selected' : 'View') + '</button> ' +
          (report.frontier.includes(i) ? 'Tradeoff' : 'Alternative') + '</td><td>' +
          c.rec_len + '</td><td>' + c.diff_deg + '</td><td>' + c.var_deg + ' / ' + c.idx_deg +
          '</td><td>' + item.first_checked_index + '–' + item.last_checked_index + '</td></tr>';
      }).join('') + '</tbody></table></div>';
    const item = items[selected];
    const c = item.complexity;
    window._lastRecPlain = item.recurrence;
    window._lastRecLatex = item.latex;
    window._lastRecJson = JSON.stringify(item.recurrence_json, null, 2);
    window._lastRecMathematica = '';
    window._lastRecSage = '';
    window._lastGeneratedRows = null;
    html += '<div class="recurrence-display" id="recurrence-katex"></div>' +
      '<div class="recurrence-copy">' +
      '<button class="copy-btn" onclick="copyRec(this,\'plain\')">Copy plaintext</button>' +
      '<button class="copy-btn" onclick="copyRec(this,\'latex\')">Copy LaTeX</button>' +
      '<button class="copy-btn" onclick="copyRec(this,\'json\')">Copy selected JSON</button></div>' +
      '<div class="recurrence-note">' + esc(
        'Rows are indexed from 1. ' + item.fit_polynomials + ' fit rows; ' +
        item.verification_polynomials + ' held-out rows. Denominator t / n degrees: ' +
        c.denom_var_deg + ' / ' + c.denom_idx_deg + '.'
      ) + '</div>';
    html += item.singular_checked_indices.length
      ? '<div class="recurrence-note">Row generation is undefined at checked indices ' +
        esc(item.singular_checked_indices.join(', ')) + ' where the denominator vanishes.</div>'
      : renderRecurrenceExtension(report.source_rows);
  } else {
    html += '<p>No candidate returned. Try more rows or different bounds.</p>';
  }
  html += '<details><summary>Search coverage and full JSON</summary>' +
    '<p class="recurrence-note">“Not found” is a search outcome, not a nonexistence proof. ' +
    'Insufficient data and fitting heuristics can also prevent a candidate.</p>' +
    '<pre class="recurrence-code">' + esc(coverage) + '</pre>' +
    '<button class="copy-btn" onclick="copyRecurrenceMenu(this)">Copy comparison JSON</button></details>';
  setResult(html, 'recurrence');
  if (items.length) {
    const slot = document.getElementById('recurrence-katex');
    try { katex.render(items[selected].latex, slot, {displayMode: true, throwOnError: false}); }
    catch (_) { slot.textContent = items[selected].recurrence; }
  }
}

function copyRecurrenceMenu(button) {
  navigator.clipboard.writeText(JSON.stringify(window._recurrenceMenu, null, 2)).then(() => {
    button.textContent = 'Copied!';
  });
}

function recurrenceSpaceExample(kind) {
  if (kind === 'default') return {
    source: {coefficients: [['0', '0', '-1'], ['0'], ['0', '2'], ['0']],
      seeds: [['1'], ['2'], ['0', '2'], ['0', '3']]},
    shapes: [
      {kind: 'box', length: 3, index_degree: 1, variable_degree: 3, leading_index_degree: 1},
      {kind: 'box', length: 2, derivative_order: 1, index_degree: 1, variable_degree: 3, leading_index_degree: 1},
      {kind: 'box', length: 1, derivative_order: 2, index_degree: 1, variable_degree: 3, leading_index_degree: 1},
    ],
  };
  if (kind === 'hermite') return {
    indexed_source: {coefficients: [[['-2'], ['-2']], [['0', '2']]], seeds: [['1'], ['0', '2']]},
    shapes: [{kind: 'box', length: 1, derivative_order: 1, variable_degree: 1}],
  };
  if (kind === 'automatic') return {
    indexed_source: {coefficients: [[['0', '1'], ['0', '1']]], seeds: [['1']]},
    shapes: [{kind: 'automatic_variable_degree', length: 1, index_degree: 1}],
  };
  if (kind === 'binomial') return {
    source: {coefficients: [['1', '1']], seeds: [['1']], origin: 0},
    coefficient_rescaling: {kind: 'factorial_power', power: 1},
    shapes: [
      {kind: 'box', length: 1, derivative_order: 1, index_degree: 1, variable_degree: 1, leading_index_degree: 1},
      {kind: 'box', length: 2, index_degree: 1, variable_degree: 1, leading_index_degree: 1},
    ],
  };
  throw new Error('Unknown recurrence example');
}

function invalidateRecurrenceSpaces() {
  if (recurrenceRunning) cancelRecurrence();
  window._recurrenceSpaces = null;
  setResult('Settings changed. Find recurrences to update the results.', 'recurrence');
  setStatus('Previous result cleared.');
}

function spaceMainInputChanged() {
  if (document.querySelector('input[name="searchMode"][value="spaces"]').checked &&
      document.getElementById('spaceSource').value === 'rows') {
    invalidateRecurrenceSpaces();
    refreshSpaceRequestPreview();
  }
}

function readSpaceInteger(id, min, max) {
  const field = document.getElementById(id);
  const value = Number(field.value);
  if (!field.value.trim() || !Number.isInteger(value) || value < min || value > max) {
    throw new Error('Choose a whole number from ' + min + ' to ' + max + ' for ' +
      document.querySelector('label[for="' + id + '"]').textContent.toLowerCase() + '.');
  }
  return value;
}

// Sources with a constant-coefficient defining recurrence admit factorial rescaling.
const SPACE_RESCALING_SOURCES = ['default', 'binomial'];

function buildGuidedSpaceRequest() {
  const mode = document.getElementById('spaceSource').value;
  const request = mode === 'rows' ? {} : recurrenceSpaceExample(mode);
  if (mode === 'rows') {
    const input = getInput();
    if (input.length > 1048576) throw new Error('Polynomial input is too large (limit: 1 MiB).');
    const parsed = JSON.parse(wasm_bindgen.parse_and_format(input));
    const errors = parsed.filter(row => row.error).map(row => row.error);
    if (errors.length) throw new Error(errors.join('; '));
    if (!parsed.length) throw new Error('Enter polynomial rows above, or choose an example.');
    request.rows = parsed.map(row => row.coefficients);
    request.row_origin = readSpaceInteger('spaceOrigin', 0, 1000000);
    request.domain_start = request.row_origin;
  }
  delete request.coefficient_rescaling;
  const rescaling = SPACE_RESCALING_SOURCES.includes(mode) ? readSpaceInteger('spaceRescaling', 0, 2) : 0;
  if (rescaling) request.coefficient_rescaling = {kind: 'factorial_power', power: rescaling};
  const length = readSpaceInteger('spaceLength', 1, 4);
  const derivative = activeWorkflow === 'integer-sequence' ? 0 : readSpaceInteger('spaceDerivative', 0, 3);
  const automatic = activeWorkflow !== 'integer-sequence' && mode !== 'rows' && document.getElementById('spaceAutomatic').checked;
  request.shapes = [];
  for (let l = 1; l <= length; l++) {
    for (let d = 0; d <= derivative; d++) {
      const shape = {kind: automatic ? 'automatic_variable_degree' : 'box',
        length: l, derivative_order: d,
        index_degree: readSpaceInteger('spaceIndexDegree', 0, 3),
        leading_index_degree: readSpaceInteger('spaceLeaderDegree', 0, 2)};
      if (!automatic) shape.variable_degree = activeWorkflow === 'integer-sequence' ? 0 :
        readSpaceInteger('spaceVariableDegree', 0, 6);
      request.shapes.push(shape);
    }
  }
  return request;
}

function updateSpaceControls() {
  const mode = document.getElementById('spaceSource').value;
  const custom = mode === 'custom';
  const integer = activeWorkflow === 'integer-sequence';
  document.getElementById('spaceGuided').hidden = custom;
  document.getElementById('spaceOriginField').hidden = mode !== 'rows';
  document.getElementById('spaceRescalingField').hidden = !SPACE_RESCALING_SOURCES.includes(mode);
  document.getElementById('spaceAutomaticLabel').hidden = mode === 'rows' || integer;
  document.getElementById('spaceDerivative').disabled = integer;
  document.getElementById('spaceVariableDegree').disabled = integer ||
    (mode !== 'rows' && document.getElementById('spaceAutomatic').checked);
  if (custom) document.getElementById('spaceAdvanced').open = true;
  const help = {
    rows: 'Matches are checked on the rows entered above. The first row index is used in every formula.',
    default: 'Uses A171608: Pₙ₊₄ = 2xPₙ₊₂ − x²Pₙ, starting with 1, 2, 2x, 3x at n = 0. The input above is not used.',
    hermite: 'Uses Hermite polynomials: Pₙ₊₂ = 2xPₙ₊₁ − 2(n + 1)Pₙ, with P₀ = 1 and P₁ = 2x. The input above is not used.',
    automatic: 'Uses Pₙ₊₁ = (n + 1)xPₙ, with P₀ = 1. The input above is not used.',
    binomial: 'Uses Pₙ₊₁ = (1 + x)Pₙ with P₀ = 1, so Pₙ = (1 + x)ⁿ, and searches Qₙ with the coefficient of xᵏ divided by k!. The input above is not used.',
    custom: 'Uses the request JSON below, including its own rows or defining recurrence. The input above is not used.',
  };
  document.getElementById('spaceSourceHelp').textContent = help[mode];
  const length = Number(document.getElementById('spaceLength').value);
  const derivative = integer ? 0 : Number(document.getElementById('spaceDerivative').value);
  const xDegree = integer ? '0' : document.getElementById('spaceVariableDegree').value;
  const automatic = !integer && mode !== 'rows' && document.getElementById('spaceAutomatic').checked;
  document.getElementById('spaceSearchSummary').textContent =
    length * (derivative + 1) + ' combinations: ' + (length === 1 ? '1 previous term' : '1–' + length + ' previous terms') + ', ' +
    (derivative ? 'derivatives up to order ' + derivative : 'no derivatives') + '. ' +
    (automatic ? 'The x-degree is chosen automatically' : 'Coefficient degrees: x ≤ ' + xDegree) +
    ', n ≤ ' + document.getElementById('spaceIndexDegree').value + '.';
}

function refreshSpaceRequestPreview() {
  if (document.getElementById('spaceSource').value === 'custom' || !wasmReady) return;
  try {
    document.getElementById('spaceRequest').value = JSON.stringify(buildGuidedSpaceRequest(), null, 2);
  } catch (_) {
    // Never leave a previous example in the editor when the current input is invalid.
    document.getElementById('spaceRequest').value = '';
  }
}

function spaceControlsChanged() {
  updateSpaceControls();
  refreshSpaceRequestPreview();
  invalidateRecurrenceSpaces();
}

function changeSpaceSource() {
  const mode = document.getElementById('spaceSource').value;
  if (mode !== 'rows' && mode !== 'custom') {
    const shape = recurrenceSpaceExample(mode).shapes[0];
    document.getElementById('spaceLength').value = mode === 'default' ? 3 : (mode === 'binomial' ? 2 : shape.length);
    document.getElementById('spaceDerivative').value = mode === 'default' ? 2 : (shape.derivative_order || 0);
    document.getElementById('spaceIndexDegree').value = shape.index_degree || 0;
    document.getElementById('spaceVariableDegree').value = shape.variable_degree || 1;
    document.getElementById('spaceLeaderDegree').value = shape.leading_index_degree || 0;
    document.getElementById('spaceAutomatic').checked = mode === 'automatic';
    const rescaling = recurrenceSpaceExample(mode).coefficient_rescaling;
    document.getElementById('spaceRescaling').value = rescaling ? rescaling.power : 0;
  } else if (mode === 'rows') {
    document.getElementById('spaceAutomatic').checked = false;
  }
  spaceControlsChanged();
}

function editSpaceRequest() {
  document.getElementById('spaceSource').value = 'custom';
  updateSpaceControls();
  invalidateRecurrenceSpaces();
}

function prepareSpaceRequest() {
  if (document.getElementById('spaceSource').value === 'custom') {
    return document.getElementById('spaceRequest').value;
  }
  const request = JSON.stringify(buildGuidedSpaceRequest());
  document.getElementById('spaceRequest').value = JSON.stringify(JSON.parse(request), null, 2);
  return request;
}

function recurrenceSpaceSummary(report) {
  const total = report.results.length;
  const finished = report.results.filter(result => result.certificate &&
    ['successor_exists', 'absent_in_shape'].includes(result.status)).length;
  return finished + ' of ' + total + ' searches finished' +
    (finished === total ? '.' : '; ' + (total - finished) + ' unfinished.');
}

function spaceFormulaHtml(formula) {
  // Format the backend's forward-index notation; coefficients and indices are unchanged.
  const tex = formula.replace(/ \+ \(-1\) /g, ' - ').replace(/\(-1\) /g, '-')
    .replace(/\(1\) /g, '')
    .replace(/D\^(\d+)\(([PQ])\[n\+(\d+)\]\)/g, (_, order, name, lag) =>
      name + '_{' + (lag === '0' ? 'n' : 'n+' + lag) + '}' +
      (Number(order) <= 3 ? "'".repeat(Number(order)) : '^{(' + order + ')}'))
    .replace(/([PQ])\[n\+(\d+)\]/g, (_, name, lag) => name + '_{' + (lag === '0' ? 'n' : 'n+' + lag) + '}')
    .replace(/X/g, 'x').replace(/\^(\d+)/g, '^{$1}')
    .replace(/(\d+)\/(\d+)/g, '\\frac{$1}{$2}');
  return '<div class="recurrence-display">' + katexInline(tex) + '</div>';
}

function renderRecurrenceSpaces(report) {
  window._recurrenceSpaces = report;
  const source = report.request.source || report.request.indexed_source || report.request.source_recurrence;
  const states = {successor_exists: 'Recurrence found', absent_in_shape: 'No recurrence within these bounds',
    incomplete_budget: 'Search limit reached', insufficient_rows: 'More rows needed'};
  let html = '<strong>Recurrence search results</strong><p>' + esc(recurrenceSpaceSummary(report)) + '</p>';
  html += '<p class="space-help">' + (source
    ? 'Identities below follow from the defining recurrence and initial values.'
    : 'Matches below fit the supplied rows exactly; they may fail on later rows.') + '</p>';
  const rescaling = report.request.coefficient_rescaling;
  if (rescaling) html += '<p class="space-help">Formulas name the rescaled rows Q: the coefficient of x^k in P is divided by (k!)^' +
    esc(rescaling.power) + '.</p>';
  const ordered = report.results.map((result, i) => ({result, i}));
  ordered.sort((a, b) => Number(b.result.status === 'successor_exists') - Number(a.result.status === 'successor_exists'));
  ordered.forEach(({result, i}) => {
    const shape = report.request.shapes && report.request.shapes[i];
    const name = shape && shape.length ? shape.length + ' previous term' + (shape.length === 1 ? '' : 's') +
      ' · ' + ((shape.derivative_order || 0) === 0 ? 'no derivatives' : 'derivatives up to order ' + shape.derivative_order) :
      'Search ' + (i + 1);
    html += '<section class="space-result"><h3>' + esc(name) + '</h3><p><strong>' +
      esc(states[result.status] || 'Search unfinished') + '</strong></p>';
    if (result.reason) html += '<p class="space-help">' + esc(result.reason) + '</p>';
    if (result.status === 'absent_in_shape') html += '<p class="space-help">Try higher coefficient bounds or allow more terms or derivatives.</p>';
    if (!['successor_exists', 'absent_in_shape'].includes(result.status)) {
      html += '<p class="space-help">This search does not tell us whether a recurrence exists.</p>';
    }
    result.representatives.forEach((rep, j) => {
      if (j === 1) html += '<details><summary>More formulas (' + (result.representatives.length - 1) + ')</summary>';
      html += spaceFormulaHtml(rep.formula) +
        '<button class="copy-btn" onclick="copySpaceFormula(this,' + i + ',' + j + ')">Copy formula</button>';
      if (rep.singular_ranks === null) html += '<p class="space-help">Check when the coefficient of the next term vanishes before using this rule to generate rows.</p>';
      else if (rep.singular_ranks.length) html += '<p class="space-help">Cannot determine the next term at base n = ' +
        esc(rep.singular_ranks.join(', ')) + '; supply these rows separately.</p>';
      if (rep.requires_x_division) html += '<p class="space-help">Generating the next row requires exact polynomial division.</p>';
    });
    if (result.representatives.length > 1) html += '</details>';
    html += '<details><summary>Search bounds</summary><p class="space-help">Base index starts at n = ' +
      esc(report.request.domain_start) + '.';
    if (shape && shape.kind !== 'explicit') {
      html += ' Coefficient degree in n: ' + esc(shape.index_degree || 0) +
        '; in x: ' + esc(result.automatic_degree_bound ? result.automatic_degree_bound.variable_degree : shape.variable_degree || 0) +
        '. Next-term factor degrees: n ≤ ' + esc(shape.leading_index_degree || 0) + ', x ≤ ' + esc(shape.leading_variable_degree || 0) + '.';
    }
    if (result.rescaling_bound) html += ' Checking ' + esc(result.rescaling_bound.annihilator_power) +
      ' base indices covers every n for the rescaled rows.';
    if (result.automatic_degree_bound) html += result.certificate
      ? ' The x-degree was chosen automatically to cover all possible x-degrees for existence within the other bounds.'
      : ' The automatic x-degree search is unfinished.';
    html += '</p></details></section>';
  });
  setResult(html, 'recurrence');
}

async function copySpaceFormula(button, resultIndex, formulaIndex) {
  const formula = window._recurrenceSpaces?.results[resultIndex]?.representatives[formulaIndex]?.formula;
  if (!formula) return;
  try {
    await navigator.clipboard.writeText(formula);
    button.textContent = 'Copied!';
    setTimeout(() => { button.textContent = 'Copy formula'; }, 1500);
  } catch (_) {
    setStatus('Could not copy the formula. Check clipboard permission and try again.');
  }
}

function createWorker() {
  const base = new URL('.', location.href).href;
  const code = `
    importScripts('${base}pkg/polytool_web.js?v=${ASSET_VERSION}');
    wasm_bindgen({ module_or_path: '${base}pkg/polytool_web_bg.wasm?v=${ASSET_VERSION}' }).then(() => {
      postMessage({type: 'ready'});
    });
    onmessage = function(e) {
      if (e.data.type === 'find_recurrence') {
        if (e.data.spaces) {
          const result = wasm_bindgen.find_recurrence_spaces(e.data.input);
          postMessage({type: 'result', data: result});
          return;
        }
        const finder = e.data.representations
          ? wasm_bindgen.find_recurrence_representations : wasm_bindgen.find_recurrence;
        const result = finder(
          e.data.input, e.data.max_rec_len, e.data.max_var_deg,
          e.data.max_idx_deg, e.data.max_diff_deg,
          e.data.try_inhomogeneous, e.data.try_denominator,
          e.data.try_alternating_sign, e.data.max_candidates,
          e.data.max_denom_var_deg, e.data.max_denom_idx_deg
        );
        postMessage({type: 'result', data: result});
      }
    };
  `;
  const blob = new Blob([code], {type: 'application/javascript'});
  return new Worker(URL.createObjectURL(blob));
}

function startRecurrence(text) {
  if (recurrenceRunning) return;
  recurrenceRunning = true;
  window._recurrenceSpaces = null;
  document.getElementById('btnRec').disabled = true;
  document.getElementById('btnCancel').style.display = '';
  setResult('Searching for recurrence\u2026', 'recurrence');
  setStatus('Computing (this may take a while for large search spaces)\u2026');

  const adaptive = document.querySelector('input[name="searchMode"][value="adaptive"]').checked;
  const readInt = (id, fallback) => {
    const value = parseInt(document.getElementById(id).value, 10);
    return Number.isNaN(value) ? fallback : value;
  };
  const params = {
    type: 'find_recurrence',
    spaces: document.querySelector('input[name="searchMode"][value="spaces"]').checked,
    representations: document.querySelector('input[name="searchMode"][value="representations"]').checked,
    max_candidates: readInt('representationBudget', 2000),
    max_denom_var_deg: readInt('representationDenomVar', 0),
    max_denom_idx_deg: readInt('representationDenomIdx', 1),
    input: text,
    max_rec_len: adaptive ? 10 : readInt('maxRecLen', 3),
    max_var_deg: activeWorkflow === 'integer-sequence'
      ? 0 : (adaptive ? 5 : readInt('maxVarDeg', 2)),
    max_idx_deg: adaptive ? 5 : readInt('maxIdxDeg', 2),
    max_diff_deg: activeWorkflow === 'integer-sequence'
      ? 0 : (adaptive ? 5 : readInt('maxDiffDeg', 1)),
    try_inhomogeneous: document.getElementById('tryInhomogeneous').checked,
    try_denominator: document.getElementById('tryDenominator').checked,
    try_alternating_sign: document.getElementById('tryAlternatingSign').checked,
  };
  const sourceRows = text.split(/\r?\n/).filter(line => line.trim()).length;

  if (worker) worker.terminate();
  worker = createWorker();
  const t0 = performance.now();

  worker.onmessage = function(e) {
    if (e.data.type === 'ready') {
      worker.postMessage(params);
    } else if (e.data.type === 'result') {
      const ms = (performance.now() - t0).toFixed(0);
      recurrenceRunning = false;
      document.getElementById('btnRec').disabled = false;
      document.getElementById('btnCancel').style.display = 'none';

      try {
        const data = JSON.parse(e.data.data);
        if (data.schema === 'polytool.recurrence-spaces.v1') {
          renderRecurrenceSpaces(data);
          setStatus(recurrenceSpaceSummary(data) + ' Returned in ' + ms + ' ms.');
        } else if (data.schema === 'polytool.recurrence-representations.v1') {
          window._recurrenceMenu = data;
          renderRecurrenceMenu(data.frontier[0] ?? 0);
          setStatus((data.status === 'budget_exhausted' ? 'Partial comparison' : 'Comparison completed') + ' in ' + ms + ' ms.');
        } else if (data.found) {
          const domainNote = activeWorkflow === 'integer-sequence'
            ? '<div class="recurrence-note">Here <strong>P(n)</strong> denotes the entered integer sequence; the search is restricted to scalar recurrences in n.</div>'
            : '';
          const latexEl = '<div class="recurrence-display" id="recurrence-katex"></div>';
          const copyEl = '<div class="recurrence-copy">' +
            '<button class="copy-btn" onclick="copyRec(this,\'plain\')">Copy plaintext</button>' +
            '<button class="copy-btn" onclick="copyRec(this,\'latex\')">Copy LaTeX</button>' +
            (data.mathematica ? '<button class="copy-btn" onclick="copyRec(this,\'mathematica\')">Copy Mathematica</button>' : '') +
            (data.sage ? '<button class="copy-btn" onclick="copyRec(this,\'sage\')">Copy Sage</button>' : '') +
            (data.recurrence_json ? '<button class="copy-btn" onclick="copyRec(this,\'json\')">Copy JSON</button>' : '') +
            (data.recurrence_json && activeWorkflow === 'polynomial-sequence'
              ? '<button class="copy-btn" onclick="openTransformWithFoundRecurrence()">Transform or check theorems</button>' : '') +
            '</div>';
          const statsEl = '<div class="recurrence-stats">' +
                          esc(data.unknowns + ' unknowns, ' +
                          data.equations + ' equations, ' +
                          data.candidates_tried + ' candidates tried') + '</div>';
          const codeEl = renderRecurrenceCodeStack(data);
          setResult(domainNote + latexEl + copyEl + statsEl + renderRecurrenceExtension(sourceRows) + codeEl, 'recurrence');
          // Store for copy
          window._lastRecPlain = data.recurrence;
          window._lastRecLatex = data.latex;
          window._lastRecMathematica = data.mathematica || '';
          window._lastRecSage = data.sage || '';
          window._lastRecJson = data.recurrence_json || '';
          // Render KaTeX (may load async)
          try {
            katex.render(data.latex, document.getElementById('recurrence-katex'), {displayMode: true, throwOnError: false});
          } catch(e) {
            document.getElementById('recurrence-katex').textContent = data.recurrence;
          }
          setStatus('Found in ' + ms + ' ms.');
        } else {
          setResult(span(data.error ? 'err' : 'no', data.error || 'No recurrence found.'), 'recurrence');
          setStatus((data.error ? 'Request rejected' : 'Search completed') + ' in ' + ms + ' ms.');
        }
      } catch (error) {
        setResult(span('err', 'Could not display recurrence result: ' + error), 'recurrence');
        setStatus('Error displaying recurrence result.');
      }
    }
  };

  worker.onerror = function(e) {
    recurrenceRunning = false;
    document.getElementById('btnRec').disabled = false;
    document.getElementById('btnCancel').style.display = 'none';
    setResult(span('err', 'Worker error: ' + e.message), 'recurrence');
    setStatus('Error.');
  };
}

function cancelRecurrence() {
  if (worker) { worker.terminate(); worker = null; }
  recurrenceRunning = false;
  document.getElementById('btnRec').disabled = false;
  document.getElementById('btnCancel').style.display = 'none';
  setResult('Recurrence search cancelled.', 'recurrence');
  setStatus('Cancelled.');
}
