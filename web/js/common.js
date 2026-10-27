function loadExample(name) {
  const ex = EXAMPLES[name] || SEQUENCE_CATALOG[name];
  if (!ex) return;
  document.getElementById('input').value = ex.lines.join('\n');
  spaceMainInputChanged();
  setStatus('Loaded: ' + ex.label);
}

function loadWorkflowExample(workflow, name) {
  const text = WORKFLOW_EXAMPLES[workflow] && WORKFLOW_EXAMPLES[workflow][name];
  if (!text) return;
  document.getElementById('input').value = text;
  workflowDrafts[workflow] = text;
  spaceMainInputChanged();
  setStatus('Example loaded.');
}

function switchTab(name, trigger) {
  const panel = document.getElementById('panel-' + name);
  const tab = trigger || document.querySelector('.tab[data-tab="' + name + '"]');
  if (!panel || !tab || tab.hidden) return;
  document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tab-panel').forEach(p => p.classList.remove('active'));
  tab.classList.add('active');
  panel.classList.add('active');
  activeTab = name;
}

function setWorkflow(name, options) {
  if (!WORKFLOW_CONFIG[name]) return;
  const opts = options || {};
  const input = document.getElementById('input');
  if (activeWorkflow && input) workflowDrafts[activeWorkflow] = input.value;
  if (recurrenceRunning) cancelRecurrence();
  if (transformRunning) cancelTransformJob(true);
  transformEntries = [];
  activeWorkflow = name;
  const config = WORKFLOW_CONFIG[name];

  document.querySelectorAll('.workflow-button').forEach(button => {
    button.classList.toggle('active', button.dataset.workflow === name);
  });
  document.querySelectorAll('.workflow-examples').forEach(group => { group.hidden = true; });
  const examplesId = name === 'polynomial' ? 'singlePolynomialExamples' :
    (name === 'polynomial-sequence' ? 'polynomialSequenceExamples' : 'integerSequenceExamples');
  document.getElementById(examplesId).hidden = false;
  document.getElementById('inputLabel').textContent = config.label;
  document.getElementById('inputHelp').textContent = config.help;
  input.placeholder = config.placeholder;
  if (!opts.keepInput) input.value = workflowDrafts[name] || '';

  document.querySelectorAll('.tab').forEach(tab => {
    tab.hidden = !tab.dataset.workflows.split(/\s+/).includes(name);
  });
  document.querySelectorAll('.tab-panel').forEach(panel => {
    const tab = document.querySelector('.tab[data-tab="' + panel.id.replace('panel-', '') + '"]');
    panel.hidden = !tab || tab.hidden;
  });
  document.getElementById('btnProperties').textContent =
    name === 'polynomial' ? 'Examine polynomial' : 'Check sequence properties';
  document.querySelectorAll('.polynomial-recurrence-option').forEach(option => {
    option.hidden = name === 'integer-sequence';
  });
  updateSpaceControls();
  refreshSpaceRequestPreview();
  window._recurrenceSpaces = null;
  switchTab(config.defaultTab);
  document.querySelectorAll('.tab-results').forEach(result => { result.innerHTML = ''; });
  setStatus('');
  if (!opts.skipHash && history.replaceState) history.replaceState(null, '', '#' + name);
}

function toggleSearchMode() {
  if (recurrenceRunning) cancelRecurrence();
  const adaptive = document.querySelector('input[name="searchMode"][value="adaptive"]').checked;
  const compare = document.querySelector('input[name="searchMode"][value="representations"]').checked;
  const spaces = document.querySelector('input[name="searchMode"][value="spaces"]').checked;
  document.getElementById('spaceOptions').hidden = !spaces;
  document.getElementById('legacyRecurrenceOptions').hidden = spaces;
  document.querySelectorAll('#legacyRecurrenceOptions input').forEach(input => { input.disabled = spaces; });
  document.getElementById('manualOpts').hidden = spaces;
  document.getElementById('manualOpts').classList.toggle('disabled', adaptive);
  updateSpaceControls();
  refreshSpaceRequestPreview();
  document.getElementById('comparisonOptions').hidden = !compare;
  document.getElementById('btnRec').textContent = spaces ? 'Find recurrences' : (compare ? 'Compare representations' : 'Find recurrence');
}

let wasmReady = false;
let worker = null;
let recurrenceRunning = false;

let check_properties, check_interlacing_pairs, check_interlacing_pair, find_recurrence,
    generate_recurrence_rows, analyze_decompositions;

async function init() {
  try {
    await wasm_bindgen({ module_or_path: 'pkg/polytool_web_bg.wasm?v=' + ASSET_VERSION });
    ({ check_properties, check_interlacing_pairs, check_interlacing_pair, find_recurrence,
       generate_recurrence_rows, analyze_decompositions } = wasm_bindgen);
    wasmReady = true;
    populateTransformCatalog();
    refreshTransformRequestPreview();
    setStatus('');
  } catch (e) {
    setStatus('Failed to load WASM: ' + e);
  }
}

function setStatus(msg) {
  document.getElementById('status').textContent = msg;
}

let activeTab = 'properties';
const OEIS_PROXY_URL = 'https://arxiv.symmetricfunctions.com/api/oeis-search';
let oeisLookupToken = 0;
let oeisInlineRequestToken = 0;
let oeisInlineSlotSerial = 0;

function setResult(html, tab = activeTab) {
  hideFloatingInfoTip();
  document.getElementById('results-' + tab).innerHTML = html;
}

function stripCommentOnlyLines(raw) {
  return String(raw || '')
    .split(/\r?\n/)
    .map(line => /^\s*#/.test(line) ? '' : line)
    .join('\n');
}

function parseCoefficientListLiteral(raw) {
  let s = stripCommentOnlyLines(raw).trim();
  const assignment = s.match(/^[A-Za-z$][A-Za-z0-9$]*\s*:?=\s*([\s\S]*)$/);
  if (assignment) s = assignment[1].trim();
  s = s.replace(/;\s*$/, '').trim();
  if (!s) return null;

  let pos = 0;
  const skipWs = () => {
    while (pos < s.length && /\s/.test(s[pos])) pos++;
  };

  const parseInteger = () => {
    skipWs();
    const start = pos;
    if (s[pos] === '+' || s[pos] === '-') pos++;
    const digitsStart = pos;
    while (pos < s.length && /\d/.test(s[pos])) pos++;
    if (pos === digitsStart) {
      pos = start;
      return null;
    }
    return s.slice(start, pos);
  };

  const parseList = (depth = 0) => {
    // A coefficient row or a list of rows needs at most two levels.
    if (depth > 1) return null;
    skipWs();
    const open = s[pos];
    const close = ({ '{': '}', '[': ']', '(': ')' })[open];
    if (!close) return null;
    pos++;
    const items = [];

    while (true) {
      skipWs();
      if (s[pos] === close) {
        pos++;
        return items;
      }

      const item = /[\{[(]/.test(s[pos] || '') ? parseList(depth + 1) : parseInteger();
      if (item === null) return null;
      items.push(item);

      skipWs();
      if (s[pos] === ',') {
        pos++;
        continue;
      }
      if (s[pos] === close) {
        pos++;
        return items;
      }
      return null;
    }
  };

  const parsed = parseList();
  if (parsed === null) return null;
  skipWs();
  if (pos !== s.length) return null;

  const isInteger = item => typeof item === 'string';
  const isIntegerRow = item => Array.isArray(item) && item.length > 0 && item.every(isInteger);
  if (parsed.length > 0 && parsed.every(isInteger)) return [parsed];
  if (parsed.length > 0 && parsed.every(isIntegerRow)) return parsed;
  return null;
}

function parseLinewiseCoefficientLists(raw) {
  // Split only outside brackets: a pasted row may wrap onto several lines.
  // Bare commas separate coefficients; commas after a closed list separate rows.
  const chunks = [];
  const stack = [];
  const closes = { '{': '}', '[': ']', '(': ')' };
  let start = 0;
  let previous = '';
  for (let i = 0; i < raw.length; i++) {
    const ch = raw[i];
    if (closes[ch]) {
      stack.push(closes[ch]);
      if (stack.length > 2) return null;
    } else if (ch === '}' || ch === ']' || ch === ')') {
      if (stack.pop() !== ch) return null;
    } else if (!stack.length && (ch === ';' || ch === '\n' ||
               (ch === ',' && /[}\])]/.test(previous)))) {
      chunks.push(raw.slice(start, i));
      start = i + 1;
    }
    if (!/\s/.test(ch)) previous = ch;
  }
  if (stack.length) return null;
  chunks.push(raw.slice(start));
  const rows = [];
  for (const chunk of chunks) {
    const trimmed = chunk.trim();
    if (!trimmed) continue;
    const parsed = parseCoefficientListLiteral(trimmed);
    if (parsed && parsed.length === 1) {
      rows.push(parsed[0]);
    } else {
      // Validate the entire row before tokenizing. Never discard bad tokens,
      // empty comma fields, decimals or ellipses, or round through JS Number.
      if (!/^[+-]?\d+(?:(?:\s*,\s*|\s+)[+-]?\d+)*\s*,?$/.test(trimmed)) return null;
      rows.push(trimmed.match(/[+-]?\d+/g));
    }
  }
  return rows.length ? rows : null;
}

function normalizePolynomialInput(raw) {
  const clean = stripCommentOnlyLines(String(raw || '')
    .replace(/\r\n?|[\u2028\u2029]/g, '\n')
    .replace(/\u2212/g, '-')
    .replace(/[\u00a0\u202f\ufeff]/g, ' '));
  const rows = parseCoefficientListLiteral(clean) || parseLinewiseCoefficientLists(clean);
  if (!rows) return clean; // Let the existing polynomial parser report the error.
  return rows.map(row => row.join(', ')).join('\n');
}

function parseIntegerSequenceInput(raw) {
  const clean = stripCommentOnlyLines(String(raw || '')
    .replace(/\r\n?|[\u2028\u2029]/g, '\n')
    .replace(/\u2212/g, '-')
    .replace(/[\u00a0\u202f\ufeff]/g, ' '))
    .trim();
  if (!clean) return [];

  const lines = clean.split('\n').map(line => line.trim()).filter(Boolean);
  const pairs = lines.map(line => line.match(/^([+-]?\d+)\s+([+-]?\d+)$/));
  if (lines.length >= 2 && pairs.every(Boolean)) {
    const indices = pairs.map(match => BigInt(match[1]));
    const consecutive = indices.every((value, index) =>
      index === 0 || value === indices[index - 1] + 1n);
    if (consecutive) return pairs.map(match => match[2]);
  }

  if (!/^[\s,;{}\[\]()+\-\d]+$/.test(clean)) return null;
  const stripped = clean.replace(/[{}\[\]()]/g, ' ');
  if (/[,;]\s*[,;]/.test(stripped) || /(?:^|[,;])\s*(?:$|\n)/m.test(stripped)) return null;
  const terms = stripped.match(/[+-]?\d+/g) || [];
  const residue = stripped.replace(/[+-]?\d+/g, '').replace(/[\s,;]/g, '');
  return residue ? null : terms;
}

function getIntegerTerms() {
  const terms = parseIntegerSequenceInput(document.getElementById('input').value);
  if (terms === null) {
    throw new Error('Integer sequences may contain only whole numbers and ordinary separators.');
  }
  return terms;
}

function getInput() {
  if (typeof activeWorkflow !== 'undefined' && activeWorkflow === 'integer-sequence') {
    return getIntegerTerms().join('\n');
  }
  return normalizePolynomialInput(document.getElementById('input').value);
}

function copyRec(btn, kind) {
  const key = ({
    plain: '_lastRecPlain',
    latex: '_lastRecLatex',
    mathematica: '_lastRecMathematica',
    sage: '_lastRecSage',
    json: '_lastRecJson',
  })[kind] || '_lastRecPlain';
  const text = window[key] || '';
  navigator.clipboard.writeText(text).then(() => {
    const orig = btn.textContent;
    btn.textContent = 'Copied!';
    setTimeout(() => btn.textContent = orig, 1500);
  });
}

function renderRecurrenceCodeCard(title, note, code, link) {
  if (!code) return '';
  const linkHtml = link
    ? ' <a href="' + escAttr(link.href) + '" target="_blank" rel="noopener">' +
      esc(link.label) + '</a>.'
    : '';
  return '<div class="recurrence-code-card">' +
    '<h3>' + esc(title) + '</h3>' +
    '<div class="recurrence-note">' + esc(note) + linkHtml + '</div>' +
    '<pre class="recurrence-code">' + esc(code) + '</pre>' +
    '</div>';
}

function renderRecurrenceCodeStack(data) {
  const cards = [];
  if (data.mathematica) {
    cards.push(renderRecurrenceCodeCard(
      'Mathematica',
      'Defines P[n, t] from the discovered recurrence together with the needed initial values.',
      data.mathematica
    ));
  }
  if (data.sage) {
    cards.push(renderRecurrenceCodeCard(
      'Sage',
      'Defines P(n) in a polynomial ring over QQ, memoized from the discovered recurrence and initial values.',
      data.sage
    ));
  }
  if (data.recurrence_json) {
    cards.push(renderRecurrenceCodeCard(
      'Recurrence JSON',
      'Machine-readable recurrence with the initial values needed to reproduce the sequence.',
      data.recurrence_json,
      {
        href: 'https://github.com/PerAlexandersson/polytool#machine-readable-recurrence-json-and-row-generation',
        label: 'Generate more rows with polytool recurrence-generate'
      }
    ));
  }
  return cards.length ? ('<div class="recurrence-code-stack">' + cards.join('') + '</div>') : '';
}

function renderRecurrenceExtension(sourceRows) {
  return '<div class="extension-card">' +
    '<h3>Extend this sequence</h3>' +
    '<div class="recurrence-note">Generate exact new terms from the recurrence, preview them, then move the extended sequence back into the input.</div>' +
    '<div class="extension-controls"><label>Add <input type="number" id="extensionCount" value="10" min="1" max="200"> rows</label>' +
    '<button class="run-small" onclick="previewRecurrenceExtension(' + sourceRows + ')">Generate preview</button></div>' +
    '<div id="extensionResult"></div></div>';
}

function previewRecurrenceExtension(sourceRows) {
  const slot = document.getElementById('extensionResult');
  const additional = parseInt(document.getElementById('extensionCount').value, 10);
  if (!Number.isInteger(additional) || additional < 1 || additional > 200) {
    slot.innerHTML = '<div class="err">Choose between 1 and 200 additional rows.</div>';
    return;
  }
  const totalRows = sourceRows + additional;
  const data = JSON.parse(generate_recurrence_rows(window._lastRecJson || '', totalRows));
  if (data.error) {
    slot.innerHTML = '<div class="err">' + esc(data.error) + '</div>';
    return;
  }
  window._lastGeneratedRows = data;
  window._lastGeneratedWorkflow = activeWorkflow;
  let previewRows = data.polynomials.slice(sourceRows);
  let preview = previewRows.map(row => row.join(', ')).join('\n');
  if (activeWorkflow === 'integer-sequence') preview = previewRows.map(row => row[0]).join(', ');
  const note = data.integral
    ? 'Generated ' + additional + ' new row' + (additional === 1 ? '' : 's') + ' exactly.'
    : 'The recurrence generated rational coefficients. They can be copied, but the current input accepts integer coefficients only.';
  slot.innerHTML = '<div class="recurrence-note">' + esc(note) + '</div>' +
    '<pre class="recurrence-code extension-preview">' + esc(preview) + '</pre>' +
    '<div class="extension-controls">' +
      (data.integral ? '<button class="run-small" onclick="useGeneratedRows()">Use extended sequence</button>' : '') +
      '<button class="run-small" onclick="copyGeneratedRows(this)">Copy full sequence</button></div>';
}

function generatedRowsText(data, workflow) {
  if (workflow === 'integer-sequence') return data.polynomials.map(row => row[0]).join(', ');
  return data.rows_text;
}

function useGeneratedRows() {
  const data = window._lastGeneratedRows;
  if (!data || !data.integral) return;
  const text = generatedRowsText(data, window._lastGeneratedWorkflow);
  document.getElementById('input').value = text;
  workflowDrafts[activeWorkflow] = text;
  setStatus('Extended sequence moved into the input; you can analyze it or search again.');
}

function copyGeneratedRows(button) {
  const data = window._lastGeneratedRows;
  if (!data) return;
  navigator.clipboard.writeText(generatedRowsText(data, window._lastGeneratedWorkflow)).then(() => {
    const label = button.textContent;
    button.textContent = 'Copied!';
    setTimeout(() => { button.textContent = label; }, 1500);
  });
}

function span(cls, text) {
  return '<span class="' + cls + '">' + esc(text) + '</span>';
}

function esc(s) {
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

function escAttr(s) {
  return esc(String(s)).replace(/"/g, '&quot;').replace(/'/g, '&#39;');
}

function katexInline(tex) {
  try {
    if (window.katex && typeof window.katex.renderToString === 'function') {
      return window.katex.renderToString(tex, { displayMode: false, throwOnError: false });
    }
  } catch (e) {}
  return '<span>' + esc(tex) + '</span>';
}

async function run(action) {
  if (!wasmReady) { setStatus('WASM not loaded yet.'); return; }
  try {
    const spaces = action === 'recurrence' && document.querySelector('input[name="searchMode"][value="spaces"]').checked;
    const text = spaces ? prepareSpaceRequest() : getInput();
    if (!text.trim()) { setStatus('No input.'); return; }
    switch (action) {
      case 'properties': renderProperties(text); break;
      case 'decomposition': renderDecompositions(text); break;
      case 'recurrence': startRecurrence(text); break;
      case 'oeis': await runOeis(); break;
    }
  } catch (e) {
    setResult(span('err', 'Error: ' + e), action === 'recurrence' ? 'recurrence' : undefined);
    if (action === 'recurrence') setStatus('Check the input and search settings.');
  }
}
