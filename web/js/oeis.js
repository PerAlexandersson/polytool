// --- OEIS lookup ---

function isZeroValue(value) {
  return normalizeIntegerString(value) === '0';
}

function stripLeadingZeros(seq) {
  let i = 0;
  while (i < seq.length && isZeroValue(seq[i])) i++;
  if (i >= seq.length) return seq.length ? ['0'] : [];
  return seq.slice(i);
}

function firstDifferences(seq) {
  const integers = seq.map(normalizeIntegerString);
  if (integers.some(value => value === null)) return [];
  const out = [];
  for (let i = 1; i < integers.length; i++) {
    out.push((BigInt(integers[i]) - BigInt(integers[i - 1])).toString());
  }
  return out;
}

function parseSeqData(raw) {
  return String(raw || '')
    .split(',')
    .map(s => s.trim())
    .filter(Boolean);
}

function oeisSearchHref(seq) {
  return 'https://oeis.org/search?q=' + encodeURIComponent(seq.join(',')) + '&language=english';
}

function flattenRows(rows) {
  return rows.flatMap(row => row && row.length ? row : ['0']);
}

function sumIntegerRow(row, alternating) {
  const safeRow = row && row.length ? row : ['0'];
  let sum = 0n;
  for (let i = 0; i < safeRow.length; i++) {
    const term = normalizeIntegerString(safeRow[i]);
    if (term === null) return null;
    const value = BigInt(term);
    sum += alternating && i % 2 === 1 ? -value : value;
  }
  return sum.toString();
}

function rowSumsOfRows(rows) {
  const sums = rows.map(row => sumIntegerRow(row, false));
  return sums.some(value => value === null) ? null : sums;
}

function copyData(btn) {
  const text = btn.dataset.copy || '';
  const htmlLabel = btn.dataset.labelHtml || '';
  navigator.clipboard.writeText(text).then(() => {
    const orig = htmlLabel ? btn.innerHTML : btn.textContent;
    btn.textContent = 'Copied!';
    setTimeout(() => {
      if (htmlLabel) {
        btn.innerHTML = htmlLabel;
      } else {
        btn.textContent = orig;
      }
    }, 1500);
  });
}

function queryLink(label, seq) {
  if (!seq.length) return '';
  return '<a class="oeis-link" href="' + oeisSearchHref(seq) + '" target="_blank" rel="noopener" title="' +
    escAttr(seq.join(', ')) + '">' + esc(label) + '</a>';
}

function queryButton(label, seq, slotId, queryLabel, queryLabelHtml) {
  if (!seq.length || !slotId) return '';
  return '<button class="run-small" data-oeis-seq="' + escAttr(seq.join(',')) +
    '" data-oeis-slot="' + escAttr(slotId) +
    '" data-oeis-label="' + escAttr(queryLabel || label) +
    '" data-oeis-label-html="' + escAttr(queryLabelHtml || '') +
    '" onclick="runOeisInlineQuery(this)">' + esc(label) + '</button>';
}

function columnCopyButton(labelTex, rows, options) {
  const usableRows = rows.filter(row => row != null);
  if (!usableRows.length) return '';
  const text = usableRows.map(serializeCoeffRow).join('\n');
  const labelHtml = 'Copy ' + katexInline(labelTex);
  const skipped = rows.length - usableRows.length;
  const baseTitle = (options && options.title) || ('Copy ' + usableRows.length + ' coefficient row(s) as newline-separated input.');
  const title = skipped > 0 ? (baseTitle + ' Unavailable rows are skipped.') : baseTitle;
  return '<button class="run-small" data-copy="' + escAttr(text) +
    '" data-label-html="' + escAttr(labelHtml) +
    '" onclick="copyData(this)" title="' + escAttr(title) + '">' + labelHtml + '</button>';
}

function renderDecompCopyBar(view, data) {
  const clean = data.filter(item => !item.error);
  if (!clean.length) return '';

  const omitted = data.length - clean.length;
  const buttons = [];

  if (view === 'core') {
    buttons.push(columnCopyButton('p(t)', clean.map(item => item.coefficients)));
    buttons.push(columnCopyButton('a(t)', clean.map(item => item.a.coefficients)));
    buttons.push(columnCopyButton('\\gamma_a(z)', clean.map(item => item.gamma_a_coefficients)));
    buttons.push(columnCopyButton('b(t)', clean.map(item => item.b.coefficients)));
    buttons.push(columnCopyButton('\\gamma_b(z)', clean.map(item => item.gamma_b_coefficients)));
  } else if (view === 'companion') {
    buttons.push(columnCopyButton('f(t)', clean.map(item => item.f_polynomial.coefficients)));
    buttons.push(columnCopyButton('R_d(f)', clean.map(item => item.r_transform_of_f.coefficients)));
    buttons.push(columnCopyButton('\\widetilde{a}(t)', clean.map(item => item.r_a.coefficients)));
    buttons.push(columnCopyButton('\\widetilde{b}(t)', clean.map(item => item.r_b.coefficients)));
  } else {
    buttons.push(columnCopyButton('p(t)', clean.map(item => item.coefficients)));
    buttons.push(columnCopyButton('c_j', clean.map(item => item.magic.coordinates)));
  }

  const note = 'Copy a displayed column as newline-separated coefficient vectors; you can paste the result back into the main input box.' +
    (omitted ? ' Rows with analysis errors are skipped.' : '');
  return '<div class="decomp-copybar">' +
    '<div class="oeis-note">' + note + '</div>' +
    '<div class="oeis-actions">' + buttons.filter(Boolean).join('') + '</div>' +
    '</div>';
}

function oeisDataPreview(data) {
  const terms = String(data || '').split(',').map(s => s.trim()).filter(Boolean);
  if (!terms.length) return '';
  const preview = terms.slice(0, 18).join(', ');
  return preview + (terms.length > 18 ? ', ...' : '');
}

function renderOeisMatches(payload) {
  const results = Array.isArray(payload && payload.results) ? payload.results : [];
  if (!results.length) {
    return '<div class="oeis-none">No automatic match returned for this query yet. Try the OEIS link above.</div>';
  }

  let html = '';
  for (const item of results) {
    const id = item.id || (item.number != null ? ('A' + String(item.number).padStart(6, '0')) : 'OEIS');
    const href = item.id ? ('https://oeis.org/' + encodeURIComponent(item.id)) : null;
    html += '<div class="oeis-match">';
    html += href
      ? '<a href="' + href + '" target="_blank" rel="noopener">' + esc(id) + '</a>'
      : '<strong>' + esc(id) + '</strong>';
    if (item.name) {
      html += '<span class="oeis-name">' + esc(item.name) + '</span>';
    }
    const preview = oeisDataPreview(item.data);
    if (preview) {
      html += '<div class="oeis-seq oeis-data">' + esc(preview) + '</div>';
    }
    html += '</div>';
  }
  if (payload.count > results.length) {
    html += '<div class="oeis-none">Showing the first ' + results.length + ' of ' + payload.count + ' OEIS matches.</div>';
  }
  return html;
}

function renderOeisMatchesWithContext(seq, label, labelHtml, payload) {
  let html = '<div class="oeis-note">Showing matches for <strong>' + (labelHtml || esc(label)) +
    '</strong>. <a href="' + oeisSearchHref(seq) + '" target="_blank" rel="noopener">Open this query on OEIS</a>.</div>';
  html += renderOeisMatches(payload);
  return html;
}

async function fetchOeisMatches(seq, slotId, token) {
  const slot = document.getElementById(slotId);
  if (!slot || !seq.length) return { ok: false, skipped: true };

  slot.innerHTML = '<div class="oeis-none">Looking up OEIS...</div>';

  try {
    const resp = await fetch(
      OEIS_PROXY_URL + '?q=' + encodeURIComponent(seq.join(',')) + '&limit=3',
      { method: 'GET', mode: 'cors' }
    );
    const payload = await resp.json();

    if (token !== oeisLookupToken) return { ok: false, stale: true };
    if (!resp.ok || !payload.ok) {
      throw new Error((payload && payload.error) || ('HTTP ' + resp.status));
    }

    slot.innerHTML = renderOeisMatches(payload);
    return { ok: true, matches: (payload.results || []).length, total: payload.count || 0 };
  } catch (e) {
    if (token !== oeisLookupToken) return { ok: false, stale: true };
    slot.innerHTML = '<div class="oeis-none">Automatic lookup is unavailable right now. Use the OEIS link above.</div>';
    return { ok: false, error: e };
  }
}

async function runOeisInlineQuery(btn) {
  const slotId = btn.dataset.oeisSlot || '';
  const seq = parseSeqData(btn.dataset.oeisSeq);
  const label = btn.dataset.oeisLabel || btn.textContent || 'OEIS query';
  const labelHtml = btn.dataset.oeisLabelHtml || esc(label);
  const slot = document.getElementById(slotId);
  if (!slot || !seq.length) return;

  const token = String(++oeisInlineRequestToken);
  slot.dataset.inlineToken = token;
  slot.innerHTML = '<div class="oeis-none">Looking up OEIS for ' + labelHtml + '...</div>';

  try {
    const resp = await fetch(
      OEIS_PROXY_URL + '?q=' + encodeURIComponent(seq.join(',')) + '&limit=3',
      { method: 'GET', mode: 'cors' }
    );
    const payload = await resp.json();

    if (slot.dataset.inlineToken !== token) return;
    if (!resp.ok || !payload.ok) {
      throw new Error((payload && payload.error) || ('HTTP ' + resp.status));
    }

    slot.innerHTML = renderOeisMatchesWithContext(seq, label, labelHtml, payload);
  } catch (e) {
    if (slot.dataset.inlineToken !== token) return;
    slot.innerHTML = '<div class="oeis-none">Automatic lookup is unavailable right now. <a href="' +
      oeisSearchHref(seq) + '" target="_blank" rel="noopener">Open this query on OEIS</a>.</div>';
  }
}

function renderOeisSection(title, seq, options) {
  if (!seq.length) return '';

  const previewLimit = options.previewLimit || 30;
  const actions = [];
  const trimmed = stripLeadingZeros(seq);
  const diffs = firstDifferences(seq);
  const slotId = options.matchId || options.slotId || '';
  const primaryButtonLabel = options.primaryButtonLabel || 'Load matches';
  const primaryQueryLabel = options.primaryQueryLabel || title;
  const primaryQueryLabelHtml = options.primaryQueryLabelHtml || '';
  const hasInlineSlot = Boolean(slotId);

  if (hasInlineSlot && !options.autoFetch) {
    actions.push(queryButton(primaryButtonLabel, seq, slotId, primaryQueryLabel, primaryQueryLabelHtml));
  } else if (!hasInlineSlot) {
    actions.push(queryLink('Search OEIS', seq));
  }
  if (trimmed.length && trimmed.length < seq.length) {
    actions.push(hasInlineSlot
      ? queryButton(
          'Drop leading 0s',
          trimmed,
          slotId,
          primaryQueryLabel + ' with leading zeros removed',
          (primaryQueryLabelHtml || esc(primaryQueryLabel)) + ' with leading zeros removed'
        )
      : queryLink('Drop leading 0s', trimmed));
  }
  if (diffs.length >= 6) {
    actions.push(hasInlineSlot
      ? queryButton(
          'First differences',
          diffs,
          slotId,
          'first differences of ' + primaryQueryLabel,
          'first differences of ' + (primaryQueryLabelHtml || esc(primaryQueryLabel))
        )
      : queryLink('First differences', diffs));
  }
  if (Array.isArray(options.extraQueries)) {
    for (const query of options.extraQueries) {
      if (!query || !query.seq || !query.seq.length) continue;
      actions.push(hasInlineSlot
        ? queryButton(
            query.label,
            query.seq,
            slotId,
            query.queryLabel || (primaryQueryLabel + ': ' + query.label),
            query.queryLabelHtml || ''
          )
        : queryLink(query.label, query.seq));
    }
  }
  actions.push(queryLink('Open on OEIS', seq));
  actions.push('<button class="run-small" data-copy="' + escAttr(seq.join(', ')) + '" onclick="copyData(this)">Copy terms</button>');

  let html = '<div class="oeis-section oeis-card">';
  html += '<h3>' + (options.titleHtml || esc(title)) + '</h3>';
  if (options.noteHtml || options.note) {
    html += '<div class="oeis-note">' + (options.noteHtml || esc(options.note)) + '</div>';
  }
  html += '<div class="oeis-seq">' + esc(seq.slice(0, previewLimit).join(', ')) +
    (seq.length > previewLimit ? ', \u2026' : '') + '</div>';
  html += '<div class="oeis-actions">' + actions.join('') + '</div>';
  if (slotId) {
    const initial = options.initialResultsHtml ||
      (options.autoFetch
        ? '<div class="oeis-none">Looking up OEIS...</div>'
        : '<div class="oeis-none">Use a query button above to load OEIS matches here.</div>');
    html += '<div id="' + escAttr(slotId) + '" class="oeis-results">' + initial + '</div>';
  }
  html += '</div>';
  return html;
}

function renderSheetOeisGroup(title, rows, options) {
  if (!rows.length) return '';

  const triangle = flattenRows(rows);
  const trimmedTriangle = rows.flatMap(stripLeadingZeros);
  const rowSums = rowSumsOfRows(rows);

  if (rowSums) {
    const extraQueries = [{
      label: 'Triangle',
      seq: triangle,
      queryLabel: title + ': triangle read by rows',
      queryLabelHtml: (options.titleHtml || esc(title)) + ': triangle read by rows',
    }];
    if (trimmedTriangle.join(',') !== triangle.join(',')) {
      extraQueries.push({
        label: 'Triangle, strip 0s',
        seq: trimmedTriangle,
        queryLabel: title + ': triangle after stripping initial zeros in each row',
        queryLabelHtml: (options.titleHtml || esc(title)) + ': triangle after stripping initial zeros in each row',
      });
    }
    return renderOeisSection(title, rowSums, {
      titleHtml: options.titleHtml,
      noteHtml: options.noteHtml,
      note: options.note || ('Row sums of the displayed ' + title + '. Use the triangle actions for the full row-by-row data.'),
      previewLimit: options.previewLimit || 24,
      slotId: options.slotId,
      primaryButtonLabel: options.primaryButtonLabel || 'Row sums',
      primaryQueryLabel: title + ': row sums',
      primaryQueryLabelHtml: (options.titleHtml || esc(title)) + ': row sums',
      extraQueries,
    });
  }

  return renderOeisSection(title, triangle, {
    titleHtml: options.titleHtml,
    noteHtml: options.fallbackNoteHtml,
    note: options.fallbackNote || ('Triangle by rows for the displayed ' + title + '.'),
    previewLimit: options.previewLimit || 42,
    slotId: options.slotId,
    primaryButtonLabel: options.primaryButtonLabel || 'Triangle',
    primaryQueryLabel: title + ': triangle read by rows',
    primaryQueryLabelHtml: (options.titleHtml || esc(title)) + ': triangle read by rows',
  });
}

function renderSheetOeisHelpers(view, data) {
  const clean = data.filter(item => !item.error);
  if (!clean.length) return '';

  const groups = [];
  const slotBase = 'sheet-oeis-' + (++oeisInlineSlotSerial);
  if (view === 'core') {
    groups.push(renderSheetOeisGroup('p(t)', clean.map(item => item.coefficients), {
      titleHtml: katexInline('p(t)'),
      noteHtml: 'Row sums of the input polynomials ' + katexInline('p(t)') + '. Triangle actions search the full coefficient triangle.',
      slotId: slotBase + '-0',
    }));
    groups.push(renderSheetOeisGroup('a(t)', clean.map(item => item.a.coefficients), {
      titleHtml: katexInline('a(t)'),
      noteHtml: 'Row sums of the symmetric parts ' + katexInline('a(t)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-1',
    }));
    groups.push(renderSheetOeisGroup('gamma(a)', clean.map(item => item.gamma_a_coefficients || ['0']), {
      titleHtml: katexInline('\\gamma_a(z)'),
      noteHtml: 'Row sums of the gamma-coefficient rows for ' + katexInline('\\gamma_a(z)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-2',
    }));
    groups.push(renderSheetOeisGroup('b(t)', clean.map(item => item.b.coefficients), {
      titleHtml: katexInline('b(t)'),
      noteHtml: 'Row sums of the companion parts ' + katexInline('b(t)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-3',
    }));
    groups.push(renderSheetOeisGroup('gamma(b)', clean.map(item => item.gamma_b_coefficients || ['0']), {
      titleHtml: katexInline('\\gamma_b(z)'),
      noteHtml: 'Row sums of the gamma-coefficient rows for ' + katexInline('\\gamma_b(z)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-4',
    }));
  } else if (view === 'companion') {
    groups.push(renderSheetOeisGroup('f(t)', clean.map(item => item.f_polynomial.coefficients), {
      titleHtml: katexInline('f(t)'),
      noteHtml: 'Row sums of the auxiliary polynomials ' + katexInline('f(t)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-0',
    }));
    groups.push(renderSheetOeisGroup('R_d(f)', clean.map(item => item.r_transform_of_f.coefficients), {
      titleHtml: katexInline('R_d(f)'),
      noteHtml: 'Row sums of ' + katexInline('R_d(f)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-1',
    }));
    groups.push(renderSheetOeisGroup('a~(t)', clean.map(item => item.r_a.coefficients), {
      titleHtml: katexInline('\\widetilde{a}(t)'),
      noteHtml: 'Row sums of the companion symmetric parts ' + katexInline('\\widetilde{a}(t)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-2',
    }));
    groups.push(renderSheetOeisGroup('b~(t)', clean.map(item => item.r_b.coefficients), {
      titleHtml: katexInline('\\widetilde{b}(t)'),
      noteHtml: 'Row sums of the companion parts ' + katexInline('\\widetilde{b}(t)') + '. Triangle actions search the displayed triangle.',
      slotId: slotBase + '-3',
    }));
  } else {
    groups.push(renderSheetOeisGroup('magic coords', clean.map(item => item.magic.coordinates), {
      titleHtml: katexInline('c_j'),
      noteHtml: 'Row sums of the displayed coordinate rows when they are integral. Triangle actions search the full coordinate triangle for ' + katexInline('c_j') + '.',
      fallbackNoteHtml: 'Triangle by rows for the basis coordinates ' + katexInline('c_j') + '.',
      slotId: slotBase + '-0',
    }));
  }

  return '<div class="oeis-note">OEIS helpers for the displayed data above. These use the same full-width card layout as the main OEIS tab; the main sequence is the row-sum sequence when available, and the triangle actions search the full displayed object row by row.</div>' +
    '<div class="sheet-oeis-stack">' + groups.join('') + '</div>';
}

async function runOeis() {
  if (!wasmReady) { setStatus('WASM not loaded yet.'); return; }
  let sections = [];
  if (activeWorkflow === 'integer-sequence') {
    const terms = getIntegerTerms();
    if (!terms.length) { setStatus('No input.'); return; }
    sections = [{
      title: 'Integer sequence',
      seq: terms,
      note: 'The terms exactly as entered. First differences and leading-zero variants are available below.',
      previewLimit: 42,
    }];
  } else {
    const text = getInput();
    if (!text.trim()) { setStatus('No input.'); return; }
    const data = JSON.parse(check_properties(text));
    const polys = data.filter(d => !d.error && d.coefficients);
    if (polys.length === 0) { setResult('<span class="err">No valid polynomials found.</span>'); return; }

    if (activeWorkflow === 'polynomial') {
      if (polys.length !== 1) {
        setResult('<span class="err">This view expects exactly one polynomial.</span>');
        return;
      }
      sections = [{
        title: 'Coefficient sequence',
        seq: polys[0].coefficients,
        note: 'Coefficients in ascending degree, from the constant term onward.',
        previewLimit: 42,
      }];
    } else {
      const rowSums = polys.map(d => sumIntegerRow(d.coefficients, false));
      const triangle = polys.flatMap(d => d.coefficients);
      const altRowSums = polys.map(d => sumIntegerRow(d.coefficients, true));
      sections = [
        {
          title: 'Row sums (evaluate at t = 1)',
          titleHtml: 'Row sums (evaluate at ' + katexInline('t=1') + ')',
          seq: rowSums,
          note: 'Often the fastest way to identify a polynomial family.',
          previewLimit: 24,
        },
        {
          title: 'Alternating row sums (evaluate at t = -1)',
          titleHtml: 'Alternating row sums (evaluate at ' + katexInline('t=-1') + ')',
          seq: altRowSums,
          note: 'Useful when the family has an especially clean signed specialization.',
          previewLimit: 24,
        },
        {
          title: 'Triangle read by rows',
          seq: triangle,
          note: 'Flatten the coefficient triangle row by row.',
          previewLimit: 42,
        },
      ];
    }
  }

  let html = '';
  const token = ++oeisLookupToken;

  html += '<div class="oeis-note">Automatic matches below are fetched through the OEIS proxy at <code>arxiv.symmetricfunctions.com</code>. The buttons still open OEIS directly with the full query.</div>';
  html += sections.map((section, idx) => renderOeisSection(section.title, section.seq, {
    titleHtml: section.titleHtml,
    note: section.note,
    previewLimit: section.previewLimit,
    matchId: 'oeis-match-' + token + '-' + idx,
    autoFetch: true,
  })).join('');

  setResult(html);
  setStatus('Looking up OEIS matches for ' + sections.length + ' sequence' + (sections.length === 1 ? '' : 's') + '...');

  const lookups = await Promise.all(sections.map((section, idx) =>
    fetchOeisMatches(section.seq, 'oeis-match-' + token + '-' + idx, token)
  ));

  if (token !== oeisLookupToken) return;

  const withMatches = lookups.filter(item => item && item.ok && item.matches > 0).length;
  const failed = lookups.filter(item => item && item.error).length;
  if (failed) {
    setStatus('Fetched OEIS matches for ' + withMatches + ' of ' + sections.length + ' sequence' + (sections.length === 1 ? '' : 's') + '; some proxy lookups failed, but direct OEIS links remain available.');
  } else {
    setStatus('Fetched OEIS matches for ' + withMatches + ' of ' + sections.length + ' sequence' + (sections.length === 1 ? '' : 's') + '.');
  }
}
