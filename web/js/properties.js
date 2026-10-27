function abbreviatePoly(full) {
  // Parse terms from formatted polynomial like "1 + 11t + 11t^2 + t^3"
  // Show first 3 terms + last term with dots if more than 5 terms
  const terms = full.match(/[+-]?[^+-]+/g);
  if (!terms) return full;
  const clean = terms.map(t => t.trim()).filter(t => t);
  if (clean.length <= 5) return full;

  const joinSignedTerms = items => items.map((raw, idx) => {
    let term = raw.trim();
    let sign = '+';
    if (term.startsWith('+') || term.startsWith('-')) {
      sign = term[0];
      term = term.slice(1).trim();
    }
    if (idx === 0) return sign === '-' ? '- ' + term : term;
    return (sign === '-' ? ' - ' : ' + ') + term;
  }).join('');

  return joinSignedTerms([...clean.slice(0, 3), '\u2026', clean[clean.length - 1]]);
}

function polyTooltip(full, abbrev) {
  if (full === abbrev) return esc(full);
  return '<span class="poly-tip" title="' + escAttr(full) + '">' + esc(abbrev) + '</span>';
}

function concavityLabel(item) {
  if (item.real_rooted) {
    return mini('RR', true, 'Real-rooted => ultra-log-concave => log-concave => unimodal');
  }
  if (item.ultra_log_concave) {
    return mini('ULC', true, 'Ultra-log-concave => log-concave => unimodal (but not real-rooted)');
  }
  if (item.log_concave) {
    return mini('LC', true, 'Log-concave => unimodal (not ultra-log-concave)');
  }
  return mini('—', null, '', '', 'Not log-concave');
}

function symmetryLabel(item) {
  if (item.gamma_positive && item.gamma_coefficients) {
    return mini('\u03b3+', true, '\u03b3-positive (implies palindromic). \u03b3 = [' + item.gamma_coefficients.join(', ') + ']');
  }
  if (item.palindromic) {
    return mini('pal', true, 'Palindromic (but not \u03b3-positive)');
  }
  return mini('—', null, '', '', 'Not palindromic');
}

let lastPropsData = null;  // stored for interlacing check
let lastDecompData = null;
let activeDecompView = 'core';

function gammaBasisLatex(item) {
  if (!item.gamma_coefficients) return null;
  let initialZeros = 0;
  while (initialZeros < item.coefficients.length && isZeroValue(item.coefficients[initialZeros])) {
    initialZeros++;
  }
  const effectiveDegree = Math.max(0, item.degree - initialZeros);
  const terms = [];
  for (let i = 0; i < item.gamma_coefficients.length; i++) {
    const coefficient = normalizeIntegerString(item.gamma_coefficients[i]);
    if (coefficient === null || coefficient === '0') continue;
    const negative = coefficient.startsWith('-');
    const abs = negative ? coefficient.slice(1) : coefficient;
    const exponent = effectiveDegree - 2 * i;
    const tExponent = initialZeros + i;
    const basis = (tExponent ? (tExponent === 1 ? 't' : 't^{' + tExponent + '}') : '') +
      (exponent ? '(1+t)^{' + exponent + '}' : '');
    const magnitude = basis && abs === '1' ? basis : abs + basis;
    terms.push({ negative, magnitude });
  }
  if (!terms.length) return '0';
  return terms.map((term, index) => {
    if (index === 0) return (term.negative ? '-' : '') + term.magnitude;
    return (term.negative ? ' - ' : ' + ') + term.magnitude;
  }).join('');
}

function renderSinglePolynomial(data) {
  if (data.length !== 1) {
    setResult('<span class="err">This view expects exactly one polynomial. Choose “Polynomial sequence” to analyze several rows.</span>');
    return;
  }
  const item = data[0];
  if (item.error) {
    setResult('<span class="err">' + esc(item.error) + '</span>');
    return;
  }
  const rootClass = item.real_rooted ? 'yes' : 'no';
  const repeated = item.squarefree_degree < item.degree;
  const rootDetail = item.degree === 0
    ? 'A nonzero constant has no zeros.'
    : item.distinct_real_roots + ' distinct real zero' + (item.distinct_real_roots === 1 ? '' : 's') +
      ': ' + item.negative_real_roots + ' negative, ' + (item.zero_is_root ? 'one at 0' : 'none at 0') +
      ', and ' + item.positive_real_roots + ' positive.' + (repeated ? ' Repeated zeros are present.' : '');
  let gammaValue = 'Not available';
  let gammaDetail = 'A gamma expansion requires palindromic coefficients.';
  if (item.gamma_coefficients) {
    gammaValue = item.gamma_positive ? 'Gamma-positive' : 'Gamma expansion found';
    gammaDetail = 'γ = [' + item.gamma_coefficients.join(', ') + ']';
    const expansion = gammaBasisLatex(item);
    if (expansion !== null) gammaDetail += '<div class="recurrence-display">' + katexInline('p(t)=' + expansion) + '</div>';
  }
  const shape = [
    item.log_concave ? 'log-concave' : 'not log-concave',
    item.ultra_log_concave ? 'ultra-log-concave' : 'not ultra-log-concave',
    item.palindromic ? 'palindromic' : 'not palindromic',
  ].join('; ');
  const html = '<div class="single-summary">' +
    '<div class="single-poly"><div class="single-poly-expression">' + katexInline('p(t)=' + item.polynomial.replace(/\^([0-9]+)/g, '^{$1}')) + '</div>' +
    '<div class="summary-detail">Degree ' + item.degree + '; coefficients [' + esc(item.coefficients.join(', ')) + ']</div></div>' +
    '<div class="summary-grid">' +
      '<div class="summary-card"><h3>Real zeros</h3><div class="summary-value ' + rootClass + '">' +
        (item.real_rooted ? 'All zeros are real' : 'Not all zeros are real') + '</div><div class="summary-detail">' + esc(rootDetail) + '</div></div>' +
      '<div class="summary-card"><h3>Gamma expansion</h3><div class="summary-value ' + (item.gamma_positive ? 'yes' : '') + '">' +
        esc(gammaValue) + '</div><div class="summary-detail">' + gammaDetail + '</div></div>' +
      '<div class="summary-card"><h3>Coefficient shape</h3><div class="summary-value">' +
        (item.palindromic ? 'Symmetric' : 'Asymmetric') + '</div><div class="summary-detail">' + esc(shape) +
        '</div></div>' +
    '</div>' +
    '<div class="oeis-actions"><button class="run-small" onclick="openOeisTab()">Look up these coefficients on OEIS</button></div>' +
    '</div>';
  setResult(html);
}

function openOeisTab() {
  switchTab('oeis');
  run('oeis');
}

function renderProperties(text) {
  const t0 = performance.now();
  const data = JSON.parse(check_properties(text));
  lastPropsData = data;
  const ms = (performance.now() - t0).toFixed(1);
  if (activeWorkflow === 'polynomial') {
    renderSinglePolynomial(data);
    setStatus('Polynomial checked exactly in ' + ms + ' ms.');
  } else {
    const counts = computeInterlacingCounts(data);
    renderPropsTable(data, counts);
    setStatus(data.length + ' polynomial(s) checked in ' + ms + ' ms, including consecutive interlacing profile.');
  }
}

function computeInterlacingCounts(data) {
  const n = data.length;
  const counts = new Array(n).fill(null);
  if (n === 0) return counts;
  counts[0] = 0;

  for (let i = 1; i < n; i++) {
    const itemI = data[i];
    if (itemI.error || !itemI.coefficients) { counts[i] = null; continue; }
    let count = 0;
    const ci = JSON.stringify(itemI.coefficients);

    for (let j = i - 1; j >= 0; j--) {
      const itemJ = data[j];
      if (itemJ.error || !itemJ.coefficients) break;
      const cj = JSON.stringify(itemJ.coefficients);
      const result = JSON.parse(check_interlacing_pair(cj, ci));
      if (result.strict || result.weak) {
        count++;
      } else {
        break;
      }
    }
    counts[i] = count;
  }

  return counts;
}

function renderPropsTable(data, interlacingCounts) {
  let rows = [];
  let idx = 1;
  for (const item of data) {
    if (item.error) {
      rows.push('<tr><td class="idx">' + idx + '</td><td colspan="4" class="err">' +
                esc(item.error) + '</td></tr>');
      idx++;
      continue;
    }
    const full = item.polynomial;
    const abbrev = abbreviatePoly(full);
    let intCell;
    if (interlacingCounts === null) {
      intCell = '<td class="icons interlace-pending">' + mini('—', null, '', '', 'Run interlacing to populate this column') + '</td>';
    } else {
      const count = interlacingCounts[idx - 1];
      if (count === undefined || count === null) {
        intCell = '<td class="icons">' + mini('—', null, '', '', 'Interlacing data unavailable') + '</td>';
      } else if (count > 0) {
        intCell = '<td class="icons">' + mini(String(count), true, 'Interlaces ' + count + ' consecutive previous polynomial(s) before first failure') + '</td>';
      } else {
        intCell = '<td class="icons">' + mini('0', null, '', '', 'No previous consecutive interlacing before first failure') + '</td>';
      }
    }
    rows.push(
      '<tr>' +
      '<td class="idx">' + idx + '</td>' +
      '<td class="poly">' + polyTooltip(full, abbrev) + '</td>' +
      '<td class="icons">' + concavityLabel(item) + '</td>' +
      '<td class="icons">' + symmetryLabel(item) + '</td>' +
      intCell +
      '</tr>'
    );
    idx++;
  }
  const intHeader = '<th>' + headerLabel('interlaces', 'Number of earlier polynomials with |\u0394deg| \u2264 1 that this one interlaces (B\u00e9zout method).') + '</th>';
  const header = '<tr class="hdr"><th></th><th>polynomial</th>' +
    '<th>' + headerLabel('concavity', 'Concavity hierarchy: RealRooted => ULC => LogConcave.') + '</th>' +
    '<th>' + headerLabel('symmetry', '\u03b3-positive => palindromic.') + '</th>' +
    intHeader + '</tr>';
  setResult('<div class="prop-table-wrap"><div class="prop-table-scroll"><table class="prop-table">' + header + rows.join('') + '</table></div></div>');
}

function serializeCoeffRow(coeffs) {
  return coeffs && coeffs.length ? coeffs.join(', ') : '0';
}

function propsForRows(rows) {
  const lines = rows.map(serializeCoeffRow).join('\n');
  return JSON.parse(check_properties(lines));
}

function normalizeIntegerString(value) {
  const match = String(value).trim().match(/^([+-]?)(\d+)$/);
  if (!match) return null;
  const digits = match[2].replace(/^0+(?=\d)/, '');
  return digits === '0' ? '0' : (match[1] === '-' ? '-' : '') + digits;
}

function coeffsToPolyText(coeffs, variable) {
  if (!coeffs || !coeffs.length) return '0';
  const terms = [];
  for (let i = 0; i < coeffs.length; i++) {
    const coefficient = normalizeIntegerString(coeffs[i]);
    if (coefficient === null || coefficient === '0') continue;
    const negative = coefficient.startsWith('-');
    const abs = negative ? coefficient.slice(1) : coefficient;
    let term = '';
    if (i === 0) {
      term = abs;
    } else {
      const coeffPart = abs === '1' ? '' : abs;
      term = coeffPart + variable + (i === 1 ? '' : '^' + i);
    }
    terms.push((negative ? '- ' : '+ ') + term);
  }
  if (!terms.length) return '0';
  return terms.map((term, idx) => idx === 0 ? term.replace(/^\+\s/, '') : term).join(' ');
}

function gammaPolyText(coeffs) {
  if (!coeffs) return null;
  return coeffsToPolyText(coeffs, 'z');
}

function infoTip(text) {
  return ' <span class="info">i<span class="tip">' + esc(text) + '</span></span>';
}

function infoTipHtml(html) {
  return ' <span class="info">i<span class="tip">' + html + '</span></span>';
}

function headerLabel(label, tip) {
  return esc(label) + (tip ? infoTip(tip) : '');
}

function mathHeader(tex, tipHtml) {
  return katexInline(tex) + (tipHtml ? infoTipHtml(tipHtml) : '');
}

let floatingInfoTipEl = null;
let activeInfoAnchor = null;

function ensureFloatingInfoTip() {
  if (floatingInfoTipEl) return floatingInfoTipEl;
  floatingInfoTipEl = document.createElement('div');
  floatingInfoTipEl.className = 'floating-info-tip';
  floatingInfoTipEl.setAttribute('aria-hidden', 'true');
  document.body.appendChild(floatingInfoTipEl);
  return floatingInfoTipEl;
}

function positionFloatingInfoTip(anchor) {
  if (!anchor || !floatingInfoTipEl) return;

  const rect = anchor.getBoundingClientRect();
  const margin = 8;

  floatingInfoTipEl.style.visibility = 'hidden';
  floatingInfoTipEl.style.display = 'block';
  floatingInfoTipEl.style.left = '0px';
  floatingInfoTipEl.style.top = '0px';

  const tipRect = floatingInfoTipEl.getBoundingClientRect();
  let left = rect.left + rect.width / 2 - tipRect.width / 2;
  left = Math.max(margin, Math.min(left, window.innerWidth - tipRect.width - margin));

  let top = rect.top - tipRect.height - 10;
  let placement = 'top';
  if (top < margin) {
    top = rect.bottom + 10;
    placement = 'bottom';
  }

  floatingInfoTipEl.dataset.placement = placement;
  floatingInfoTipEl.style.left = left + 'px';
  floatingInfoTipEl.style.top = top + 'px';
  floatingInfoTipEl.style.visibility = 'visible';
}

function showFloatingInfoTip(anchor) {
  const tip = anchor && anchor.querySelector('.tip');
  if (!tip) return;
  activeInfoAnchor = anchor;
  const layer = ensureFloatingInfoTip();
  layer.innerHTML = tip.innerHTML;
  renderInlineMath(layer);
  positionFloatingInfoTip(anchor);
}

function hideFloatingInfoTip(anchor) {
  if (anchor && activeInfoAnchor && anchor !== activeInfoAnchor) return;
  activeInfoAnchor = null;
  if (!floatingInfoTipEl) return;
  floatingInfoTipEl.style.display = 'none';
  floatingInfoTipEl.innerHTML = '';
}

function installInfoTooltipLayer() {
  document.addEventListener('mouseover', event => {
    const info = event.target.closest('.info');
    if (!info || info.contains(event.relatedTarget)) return;
    showFloatingInfoTip(info);
  }, true);

  document.addEventListener('mouseout', event => {
    const info = event.target.closest('.info');
    if (!info || info.contains(event.relatedTarget)) return;
    hideFloatingInfoTip(info);
  }, true);

  window.addEventListener('scroll', () => {
    if (activeInfoAnchor) positionFloatingInfoTip(activeInfoAnchor);
  }, true);

  window.addEventListener('resize', () => {
    if (activeInfoAnchor) positionFloatingInfoTip(activeInfoAnchor);
  });

  document.addEventListener('DOMContentLoaded', () => renderInlineMath(document));
}

function renderInlineMath(root) {
  if (!root || !window.katex || typeof window.katex.render !== 'function') return;
  root.querySelectorAll('[data-tex]').forEach(node => {
    if (node.dataset.katexRendered === 'true') return;
    window.katex.render(node.dataset.tex, node, {
      displayMode: false,
      throwOnError: false,
    });
    node.dataset.katexRendered = 'true';
  });
}

function polyTextToLatex(text) {
  return String(text)
    .replace(/…/g, '\\cdots')
    .replace(/([A-Za-z])\^(\d+)/g, '$1^{$2}');
}

function allZeroSeq(values) {
  return !values || !values.length || values.every(isZeroValue);
}

function polyTextInner(text) {
  if (!text) return '\u2014';
  const abbrev = abbreviatePoly(text);
  const rendered = katexInline(polyTextToLatex(abbrev));
  if (text === abbrev) return rendered;
  return '<span class="poly-tip" title="' + escAttr(text) + '">' + rendered + '</span>';
}

function polyTextCell(text, extraClass, compact) {
  const cls = 'poly' + (compact ? '' : ' wide') + (extraClass ? ' ' + extraClass : '');
  return '<td class="' + cls + '">' + polyTextInner(text) + '</td>';
}

function stackedPolyCell(mainLabelTex, mainText, subLabelTex, subText, compact) {
  const cls = 'poly' + (compact ? '' : ' wide doublewide');
  return '<td class="' + cls + '"><div class="decomp-stack">' +
    '<div class="decomp-stack-line">' +
    '<div class="decomp-stack-label">' + katexInline(mainLabelTex) + '</div>' +
    '<div class="decomp-stack-value">' + polyTextInner(mainText) + '</div>' +
    '</div>' +
    '<div class="decomp-stack-line sub">' +
    '<div class="decomp-stack-label">' + katexInline(subLabelTex) + '</div>' +
    '<div class="decomp-stack-value">' + polyTextInner(subText) + '</div>' +
    '</div>' +
    '</div></td>';
}

function vectorInner(values, compact) {
  if (values == null) {
    return '<div class="coords unavailable' + (compact ? ' compact' : '') + '">\u2014</div>';
  }
  const text = values.length ? values.join(', ') : '0';
  return '<div class="coords' + (compact ? ' compact' : '') + '" title="' + escAttr(text) + '">' + esc(text) + '</div>';
}

function vectorCell(values, compact) {
  return '<td>' + vectorInner(values, compact) + '</td>';
}

function stackedVectorCell(mainLabelTex, mainValues, subLabelTex, subValues, compact) {
  return '<td><div class="decomp-stack">' +
    '<div class="decomp-stack-line">' +
    '<div class="decomp-stack-label">' + katexInline(mainLabelTex) + '</div>' +
    '<div class="decomp-stack-value">' + vectorInner(mainValues, compact) + '</div>' +
    '</div>' +
    '<div class="decomp-stack-line sub">' +
    '<div class="decomp-stack-label">' + katexInline(subLabelTex) + '</div>' +
    '<div class="decomp-stack-value">' + vectorInner(subValues, compact) + '</div>' +
    '</div>' +
    '</div></td>';
}

function partialSumTooltip(magic) {
  if (!magic || !Array.isArray(magic.partial_sum_checks) || !magic.partial_sum_checks.length) {
    return 'partial-sum data unavailable';
  }
  const failAt = magic.partial_sum_checks.findIndex(flag => !flag);
  if (failAt === -1) {
    return 'The inequalities L_j <= R_j hold for all checked indices j = 0..' + (magic.partial_sum_checks.length - 1) + '.';
  }
  const left = Array.isArray(magic.left_partial_sums) ? magic.left_partial_sums[failAt] : '?';
  const right = Array.isArray(magic.right_partial_sums) ? magic.right_partial_sums[failAt] : '?';
  return 'Fails first at j = ' + failAt + ': L_' + failAt + ' = ' + left + ' > R_' + failAt + ' = ' + right + '.';
}

function mini(code, flag, yesTitle, noTitle, neutralTitle) {
  let cls = 'neutral';
  let title = neutralTitle || (code + ': not available');
  if (flag === true) {
    cls = 'yes';
    title = yesTitle;
  } else if (flag === false) {
    cls = 'no';
    title = noTitle;
  }
  return '<span class="mini ' + cls + '">' + esc(code) + '<span class="tip">' + esc(title) + '</span></span>';
}

function iconRow(items) {
  return '<div class="decomp-icons">' + items.join('') + '</div>';
}

function mergeDecompositionData(decompRows, inputText) {
  const inputProps = JSON.parse(check_properties(inputText));
  const aProps = propsForRows(decompRows.map(item => item.error ? ['0'] : (item.a.coefficients || ['0'])));
  const bProps = propsForRows(decompRows.map(item => item.error ? ['0'] : (item.b.coefficients || ['0'])));

  return decompRows.map((item, idx) => {
    if (item.error) return item;
    const inputProp = inputProps[idx] || {};
    const aProp = aProps[idx] || {};
    const bProp = bProps[idx] || {};
    return {
      ...item,
      gamma_input_coefficients: inputProp.gamma_coefficients || null,
      gamma_input_text: gammaPolyText(inputProp.gamma_coefficients || null),
      gamma_a_coefficients: aProp.gamma_coefficients || null,
      gamma_a_text: gammaPolyText(aProp.gamma_coefficients || null),
      gamma_a_positive: aProp.gamma_positive,
      gamma_b_coefficients: bProp.gamma_coefficients || null,
      gamma_b_text: gammaPolyText(bProp.gamma_coefficients || null),
      gamma_b_positive: bProp.gamma_positive,
    };
  });
}

function renderDecompositions(text) {
  const t0 = performance.now();
  lastDecompData = mergeDecompositionData(JSON.parse(analyze_decompositions(text)), text);
  const ms = (performance.now() - t0).toFixed(1);
  renderDecompositionView();
  setStatus(lastDecompData.length + ' polynomial(s) analyzed in ' + ms + ' ms.');
}

function setDecompView(view) {
  activeDecompView = view;
  if (lastDecompData) renderDecompositionView();
}

function renderDecompositionView() {
  if (!lastDecompData) {
    setResult('<div class="decomp-empty">Run the decomposition analysis to populate the sheet.</div>');
    return;
  }

  const subtabs = [
    ['core', 'Symmetric Split'],
    ['companion', 'R_d Split'],
    ['magic', 'Magic basis'],
  ].map(([key, label]) =>
    '<button class="decomp-subtab' + (activeDecompView === key ? ' active' : '') +
    '" onclick="setDecompView(\'' + key + '\')">' + esc(label) + '</button>'
  ).join('');

  const helpHtml = {
    core: '<strong>Symmetric decomposition</strong>For a polynomial ' + katexInline('p(t)') + ' of degree at most ' + katexInline('d') + ', define ' + katexInline('I_d(p)(t)=t^d p(1/t)') + '. Then there are unique polynomials ' + katexInline('a(t)') + ' and ' + katexInline('b(t)') + ' such that ' + katexInline('p(t)=a(t)+tb(t)') + ', with ' + katexInline('I_d(a)=a') + ' and ' + katexInline('I_{d-1}(b)=b') + '; equivalently ' + katexInline('I_d(p)=a+b') + '. This table lists the coefficient vectors of ' + katexInline('p(t)') + ', ' + katexInline('a(t)') + ', and ' + katexInline('b(t)') + ', with the corresponding gamma-coefficient vectors directly underneath. Reference: Petter Br&auml;nd&eacute;n and Liam Solus, <a href="https://academic.oup.com/imrn/article/2021/10/7764/5425553" target="_blank" rel="noopener">Symmetric decompositions and real-rootedness</a>, <i>International Mathematics Research Notices</i> 2021:10, 7764-7810.',
    companion: '<strong>Companion ' + katexInline('R_d') + ' decomposition</strong>Starting from the input polynomial ' + katexInline('p(t)=\\sum_{k=0}^d p_k t^k') + ', define the auxiliary polynomial ' + katexInline('f(t)=(1+t)^d p\\!\\left(\\frac{t}{1+t}\\right)=\\sum_{k=0}^d p_k t^k(1+t)^{d-k}') + '. Then define ' + katexInline('R_d(f)(t)=(-1)^d f(-1-t)') + '. This tab shows the decomposition ' + katexInline('f(t)=\\widetilde{a}(t)+t\\widetilde{b}(t)') + ', written as coefficient vectors, where ' + katexInline('R_d(\\widetilde{a})=\\widetilde{a}') + ' and ' + katexInline('R_{d-1}(\\widetilde{b})=\\widetilde{b}') + '.',
    magic: '<strong>Magic-basis coordinates</strong>This tab records the coordinates ' + katexInline('c_j') + ' in the basis ' + katexInline('t^j(1+t)^{d-j}') + '. The flag ' + katexInline('S\\le') + ' checks the partial-sum inequalities ' + katexInline('L_j=c_0+\\cdots+c_j\\le R_j=c_d+\\cdots+c_{d-j}') + '; hover it to see the first failing index when one occurs.',
  }[activeDecompView];
  const clean = lastDecompData.filter(item => !item.error);
  let header = '';
  let rows = [];
  let colspan = 7;

  if (activeDecompView === 'core') {
    colspan = 5;
    const compactP = clean.every(item => allZeroSeq(item.coefficients));
    const compactA = clean.every(item => allZeroSeq(item.a.coefficients));
    const compactGammaA = clean.every(item => allZeroSeq(item.gamma_a_coefficients));
    const compactB = clean.every(item => allZeroSeq(item.b.coefficients));
    const compactGammaB = clean.every(item => allZeroSeq(item.gamma_b_coefficients));
    const compactAStack = compactA && compactGammaA;
    const compactBStack = compactB && compactGammaB;
    header = '<tr>' +
      '<th></th>' +
      '<th>' + headerLabel('deg', 'Degree after trimming trailing zeros.') + '</th>' +
      '<th>' + mathHeader('p(t)', 'Coefficient vector of the input polynomial ' + katexInline('p(t)') + '.') + '</th>' +
      '<th>' + mathHeader('a(t),\\ \\gamma_a(z)', 'The first line gives the coefficient vector of the unique ' + katexInline('I_d') + '-symmetric part ' + katexInline('a(t)') + ', and the second line gives the gamma-coefficient vector of ' + katexInline('\\gamma_a(z)=\\sum_i \\gamma_i z^i') + '.') + '</th>' +
      '<th>' + mathHeader('b(t),\\ \\gamma_b(z)', 'The first line gives the coefficient vector of the unique ' + katexInline('I_{d-1}') + '-symmetric companion ' + katexInline('b(t)') + ', and the second line gives the gamma-coefficient vector of ' + katexInline('\\gamma_b(z)=\\sum_i \\gamma_i z^i') + '.') + '</th>' +
      '<th>' + headerLabel('Properties', 'aRR: a real-rooted. ga+: gamma(a) has nonnegative coefficients. bRR: b real-rooted. gb+: gamma(b) has nonnegative coefficients. b≺a: b interlaces a. AI: alternatingly increasing.') + '</th>' +
      '</tr>';
    rows = lastDecompData.map((item, idx) => {
      if (item.error) {
        return '<tr><td class="idx">' + (idx + 1) + '</td><td colspan="' + colspan + '" class="err">' + esc(item.error) + '</td></tr>';
      }
      return '<tr>' +
        '<td class="idx">' + (idx + 1) + '</td>' +
        '<td class="compact">' + item.degree + '</td>' +
        vectorCell(item.coefficients, compactP) +
        stackedVectorCell('a(t)', item.a.coefficients, '\\gamma_a(z)', item.gamma_a_coefficients, compactAStack) +
        stackedVectorCell('b(t)', item.b.coefficients, '\\gamma_b(z)', item.gamma_b_coefficients, compactBStack) +
        '<td>' + iconRow([
          mini('aRR', item.a_real_rooted, 'a(t) is real-rooted', 'a(t) is not real-rooted'),
          mini('ga+', item.gamma_a_positive, 'gamma(a) has nonnegative coefficients', 'gamma(a) has a negative coefficient', 'gamma(a): unavailable'),
          mini('bRR', item.b_real_rooted, 'b(t) is real-rooted', 'b(t) is not real-rooted'),
          mini('gb+', item.gamma_b_positive, 'gamma(b) has nonnegative coefficients', 'gamma(b) has a negative coefficient', 'gamma(b): unavailable'),
          mini('b\u226aa', item.b_interlaces_a, 'b interlaces a', 'b does not interlace a', 'b versus a: unavailable'),
          mini('AI', item.alternatingly_increasing, 'alternatingly increasing', 'not alternatingly increasing'),
        ]) + '</td>' +
        '</tr>';
    });
  } else if (activeDecompView === 'companion') {
    colspan = 6;
    const compactF = clean.every(item => allZeroSeq(item.f_polynomial.coefficients));
    const compactR = clean.every(item => allZeroSeq(item.r_transform_of_f.coefficients));
    const compactRA = clean.every(item => allZeroSeq(item.r_a.coefficients));
    const compactRB = clean.every(item => allZeroSeq(item.r_b.coefficients));
    header = '<tr>' +
      '<th></th>' +
      '<th>' + headerLabel('deg', 'Degree after trimming trailing zeros.') + '</th>' +
      '<th>' + mathHeader('f(t)', 'The auxiliary polynomial attached to the input ' + katexInline('p(t)') + ', namely ' + katexInline('f(t)=(1+t)^d p\\!\\left(\\frac{t}{1+t}\\right)=\\sum_{k=0}^d p_k t^k(1+t)^{d-k}') + '.') + '</th>' +
      '<th>' + mathHeader('R_d(f)', 'Coefficient vector of ' + katexInline('R_d(f)(t)=(-1)^d f(-1-t)') + '.') + '</th>' +
      '<th>' + mathHeader('\\widetilde{a}(t)', 'Coefficient vector of the unique ' + katexInline('R_d') + '-symmetric part in ' + katexInline('f(t)=\\widetilde{a}(t)+t\\widetilde{b}(t)') + ', so ' + katexInline('R_d(\\widetilde{a})=\\widetilde{a}') + '.') + '</th>' +
      '<th>' + mathHeader('\\widetilde{b}(t)', 'Coefficient vector of the unique companion part in ' + katexInline('f(t)=\\widetilde{a}(t)+t\\widetilde{b}(t)') + ', so ' + katexInline('R_{d-1}(\\widetilde{b})=\\widetilde{b}') + '.') + '</th>' +
      '<th>' + headerLabel('Properties', 'R≺f: R_d(f) interlaces f. M+: all magic-basis coordinates are nonnegative. S≤: all left partial sums are at most the corresponding right partial sums.') + '</th>' +
      '</tr>';
    rows = lastDecompData.map((item, idx) => {
      if (item.error) {
        return '<tr><td class="idx">' + (idx + 1) + '</td><td colspan="' + colspan + '" class="err">' + esc(item.error) + '</td></tr>';
      }
      return '<tr>' +
        '<td class="idx">' + (idx + 1) + '</td>' +
        '<td class="compact">' + item.degree + '</td>' +
        vectorCell(item.f_polynomial.coefficients, compactF) +
        vectorCell(item.r_transform_of_f.coefficients, compactR) +
        vectorCell(item.r_a.coefficients, compactRA) +
        vectorCell(item.r_b.coefficients, compactRB) +
        '<td>' + iconRow([
          mini('R\u226af', item.r_interlaces_f, 'R_d(f) interlaces f', 'R_d(f) does not interlace f', 'R_d(f) versus f: unavailable'),
          mini('M+', item.magic.all_nonnegative, 'all magic-basis coordinates are nonnegative', 'some magic-basis coordinates are negative'),
          mini('S\u2264', item.magic.left_leq_right, partialSumTooltip(item.magic), partialSumTooltip(item.magic), 'partial-sum data unavailable'),
        ]) + '</td>' +
        '</tr>';
    });
  } else {
    colspan = 4;
    const compactP = clean.every(item => allZeroSeq(item.coefficients));
    const compactCoords = clean.every(item => allZeroSeq(item.magic.coordinates));
    header = '<tr>' +
      '<th></th>' +
      '<th>' + headerLabel('deg', 'Degree after trimming trailing zeros.') + '</th>' +
      '<th>' + mathHeader('p(t)', 'Coefficient vector of the input polynomial ' + katexInline('p(t)') + '.') + '</th>' +
      '<th>' + mathHeader('c_j', 'Coordinates ' + katexInline('c_j') + ' in the basis ' + katexInline('t^j(1+t)^{d-j}') + '.') + '</th>' +
      '<th>' + headerLabel('Properties', 'M+: all magic-basis coordinates are nonnegative. S≤: the partial-sum inequalities L_j ≤ R_j hold; hover to see the first failing index when they do not. AI: alternatingly increasing.') + '</th>' +
      '</tr>';
    rows = lastDecompData.map((item, idx) => {
      if (item.error) {
        return '<tr><td class="idx">' + (idx + 1) + '</td><td colspan="' + colspan + '" class="err">' + esc(item.error) + '</td></tr>';
      }
      return '<tr>' +
        '<td class="idx">' + (idx + 1) + '</td>' +
        '<td class="compact">' + item.degree + '</td>' +
        vectorCell(item.coefficients, compactP) +
        vectorCell(item.magic.coordinates, compactCoords) +
        '<td>' + iconRow([
          mini('M+', item.magic.all_nonnegative, 'all magic-basis coordinates are nonnegative', 'some magic-basis coordinates are negative'),
          mini('S\u2264', item.magic.left_leq_right, partialSumTooltip(item.magic), partialSumTooltip(item.magic), 'partial-sum data unavailable'),
          mini('AI', item.alternatingly_increasing, 'alternatingly increasing', 'not alternatingly increasing'),
        ]) + '</td>' +
        '</tr>';
    });
  }

  const html = '<div class="decomp-subtabs">' + subtabs + '</div>' +
    '<div class="decomp-help">' + helpHtml + '</div>' +
    '<div class="decomp-sheet-wrap"><div class="decomp-sheet-scroll"><table class="decomp-sheet">' + header + rows.join('') + '</table></div></div>' +
    renderDecompCopyBar(activeDecompView, lastDecompData) +
    renderSheetOeisHelpers(activeDecompView, lastDecompData);
  setResult(html);
}
