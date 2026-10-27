// --- Transform tab: coordinate changes, implicit relations, roots, theorems ---

// Evidence labels used by polytool reports: [short label, tone, explanation].
// Tones order the strength of support; unknown labels are shown verbatim.
const EVIDENCE_VOCABULARY = {
  certified_all_indices: ['Certified for all indices', 'certified', 'Proved for every index of the domain from the defining recurrence.'],
  completed_recurrence_induction: ['Proved by induction', 'certified', 'A closed induction over the transported recurrence covers every index.'],
  certified_family_invariant: ['Certified for the family', 'certified', 'Transferred to every index by a map theorem whose conditions are certified.'],
  certified_by_reconstruction: ['Proved by reconstruction', 'certified', 'An exact reconstruction certificate proves the fitted rule for every index.'],
  equivalent_description: ['Equivalent description', 'certified', 'The seeds and relations determine exactly the transformed sequence.'],
  equivalent_generating_description: ['Equivalent generating description', 'certified', 'Certified annihilation, sufficient initial or boundary data and a certified uniqueness argument: the description determines exactly this sequence.'],
  certified_complete_space: ['Certified complete space', 'certified', 'Proved for every index as a member of a certified complete space of relations.'],
  certified_annihilating_relation: ['Certified annihilating relation', 'transported', 'The relation holds for every index of the domain, but it is not certified to determine the sequence.'],
  defining: ['Defining rule', 'transported', 'The source rule itself; every other result is relative to it.'],
  certified: ['Certified', 'certified', 'Proved exactly.'],
  exact_hypothesis: ['Checked exactly', 'certified', 'This hypothesis of the theorem is verified exactly for the recurrence.'],
  transported_identity: ['Transported identity', 'transported', 'An exact consequence of the defining recurrence: it holds for every index whenever that rule does.'],
  conditional_on_source: ['Follows from the defining recurrence', 'transported', 'Holds for every index whenever the defining recurrence does.'],
  certified_consequence: ['Certified consequence', 'transported', 'Every relation holds for all indices, but together they need not determine the sequence.'],
  conditional_theorem_match: ['Conditional theorem match', 'theorem', 'The hypotheses of a published theorem are checked exactly. The conclusion applies once the recurrence is known to hold for the intended sequence.'],
  paper_proof_not_lean_checked: ['Paper proof, not Lean-checked', 'theorem', 'The theorem has a paper proof that has not been formally verified.'],
  checked_supplied_rows: ['Checked on supplied rows', 'checked', 'Exact on the listed rows only; it says nothing about later indices.'],
  finite_exact_observation: ['Exact on listed rows', 'checked', 'Exact root isolation of the listed rows only, never a statement about the whole family.'],
  fixed_cap: ['Fixed degree cap only', 'checked', 'A finite check up to the stated degree cap. It is not a proof for larger degrees and does not make a claim certified.'],
  preserves_up_to_cap: ['Preserved up to the cap', 'checked', 'Every real-rooted polynomial of degree at most the cap keeps real-rootedness; nothing is said about larger degrees.'],
  fails_at_cap: ['Fails at the cap', 'fail', 'An exact witness: a real-rooted polynomial of the cap degree whose image is not real-rooted.'],
  inconclusive: ['Inconclusive', 'neutral', 'The test draws no conclusion.'],
  exact_fixed_polynomial: ['Exact for a fixed polynomial', 'certified', 'A fixed polynomial checked exactly: a proof for that polynomial.'],
  assumed: ['Assumed', 'assumption', 'Supplied in the request and not proved.'],
  not_checked: ['Not checked', 'assumption', 'This hypothesis was not checked.'],
  violated: ['Violated', 'fail', 'Shown false exactly; the implication does not apply.'],
  attempt_failed: ['Attempt failed', 'fail', 'This interval did not close the induction.'],
  consistent_with_rows: ['Consistent with rows', 'checked', 'Holds on the observed rows only.'],
  refuted_on_rows: ['Refuted on rows', 'fail', 'An exact check failed on an observed row.'],
  finite_fit: ['Finite fit', 'fit', 'Fitted to finitely many rows and checked on the rest. It may fail on later rows.'],
  exact_fit: ['Finite fit', 'fit', 'Fits the supplied rows exactly. It may fail on later rows.'],
  finite_evidence_only: ['Finite evidence only', 'fit', 'Supported by finitely many rows; not a proof for all n.'],
  finite_evidence: ['Finite evidence', 'fit', 'Supported by finitely many rows; not a proof for all n.'],
  supplied_assumption: ['Your assumption', 'assumption', 'Taken from the request and checked on the listed rows only.'],
  proposed_interval: ['Proposed window', 'assumption', 'A candidate interval for a proof, not itself proved.'],
  failed_on_rows: ['Fails on rows', 'fail', 'An exact check failed on a supplied row.'],
  inconsistent_seeds: ['Inconsistent seeds', 'fail', 'The seeds contradict the relations.'],
};
const EVIDENCE_TONES = [
  ['certified', 'certified_all_indices', 'proved for every index.'],
  ['transported', 'transported_identity', 'exact consequence of the defining recurrence; as reliable as that rule.'],
  ['theorem', 'conditional_theorem_match', 'hypotheses checked exactly; the conclusion still needs the recurrence for all n.'],
  ['checked', 'checked_supplied_rows', 'exact on the listed rows only.'],
  ['fit', 'finite_fit', 'fitted to finitely many rows; may fail later.'],
  ['assumption', 'proposed_interval', 'supplied or proposed, not proved.'],
  ['fail', 'failed_on_rows', 'an exact check failed.'],
];

function humanizeCode(code) {
  const text = String(code == null ? 'unlabelled' : code).replace(/_/g, ' ');
  return text.charAt(0).toUpperCase() + text.slice(1);
}

function evidenceInfo(code) {
  const entry = Object.prototype.hasOwnProperty.call(EVIDENCE_VOCABULARY, code) ? EVIDENCE_VOCABULARY[code] : null;
  if (entry) return {code, label: entry[0], tone: entry[1], help: entry[2]};
  return {code, label: humanizeCode(code), tone: 'neutral', help: 'Label reported by polytool: ' + String(code)};
}

function evidenceBadge(code) {
  const info = evidenceInfo(code);
  return '<span class="evidence evidence-' + info.tone + '" data-evidence="' + escAttr(code) +
    '" title="' + escAttr(info.help) + '">' + esc(info.label) + '</span>';
}

function evidenceLegendHtml() {
  return '<details class="evidence-legend"><summary>What the labels mean</summary><ul>' +
    EVIDENCE_TONES.map(([, code, text]) => '<li>' + evidenceBadge(code) + esc(text) + '</li>').join('') +
    '</ul></details>';
}

// --- Formulas: polytool's plain-text rules to TeX (coefficients unchanged) ---

function texIndex(lag, sign) {
  if (!lag || lag === '0') return 'n';
  return 'n' + (sign || '+') + lag;
}

function texDerivative(name, index, order) {
  return name + (order <= 3 ? "'".repeat(order) : '^{(' + order + ')}') + '_{' + index + '}';
}

function ruleTextTex(text, options) {
  const opts = options || {};
  let tex = String(text == null ? '' : text)
    .replace(/\(1\) (?=[A-Za-z])/g, '')
    .replace(/D\^(\d+)\s*\(?([A-Za-z])\[n(?:\+(\d+))?\]\)?/g, (_, order, name, lag) =>
      texDerivative(name, texIndex(lag), Number(order)))
    .replace(/\bD\s+([A-Za-z])\[n(?:\+(\d+))?\]/g, (_, name, lag) => texDerivative(name, texIndex(lag), 1))
    .replace(/([A-Za-z])\[n(?:\+(\d+))?\]/g, (_, name, lag) => name + '_{' + texIndex(lag) + '}')
    .replace(/([A-Za-z])\^\{?\((\d+)\)\}?\((n(?:[+-]\d+)?)\)/g, '$1^{($2)}_{$3}')
    .replace(/([A-Za-z])('+)\((n(?:[+-]\d+)?)\)/g, '$1$2_{$3}')
    .replace(/([A-Za-z])\((n(?:[+-]\d+)?)\)/g, '$1_{$2}')
    .replace(/\^(\d+)/g, '^{$1}')
    .replace(/(\d+)\/(\d+)/g, '\\frac{$1}{$2}')
    .replace(/X/g, 'x');
  if (opts.name) tex = tex.replace(/\bP(?=['^_])/g, opts.name);
  if (opts.variable) {
    // Rename the variable t outside TeX control words such as \left.
    tex = tex.split(/(\\[A-Za-z]+)/).map((part, i) => i % 2 ? part : part.replace(/t/g, opts.variable)).join('');
  }
  return tex;
}

function formulaHtml(tex) {
  return '<div class="recurrence-display">' + katexInline(tex) + '</div>';
}

// A scalar polytool.recurrence.v1 record, typeset with the library's LaTeX.
function recordFormulaTex(record, fallbackText, name, variable) {
  let latex = null;
  try {
    if (typeof wasm_bindgen !== 'undefined' && wasm_bindgen.recurrence_latex) {
      const shown = JSON.parse(wasm_bindgen.recurrence_latex(JSON.stringify(record)));
      if (shown.latex) latex = shown.latex;
    }
  } catch (_) {}
  return ruleTextTex(latex || fallbackText || '', {name: name || 'P', variable: variable || 'x'});
}

function ruleFormulaHtml(rule, name, variable) {
  if (!rule) return '';
  if (rule.scalar) return formulaHtml(recordFormulaTex(rule.scalar, rule.text, name, variable));
  if (rule.text) {
    return String(rule.text).split('; ').map(line => formulaHtml(ruleTextTex(line, {variable: variable || 'x'}))).join('');
  }
  return '';
}

// Gamma polynomials are written in z = x/(1 + x)^2, everything else in x.
function transformVariable(report) {
  return report && report.transformation && report.transformation.kind === 'gamma' ? 'z' : 'x';
}

function rationalTex(value) {
  const text = String(value).trim();
  const match = /^(\d+)\/(\d+)$/.exec(text);
  if (match) return '\\frac{' + match[1] + '}{' + match[2] + '}';
  return /^\d+$/.test(text) ? text : '\\text{?}';
}

function rowPolyTex(coeffs, variable) {
  const v = variable || 'x';
  const terms = [];
  (coeffs || []).forEach((raw, k) => {
    let text = String(raw).trim();
    const negative = text.startsWith('-');
    if (negative) text = text.slice(1);
    if (/^0+(\/\d+)?$/.test(text)) return;
    const monomial = k === 0 ? '' : (k === 1 ? v : v + '^{' + k + '}');
    const body = k > 0 && text === '1' ? monomial : rationalTex(text) + monomial;
    terms.push((terms.length ? (negative ? ' - ' : ' + ') : (negative ? '-' : '')) + body);
  });
  return terms.length ? terms.join('') : '0';
}

// An endpoint `value + per_index * n`, written as -1-n, n/2, 1/2+3n/2 or -n.
function endpointText(endpoint) {
  const value = String(endpoint.value).trim();
  const slope = String(endpoint.per_index == null ? '0' : endpoint.per_index).trim();
  if (!/^-?\d+(\/\d+)?$/.test(value)) return '?';
  const match = /^(-?)(\d+)(?:\/(\d+))?$/.exec(slope);
  if (!match) return '?';
  if (/^0+$/.test(match[2])) return value;
  const numerator = /^0*1$/.test(match[2]) ? '' : match[2].replace(/^0+(?=\d)/, '');
  const denominator = match[3] && !/^0*1$/.test(match[3]) ? '/' + match[3].replace(/^0+(?=\d)/, '') : '';
  const term = match[1] + numerator + 'n' + denominator;
  if (/^-?0+$/.test(value)) return term;
  return value + (term.startsWith('-') ? '' : '+') + term;
}

function intervalText(interval) {
  if (!interval) return '(−∞, ∞)';
  const lower = interval.lower ? (interval.lower.closed === false ? '(' : '[') + endpointText(interval.lower) : '(−∞';
  const upper = interval.upper ? endpointText(interval.upper) + (interval.upper.closed === false ? ')' : ']') : '∞)';
  return lower + ', ' + upper;
}

function indexRangeText(first, last) {
  if (first == null) return '';
  if (last == null) return 'for every n ≥ ' + first;
  return first === last ? 'at n = ' + first : 'for n = ' + first + '…' + last;
}

function rootPropertyText(property, interval) {
  if (property === 'real_rooted') return 'is real-rooted';
  if (property === 'real_rooted_in') return 'is real-rooted with every root in ' + intervalText(interval);
  if (property === 'real_roots_in') return 'has every real root in ' + intervalText(interval);
  return humanizeCode(property).toLowerCase();
}

function componentNames(report) {
  const components = report.transformation && report.transformation.components;
  const names = (components || []).map(c => /^[A-Za-z]$/.test(String(c.name)) ? String(c.name) : 'Q');
  return names.length ? names : ['Q'];
}

function subjectTex(subject, component, names) {
  return subject === 'source' ? 'P_n' : ((names && names[component || 0]) || 'Q') + '_n';
}

function metaHtml(text) {
  return text ? '<p class="report-meta">' + esc(text) + '</p>' : '';
}

function listItems(items, render) {
  return (items || []).map(item => '<div class="report-item">' + render(item) + '</div>').join('');
}

function rootEnclosureText(root) {
  const multiplicity = root.multiplicity > 1 ? ' ×' + root.multiplicity : '';
  return (root.lower === root.upper ? root.lower : '(' + root.lower + ', ' + root.upper + ')') + multiplicity;
}

function rowsTableHtml(title, columns, rows, firstIndex, limit, variables) {
  if (!rows || !rows.length) return '';
  const shown = rows.slice(0, limit || 6);
  let html = '<details><summary>' + esc(title) + ' (' + rows.length + ' rows)</summary><div class="report-rows"><table><thead><tr><th>n</th>' +
    columns.map(c => '<th>' + katexInline(c) + '</th>').join('') + '</tr></thead><tbody>';
  shown.forEach((row, i) => {
    html += '<tr><td>' + esc(firstIndex + i) + '</td>' + row.map((cell, j) => '<td>' + katexInline(rowPolyTex(cell, variables && variables[j])) + '</td>').join('') + '</tr>';
  });
  html += '</tbody></table></div>';
  if (rows.length > shown.length) html += metaHtml('Showing the first ' + shown.length + ' rows; the exported JSON has all of them.');
  return html + '</details>';
}

function conditionsHtml(conditions) {
  return listItems(conditions, c => evidenceBadge(c.status) + '<strong>' + esc(c.condition) + '</strong>' + metaHtml(c.detail));
}

function reportActionButton(label, call) {
  return '<button class="run-small" onclick="' + escAttr(call) + '">' + esc(label) + '</button>';
}

// --- Report renderers ---

// The Jensen test of an index-free diagonal weight: evidence `fixed_cap`, never a proof for larger degrees.
function jensenCapHtml(jensen) {
  const cap = Number(jensen.cap);
  const verdicts = {
    preserves_up_to_cap: 'The weight keeps every real-rooted polynomial of degree at most ' + cap + ' real-rooted.',
    fails_at_cap: '(1 + x)^' + cap + ' is real-rooted but its image is not: the weight does not preserve real-rootedness.',
    inconclusive: 'The test is inconclusive at degree ' + cap + '.',
  };
  const known = Object.prototype.hasOwnProperty.call(verdicts, jensen.status);
  return '<div class="cap-note" data-evidence-block="fixed_cap">' + evidenceBadge(jensen.evidence || 'fixed_cap') +
    evidenceBadge(jensen.status) + '<strong>Jensen test at degree ' + esc(cap) + '</strong>' +
    '<p class="cap-warning">Finite check only: not a proof for degree ' + esc(cap + 1) + ' or more.</p>' +
    metaHtml(known ? verdicts[jensen.status] : '') + metaHtml(jensen.statement) +
    '<details><summary>Jensen polynomial and its roots</summary>' +
    formulaHtml('J_{' + cap + '}(x) = ' + rowPolyTex(jensen.jensen_polynomial)) +
    '<p class="report-text">' + esc((jensen.roots || []).map(rootEnclosureText).join(', ') || 'no real roots') + '</p></details></div>';
}

function renderRootDomainSection(rd, names, index) {
  const rank = {certified: 0, transported: 1, theorem: 2, checked: 3, fit: 4, assumption: 5, fail: 6, neutral: 7};
  const claims = (rd.claims || []).slice().sort((a, b) =>
    rank[evidenceInfo(a.evidence).tone] - rank[evidenceInfo(b.evidence).tone]);
  let html = '<h4>Root location</h4>' +
    metaHtml('Optional analysis. The algebraic results above never depend on it.');
  html += listItems(claims, claim =>
    evidenceBadge(claim.evidence) + '<strong>' + katexInline(subjectTex(claim.subject, claim.component, names)) + ' ' +
    esc(rootPropertyText(claim.property, claim.interval)) + '</strong> ' + esc(indexRangeText(claim.first_index, claim.last_index)) +
    metaHtml(claim.derivation) +
    ((claim.conditional_on || []).length ? '<p class="report-meta">Depends on: ' + esc(claim.conditional_on.join('; ')) + '</p>' : ''));
  if (!claims.length) html += metaHtml('No root-location claims.');
  [['Recurrence step', rd.recurrence_theorem], ['Defining rule of Pₙ', rd.source_rule_theorem],
    ['Certified successor rule of ' + names[0] + 'ₙ', rd.fitted_rule_theorem]].forEach(([title, theorem]) => {
    if (!theorem) return;
    html += '<details><summary>' + esc(title) + ': ' +
      esc(theorem.matched ? 'matches ' + humanizeCode(theorem.theorem).toLowerCase() : 'no match') + '</summary>' +
      (theorem.form ? '<p class="report-text">' + esc(theorem.form) + '</p>' : '') + metaHtml(theorem.statement) +
      listItems(theorem.attempts, attempt => evidenceBadge(attempt.success ? 'completed_recurrence_induction' : 'attempt_failed') +
        esc(intervalText(attempt.interval) + ' (' + humanizeCode(attempt.origin).toLowerCase() + ')') + metaHtml(attempt.detail)) +
      (theorem.reason ? metaHtml(theorem.reason) : '') + '</details>';
  });
  if (rd.jensen_fixed_cap) html += jensenCapHtml(rd.jensen_fixed_cap);
  if ((rd.implications || []).length) {
    html += '<details><summary>Map theorems (' + rd.implications.length + ')</summary>' +
      listItems(rd.implications, imp => '<strong>' + esc(humanizeCode(imp.theorem)) + '</strong> · ' + esc(humanizeCode(imp.direction).toLowerCase()) +
        '<p>' + esc((imp.source_property || '') + (imp.target_property ? ' ⇔ ' + imp.target_property : '')) + '</p>' +
        metaHtml(imp.statement) +
        listItems(imp.hypotheses, h => evidenceBadge(h.status) + esc(h.condition) + metaHtml(h.detail)) +
        (imp.converse ? metaHtml('Converse: ' + imp.converse) : '')) + '</details>';
  }
  if ((rd.assumption_checks || []).length) {
    html += '<details><summary>Assumption checks (' + rd.assumption_checks.length + ')</summary><pre class="report-text">' +
      esc(JSON.stringify(rd.assumption_checks, null, 2)) + '</pre></details>';
  }
  if ((rd.observations || []).length) {
    html += '<details><summary>Exact root enclosures of the observed rows (' + rd.observations.length + ')</summary>' +
      '<div class="report-rows"><table><thead><tr><th>Row</th><th>Real-rooted</th><th>Roots</th></tr></thead><tbody>' +
      rd.observations.map(o => '<tr><td>' + katexInline(subjectTex(o.subject, o.component, names).replace('_n', '_{' + o.index + '}')) +
        '</td><td>' + (o.real_rooted ? 'yes' : 'no') + '</td><td class="report-text">' +
        esc((o.roots || []).map(rootEnclosureText).join(', ') || 'none') + '</td></tr>').join('') +
      '</tbody></table></div>' + metaHtml('Each enclosure contains exactly one root; equal endpoints are exact roots.') + '</details>';
  }
  (rd.notes || []).forEach(note => { html += metaHtml(note); });
  if ((rd.conventions || []).length) {
    html += '<details><summary>Conventions</summary>' + rd.conventions.map(metaHtml).join('') + '</details>';
  }
  return html;
}

function renderTransformReport(report, index) {
  const t = report.transformation || {};
  const names = componentNames(report);
  let html = '<h4>Coordinate change</h4><p class="report-text">' + esc(t.row_map || t.kind) + '</p>' +
    metaHtml('Inverse: ' + (t.inverse || 'none')) +
    metaHtml('Index domain: ' + (t.index_domain || '') + (t.degree_schedule ? '; ' + t.degree_schedule : ''));
  if (transformVariable(report) === 'z') {
    html += metaHtml(names[0] + 'ₙ is the gamma polynomial in z = x/(1 + x)², so the transported rule and the rows below use z.');
  }
  html += '<div class="report-item">' + evidenceBadge(report.reconstruction_verified ? 'checked_supplied_rows' : 'failed_on_rows') +
    esc(report.reconstruction_verified ? 'The inverse reproduces all ' + (report.source_rows || []).length + ' source rows exactly.'
      : 'The inverse did not reproduce every source row.') + '</div>';
  html += conditionsHtml(report.conditions);

  html += '<h4>Transported recurrence</h4>';
  const tr = report.transported;
  if (tr) {
    html += evidenceBadge(tr.evidence) + ruleFormulaHtml(tr, names[0], transformVariable(report)) +
      metaHtml('Holds from index ' + tr.valid_from_index + (tr.valid_to_index != null ? ' to ' + tr.valid_to_index : '') +
        '; checked on ' + tr.checked_targets + ' transformed rows' + (tr.holds_on_transformed_rows ? '.' : ', where it FAILED.'));
    if (tr.scalar) html += '<div class="report-actions">' + reportActionButton('Check interlacing theorems for this recurrence', 'checkReportTheorems(' + index + ",'transported')") + '</div>';
  } else if (report.transport_error) {
    html += metaHtml('Not transported: ' + report.transport_error);
  } else {
    html += metaHtml('No defining recurrence was supplied, so only the rows were mapped.');
  }

  const search = report.independent_search;
  if (search) {
    html += '<h4>Independent fit in the new coordinates</h4>' + evidenceBadge(search.evidence) +
      '<span class="report-meta">' + esc(humanizeCode(search.status)) + '</span>' + ruleFormulaHtml(search, names[0], transformVariable(report)) +
      metaHtml(search.reason) +
      (search.holds_on_transformed_rows === false ? metaHtml('The fit fails on a transformed row.') : '');
    if (search.scalar) html += '<div class="report-actions">' + reportActionButton('Check interlacing theorems for the fit', 'checkReportTheorems(' + index + ",'fit')") + '</div>';
  }
  if (report.root_domain) html += renderRootDomainSection(report.root_domain, names, index);

  const origin = report.row_origin || 0;
  const states = report.transformed_states || [];
  html += rowsTableHtml('Rows', ['P_n'].concat(names.map(n => n + '_n')),
    (report.source_rows || []).map((row, i) => [row].concat(states[i] || [])), origin, undefined,
    ['x'].concat(names.map(() => transformVariable(report))));
  html += '<div class="report-actions">' +
    (report.source_rule ? reportActionButton('Check theorems for the source recurrence', 'checkReportTheorems(' + index + ",'source')") : '') +
    reportActionButton('Exact roots of the source rows', 'isolateReportRows(' + index + ",'source')") +
    reportActionButton('Exact roots of ' + names[0] + 'ₙ', 'isolateReportRows(' + index + ",'transformed')") + '</div>';
  return html;
}

function relationHtml(relation) {
  return evidenceBadge(relation.evidence) + '<span class="report-meta">' + esc(humanizeCode(relation.role || '')) +
    (relation.first_base != null ? ', base n ≥ ' + relation.first_base : '') + '</span>' +
    formulaHtml(ruleTextTex(relation.text)) + metaHtml(relation.note);
}

function renderImplicitReport(report, index) {
  const t = report.transformation || {};
  const assessment = report.assessment || {};
  const description = report.description || {};
  let html = '<h4>Coefficient map</h4><p class="report-text">' + esc(t.row_map || t.kind) + '</p>' +
    metaHtml('Inverse: ' + (t.inverse || 'none')) + conditionsHtml(report.map_conditions);
  html += '<h4>Generating description</h4>' + evidenceBadge(assessment.status) +
    '<span class="report-meta">' + esc((description.seeds || []).length + ' seed rows from n = ' + description.first_index +
      (assessment.seeds_consistent === false ? '; the seeds are inconsistent' : '')) + '</span>';
  html += listItems(description.relations, relationHtml);
  if ((description.exceptional_index_conditions || []).length) {
    html += '<p class="report-meta">Coefficient conditions at exceptional indices:</p>' +
      listItems(description.exceptional_index_conditions, relationHtml);
  }
  if ((description.boundary_conditions || []).length) {
    html += '<p class="report-meta">Supplied boundary conditions:</p>' + listItems(description.boundary_conditions, relationHtml);
  }
  if (assessment.uniqueness) {
    const u = assessment.uniqueness;
    html += '<details><summary>Why the description is unique</summary>' + evidenceBadge(u.status) + metaHtml(u.statement) +
      (u.reasons || []).map(metaHtml).join('') + '</details>';
  }
  html += '<h4>Transported relations</h4>';
  html += listItems(report.transports, transport => relationHtml(transport.relation) +
    metaHtml('From ' + humanizeCode(transport.source).toLowerCase() + '; checked at ' + transport.checked_bases + ' base indices' +
      (transport.holds_on_transformed_rows ? '.' : ', where it FAILED.')) +
    (transport.clearing ? '<p>' + evidenceBadge(transport.clearing.exceptional && transport.clearing.exceptional.status) +
      esc('Common multiplier H = ' + transport.clearing.multiplier_text) + '</p>' +
      '<details><summary>Clearing statement</summary>' + metaHtml(transport.clearing.statement) + '</details>' : '') +
    (transport.reconstruction ? '<details><summary>Reconstruction certificate</summary>' + metaHtml(transport.reconstruction.statement) + '</details>' : ''));
  if ((report.additional_relations || []).length) {
    html += '<h4>Additional row relations</h4>' + listItems(report.additional_relations, relationHtml);
  }
  const search = report.independent_search;
  if (search) {
    html += '<h4>Fitted successor rule</h4>' + evidenceBadge(search.evidence) +
      '<span class="report-meta">' + esc(humanizeCode(search.status)) + '</span>' +
      (search.relation ? formulaHtml(ruleTextTex(search.relation.text)) : '') + metaHtml(search.reason);
    if (search.scalar) html += '<div class="report-actions">' + reportActionButton('Check interlacing theorems for the fit', 'checkReportTheorems(' + index + ",'fit')") + '</div>';
  }
  if (report.root_domain) html += renderRootDomainSection(report.root_domain, ['Q'], index);
  const origin = report.row_origin || 0;
  const transformed = report.transformed_rows || [];
  html += rowsTableHtml('Rows', ['P_n', 'Q_n'], (report.source_rows || []).map((row, i) => [row, transformed[i] || []]), origin);
  html += '<div class="report-actions">' +
    (report.source_rule ? reportActionButton('Check theorems for the source recurrence', 'checkReportTheorems(' + index + ",'source')") : '') +
    reportActionButton('Exact roots of the source rows', 'isolateReportRows(' + index + ",'source')") +
    reportActionButton('Exact roots of Qₙ', 'isolateReportRows(' + index + ",'transformed')") + '</div>';
  return html;
}

function theoremName(code) {
  return {weighted_renewal_window: 'Weighted renewal window',
    fixed_delay_predecessor_window: 'Fixed-delay predecessor window'}[code] || humanizeCode(code);
}

function referenceHtml(reference) {
  const text = String(reference || '');
  if (/^https:\/\/[^\s"'<>]+$/.test(text)) return '<a href="' + escAttr(text) + '" target="_blank" rel="noopener">' + esc(text) + '</a>';
  return esc(text);
}

function renderTheoremReport(report) {
  const evidence = report.finite_evidence || {};
  let html = evidenceBadge(evidence.status) + '<span class="report-meta">' + esc(humanizeCode(evidence.source)) + '</span>' +
    formulaHtml(recordFormulaTex(report.recurrence_json, report.recurrence_text, 'P')) + metaHtml(evidence.note);
  const matches = report.theorem_matches || [];
  if (!matches.length) html += '<p><strong>No supported theorem matches this recurrence.</strong></p>';
  html += listItems(matches, m => evidenceBadge(m.status) + evidenceBadge(m.proof_status) +
    '<p><strong>' + esc(theoremName(m.theorem)) + '</strong></p><p class="report-text">' + esc(m.family) + '</p>' +
    '<p>' + esc(m.conclusion) + '</p>' +
    (m.window ? '<p>Interlacing window: <strong>' + esc(m.window.inequality) + '</strong> · ' + esc(m.window.statement) + '</p>' : '') +
    '<p class="report-meta">Still to show:</p><ul>' + (m.remaining_obligations || []).map(o => '<li>' + esc(o) + '</li>').join('') + '</ul>' +
    '<details><summary>Hypotheses and boundary rows</summary>' +
    listItems(m.hypotheses, h => evidenceBadge('exact_hypothesis') + esc(h.statement) + metaHtml(humanizeCode(h.check))) +
    (m.boundary ? '<div class="report-item">' + evidenceBadge(m.boundary.rows_match_canonical ? 'checked_supplied_rows' : 'failed_on_rows') +
      esc(humanizeCode(m.boundary.status)) + metaHtml(m.boundary.contract) + '</div>' : '') +
    '<p class="report-meta">Reference: ' + referenceHtml(m.reference) + '</p></details>');
  if ((report.rejected_theorems || []).length) {
    html += '<details><summary>Theorems that do not apply (' + report.rejected_theorems.length + ')</summary>' +
      listItems(report.rejected_theorems, r => '<strong>' + esc(theoremName(r.theorem)) + '</strong>' + metaHtml(r.reason)) + '</details>';
  }
  return html;
}

function renderRootEnclosureReport(report) {
  const items = report.items || [];
  const interval = report.request && report.request.interval;
  let html = evidenceBadge('finite_exact_observation') +
    '<span class="report-meta">Exact for these polynomials only.' + (interval ? ' Interval: ' + esc(intervalText(interval)) + '.' : '') + '</span>';
  html += '<div class="report-rows"><table><thead><tr><th>n</th><th>Degree</th><th>Real-rooted</th><th>Roots (enclosures)</th>' +
    (interval ? '<th>All roots in interval</th>' : '') + '</tr></thead><tbody>' +
    items.map(item => '<tr><td>' + esc(item.index) + '</td><td>' + esc(item.degree == null ? 'zero' : item.degree) +
      '</td><td>' + (item.real_rooted ? 'yes' : 'no') + ' (' + esc(item.real_root_count) + ' real)</td><td class="report-text">' +
      esc((item.roots || []).map(rootEnclosureText).join(', ') || 'none') + '</td>' +
      (interval ? '<td>' + (item.real_rooted_in_interval ? 'yes' : 'no') + '</td>' : '') + '</tr>').join('') +
    '</tbody></table></div>' + metaHtml('Each enclosure contains exactly one root and neither endpoint is a root; equal endpoints are exact rational roots.');
  return html;
}

function renderRecurrenceSpacesReport(report) {
  return listItems(report.results, result => evidenceBadge(result.evidence) +
    '<span class="report-meta">' + esc(humanizeCode(result.status)) + '</span>' +
    (result.representatives || []).map(rep => formulaHtml(ruleTextTex(rep.formula))).join('') + metaHtml(result.reason));
}

// --- Representation search report (polytool.recurrence-representation-search.v1) ---

function weightLabel(weight) {
  const w = weight || {};
  if (w.kind === 'factorial_power') return 'coefficient of xᵏ ÷ ' + (w.power === 1 ? 'k!' : '(k!)^' + w.power);
  if (w.kind === 'binomial') return 'coefficient of xᵏ × binom(k+' + w.m + ', ' + w.m + ')';
  if (w.kind === 'inverse_binomial') return 'coefficient of xᵏ ÷ binom(k+' + w.m + ', ' + w.m + ')';
  return humanizeCode(w.kind || 'weight');
}

// A short name for one allowed map of the request; all text is escaped by the caller.
function transformationShortLabel(t) {
  const kind = String((t && t.kind) || 'map');
  if (kind === 'diagonal') return 'Diagonal: ' + weightLabel(t.weight);
  if (kind === 'affine') return 'Substitute x ↦ ' + t.a + 'x + ' + t.b;
  if (kind === 'multiply_power') return 'Multiply by x^' + t.power;
  if (kind === 'divide_power') return 'Divide by x^' + t.power;
  if (kind === 'derivative_antecedent') return 'Derivative antecedent, g = ' + (Array.isArray(t.g) ? t.g : []).join(', ');
  return humanizeCode(kind);
}

function representationMapLabel(report, index) {
  const list = (report.request && report.request.transformations) || [];
  return Object.prototype.hasOwnProperty.call(list, index) ? transformationShortLabel(list[index]) : 'map ' + index;
}

function representationChainLabel(report, rep) {
  if (!rep.chain || !rep.chain.length) return 'source';
  return rep.chain.map(i => representationMapLabel(report, i)).join(' → ');
}

function descriptionCostText(cost) {
  const c = cost || {};
  return 'shift ' + c.shift_length + ' · derivative order ' + c.derivative_order + ' · coefficient degree ' +
    c.coefficient_degree + ' · index degree ' + c.index_degree + ' · ' + c.component_count +
    (c.component_count === 1 ? ' component' : ' components') + ' · reconstruction cost ' + c.reconstruction_cost;
}

function descriptionFormulaHtml(description, name) {
  if (description.scalar) return formulaHtml(recordFormulaTex(description.scalar, description.text, name));
  return String(description.text || '').split('; ')
    .map(line => formulaHtml(ruleTextTex(line, {variable: 'x', name}))).join('');
}

function repComponentName(rep) {
  const first = rep && rep.components && rep.components[0];
  return /^[A-Za-z]$/.test(String(first)) ? String(first) : 'Q';
}

function showRepresentation(card, id) {
  const element = document.getElementById('rep-' + card + '-' + id);
  if (!element) return;
  element.open = true;
  if (element.scrollIntoView) element.scrollIntoView({block: 'nearest'});
}

function representationLinkButton(card, id) {
  return '<button class="rep-link" onclick="showRepresentation(' + Number(card) + ',' + Number(id) + ')">#' + esc(id) + '</button>';
}

function repEdgeOutcome(edge, report, card) {
  const label = representationMapLabel(report, edge.transformation);
  const detail = edge.detail ? '<span class="report-meta"> ' + esc(edge.detail) + '</span>' : '';
  if (edge.outcome === 'duplicate_map') {
    return '<li class="rep-edge rep-edge-duplicate"><span class="rep-outcome">Duplicate map</span> ' + esc(label) +
      (edge.representation != null ? ' has the same normal form as ' + representationLinkButton(card, edge.representation) : '') + detail + '</li>';
  }
  if (edge.outcome === 'rejected') {
    return '<li class="rep-edge rep-edge-rejected"><span class="rep-outcome">Rejected</span> ' + esc(label) + detail + '</li>';
  }
  if (edge.outcome === 'skipped_budget') {
    return '<li class="rep-edge rep-edge-skipped"><span class="rep-outcome">Skipped: budget</span> ' + esc(label) + detail + '</li>';
  }
  return '<li class="rep-edge"><span class="rep-outcome">' + esc(humanizeCode(edge.outcome)) + '</span> ' + esc(label) + detail + '</li>';
}

// The tree of maps: source, then every attempted map application in queue order.
function representationTreeHtml(report, card) {
  const reps = new Map((report.representations || []).map(rep => [rep.id, rep]));
  const children = new Map();
  (report.edges || []).forEach(edge => {
    if (!children.has(edge.parent)) children.set(edge.parent, []);
    children.get(edge.parent).push(edge);
  });
  const frontier = new Set(report.depth_frontier || []);
  const visited = new Set();
  const nodeHtml = (id, edge, level) => {
    const rep = reps.get(id);
    if (!rep || visited.has(id) || level > 12) return '';
    visited.add(id);
    const best = (rep.descriptions || [])[0];
    let html = '<li class="rep-node">' + representationLinkButton(card, id) + ' <strong>' +
      esc(edge ? representationMapLabel(report, edge.transformation) : 'Source') + '</strong> ';
    if (best) html += evidenceBadge(best.label);
    if (rep.equal_sequence_as) {
      html += '<span class="rep-tag rep-tag-equal">same sequence as #' + esc(rep.equal_sequence_as.representation) + '</span>';
    }
    if (frontier.has(id)) html += '<span class="rep-tag rep-tag-frontier">depth limit: not expanded</span>';
    else if (!rep.expanded && !rep.equal_sequence_as) html += '<span class="rep-tag rep-tag-frontier">not expanded</span>';
    const kids = (children.get(id) || []).map(next => next.outcome === 'new_representation' && next.representation != null
      ? nodeHtml(next.representation, next, level + 1) : repEdgeOutcome(next, report, card)).join('');
    return html + (kids ? '<ul>' + kids + '</ul>' : '') + '</li>';
  };
  return '<ul class="rep-tree">' + nodeHtml(0, null, 0) + '</ul>';
}

function repBudgetHtml(report) {
  const events = report.budget || [];
  const frontier = report.depth_frontier || [];
  let html = '';
  if (events.length || report.status === 'budget_exhausted') {
    html += '<div class="rep-budget" role="note"><strong>Budget exhausted: this search is incomplete.</strong>' +
      '<p class="report-meta">A missing description or representation below is not evidence that it does not exist.</p><ul>' +
      events.map(event => '<li><span class="rep-outcome">' + esc(humanizeCode(event.limit)) + '</span>' +
        (event.representation != null ? ' at representation ' + esc(event.representation) : '') + ': ' + esc(event.detail) + '</li>').join('') +
      '</ul></div>';
  }
  if (frontier.length) {
    const depth = report.request && report.request.limits && report.request.limits.max_depth;
    html += '<div class="rep-budget rep-depth" role="note">Depth limit' + (depth != null ? ' (' + esc(depth) + ')' : '') +
      ' reached: the allowed maps were not applied to ' + frontier.map(id => '#' + esc(id)).join(', ') + '.</div>';
  }
  return html;
}

function renderRepresentationDescription(description, rep) {
  const sameCode = description.label === description.evidence;
  let html = '<div class="report-item rep-description">' + evidenceBadge(description.label) +
    (sameCode ? '' : evidenceBadge(description.evidence)) +
    '<span class="report-meta">' + esc(humanizeCode(description.kind)) + '</span>' +
    descriptionFormulaHtml(description, repComponentName(rep)) +
    '<p class="report-meta">' + esc(descriptionCostText(description.cost)) + '</p>';
  const found = [description.provenance].concat(description.also_found_by || []).filter(Boolean);
  if (found.length) html += '<p class="report-meta">Found by: ' + esc(found.join(', ')) + '</p>';
  if (description.description) {
    const d = description.description;
    html += '<p class="report-meta">' + esc((d.seeds || []).length + ' seed rows from n = ' + d.first_index +
      ((d.exceptional_index_conditions || []).length ? '; ' + d.exceptional_index_conditions.length + ' exceptional-index condition(s)' : '') +
      ((d.boundary_conditions || []).length ? '; ' + d.boundary_conditions.length + ' boundary condition(s)' : '')) + '</p>';
  }
  html += metaHtml(description.note);
  const u = description.assessment && description.assessment.uniqueness;
  if (u) {
    html += '<details><summary>Why the description is unique</summary>' + evidenceBadge(u.status) + metaHtml(u.statement) +
      (u.reasons || []).map(metaHtml).join('') + '</details>';
  } else if (description.assessment) {
    html += '<p class="report-meta">' + evidenceBadge(description.assessment.status) + '</p>';
  }
  return html + '</div>';
}

// Small searches open every representation; larger ones open on request (tree links).
const REPRESENTATION_OPEN_LIMIT = 6;

function renderRepresentationNode(report, rep, card) {
  const best = (rep.descriptions || [])[0];
  const origin = rep.row_origin || 0;
  const columns = (rep.components || ['Q']).map(c => c + '_n');
  const open = (report.representations || []).length <= REPRESENTATION_OPEN_LIMIT ? ' open' : '';
  let html = '<details class="rep-detail" id="rep-' + Number(card) + '-' + Number(rep.id) + '"' + open + '><summary><strong>Representation ' +
    esc(rep.id) + '</strong> · ' + esc(representationChainLabel(report, rep)) + ' · ' + esc((rep.descriptions || []).length) +
    ' description(s) ' + (best ? evidenceBadge(best.label) : '') + '</summary>';
  const segments = (rep.map && rep.map.segments) || [];
  html += '<h4>Coordinate map</h4>' + (rep.map && rep.map.identity
    ? '<p class="report-meta">Identity: the source rows.</p>'
    : segments.map(segment => '<p class="report-text">' + esc(segment) + '</p>').join('')) +
    metaHtml('Index domain: ' + ((rep.map && rep.map.index_domain) || '') + '; components: ' + (rep.components || []).join(', '));
  if (rep.depth > 0) {
    html += '<div class="report-item">' + evidenceBadge(rep.reconstruction_verified ? 'checked_supplied_rows' : 'failed_on_rows') +
      esc(rep.reconstruction_verified ? 'Inverting the whole chain reproduces all source rows exactly.'
        : 'Inverting the chain did not reproduce every source row.') + '</div>';
  }
  if (rep.equal_sequence_as) {
    html += '<div class="rep-equal"><strong>Same sequence as representation ' + esc(rep.equal_sequence_as.representation) +
      '.</strong> ' + esc(rep.equal_sequence_as.certificate) + '</div>';
  }
  html += conditionsHtml(rep.map_conditions);
  (rep.notes || []).forEach(note => { html += metaHtml(note); });
  html += '<h4>Descriptions, simplest first</h4>';
  html += (rep.descriptions || []).length
    ? rep.descriptions.map(description => renderRepresentationDescription(description, rep)).join('')
    : metaHtml('No description was found in this representation.');
  if (rep.fits) {
    html += metaHtml('Fit menu: ' + humanizeCode(rep.fits.status).toLowerCase() + '; ' + rep.fits.configurations_considered +
      ' configurations considered, ' + rep.fits.fits_found + ' fit(s) found' + (rep.fits.reason ? '. ' + rep.fits.reason : '.'));
  }
  (rep.steps || []).forEach(step => { if (step.error) html += metaHtml('A transport step failed: ' + step.error); });
  if ((rep.steps || []).length || (rep.spaces || []).length) {
    html += metaHtml('Embedded certificates: ' + (rep.steps || []).length + ' transport step(s) and ' + (rep.spaces || []).length +
      ' complete-space report(s); Replay checks each with its own verifier.');
  }
  html += rowsTableHtml('Rows', columns, rep.states || [], origin);
  return html + '</details>';
}

function renderRepresentationSearchReport(report, index) {
  const reps = report.representations || [];
  const source = report.source || {};
  const complete = report.status === 'complete';
  let html = '<p><span class="rep-status ' + (complete ? 'rep-status-complete' : 'rep-status-partial') + '">' +
    (complete ? 'Complete within the requested limits' : 'Budget exhausted: incomplete') + '</span> ' +
    '<span class="report-meta">' + esc(reps.length) + ' representation(s), ' + esc((report.ranked || []).length) +
    ' description(s), ' + esc(report.steps_used == null ? '?' : report.steps_used) + ' backend call(s)</span></p>';
  html += repBudgetHtml(report);
  html += '<h4>Source</h4><p class="report-meta">' + esc({defining_rule: 'Supplied defining rule', oeis: 'Bundled OEIS family',
    rows_only: 'Rows only'}[source.kind] || humanizeCode(source.kind)) + '; ' + esc(source.row_count == null ? '?' : source.row_count) +
    ' rows from n = ' + esc(source.row_origin == null ? '?' : source.row_origin) + '.</p>';
  if (source.kind === 'rows_only') {
    html += metaHtml('There is no defining rule, so every description is a finite fit to these rows.');
  } else if (source.rule) {
    html += formulaHtml(recordFormulaTex(source.rule, null, 'P'));
  }
  const maps = (report.request && report.request.transformations) || [];
  html += '<p class="report-meta">Allowed maps: ' + esc(maps.map(transformationShortLabel).join('; ') || 'none') + '.</p>';
  html += '<h4>Representations and maps</h4>' +
    '<p class="report-meta">Each representation applies every allowed map. Maps are merged only when their normal forms are certified equal, never by comparing rows.</p>' +
    representationTreeHtml(report, index);
  const ranked = (report.ranked || []).slice(0, 6);
  if (ranked.length) {
    html += '<h4>Simplest descriptions overall</h4>' + listItems(ranked, item => {
      const rep = reps.find(r => r.id === item.representation);
      const description = rep && (rep.descriptions || [])[item.description];
      return representationLinkButton(index, item.representation) + evidenceBadge(item.label) +
        (description ? descriptionFormulaHtml(description, repComponentName(rep)) : '<p class="report-text">' + esc(item.text) + '</p>') +
        '<p class="report-meta">' + esc(descriptionCostText(item.cost)) + '</p>';
    });
    if ((report.ranked || []).length > ranked.length) html += metaHtml('Showing the first ' + ranked.length + ' of ' + report.ranked.length + '; every representation below lists its own.');
  }
  html += '<h4>Each representation</h4>' + reps.map(rep => renderRepresentationNode(report, rep, index)).join('');
  html += metaHtml('Certification is relative to the defining source, the maps actually enumerated and the finite supports searched.');
  if ((report.conventions || []).length) {
    html += '<details><summary>Conventions</summary>' + report.conventions.map(metaHtml).join('') + '</details>';
  }
  return html;
}

// Fallback for report types without a dedicated view: every object carrying
// an `evidence` label is listed with its text, so no label is lost.
function renderGenericEvidenceReport(report) {
  const found = [];
  const visit = (value, path) => {
    if (found.length >= 200 || value == null || typeof value !== 'object') return;
    if (Array.isArray(value)) { value.forEach((item, i) => visit(item, path + '[' + i + ']')); return; }
    if (typeof value.evidence === 'string') {
      found.push({path, evidence: value.evidence,
        text: value.text || value.formula || value.statement || value.recurrence_text || value.status || ''});
    }
    Object.keys(value).forEach(key => { if (key !== 'request') visit(value[key], path ? path + '.' + key : key); });
  };
  visit(report, '');
  let html = metaHtml('No dedicated view for ' + (report.schema || 'this report') + ' yet; its evidence labels are listed as reported.');
  html += listItems(found, item => evidenceBadge(item.evidence) + '<span class="report-meta">' + esc(item.path) + '</span>' +
    (item.text ? '<p class="report-text">' + esc(item.text) + '</p>' : ''));
  if (!found.length) html += metaHtml('The report contains no evidence labels.');
  return html;
}

// Add a renderer here when a new report schema is supported.
const REPORT_RENDERERS = {
  'polytool.recurrence-transform.v1': ['Coordinate transformation', renderTransformReport],
  'polytool.recurrence-implicit.v1': ['Implicit relations', renderImplicitReport],
  'polytool.recurrence-theorems.v1': ['Interlacing theorems', renderTheoremReport],
  'polytool.root-enclosures.v1': ['Exact root enclosures', renderRootEnclosureReport],
  'polytool.recurrence-spaces.v1': ['Recurrence search', renderRecurrenceSpacesReport],
  'polytool.recurrence-representation-search.v1': ['Representation search', renderRepresentationSearchReport],
};

function replayBannerHtml(entry) {
  const replay = entry.replay;
  if (replay === 'running') return '<p class="replay-pending">Replaying the report&hellip;</p>';
  if (replay && replay.ok) return '<p class="replay-ok">Replay passed: polytool recomputed this report exactly from its embedded request. The labels below keep their meaning; replay does not strengthen them.</p>';
  if (replay && replay.error) return '<p class="replay-fail">Replay failed: ' + esc(replay.error) + '. Do not rely on this report.</p>';
  if (entry.origin === 'imported') return '<p class="replay-pending">Imported report, not yet replayed: its content is shown as written in the file.</p>';
  return '';
}

function renderReportCard(entry, index) {
  const report = entry.report || {};
  const known = Object.prototype.hasOwnProperty.call(REPORT_RENDERERS, report.schema) ? REPORT_RENDERERS[report.schema] : null;
  const title = entry.title || (known ? known[0] : humanizeCode(report.schema || 'report'));
  let body = known ? known[1](report, index) : renderGenericEvidenceReport(report);
  if (entry.replay && entry.replay.error) body = '<details><summary>Show the unverified content</summary>' + body + '</details>';
  return '<section class="report-card" data-schema="' + escAttr(report.schema || '') + '"><h3>' + esc(title) + '</h3>' +
    replayBannerHtml(entry) + body +
    '<div class="report-actions">' +
    reportActionButton('Copy report JSON', 'copyTransformReport(this,' + index + ')') +
    reportActionButton('Download JSON', 'downloadTransformReport(' + index + ')') +
    reportActionButton('Replay (verify)', 'replayTransformReport(' + index + ')') + '</div></section>';
}

function renderTransformResults() {
  if (!transformEntries.length) { setResult('', 'transform'); return; }
  setResult(evidenceLegendHtml() + transformEntries.map(renderReportCard).join(''), 'transform');
}

function transformReportSummary(report) {
  if (report.schema === 'polytool.recurrence-transform.v1') {
    return report.transported ? 'Recurrence transported.' : (report.transport_error ? 'Rows mapped; recurrence not transported.' : 'Rows mapped.');
  }
  if (report.schema === 'polytool.recurrence-implicit.v1') return evidenceInfo(report.assessment && report.assessment.status).label + '.';
  if (report.schema === 'polytool.recurrence-theorems.v1') return (report.theorem_matches || []).length + ' conditional theorem match(es).';
  if (report.schema === 'polytool.root-enclosures.v1') return (report.items || []).length + ' polynomial(s) isolated.';
  if (report.schema === 'polytool.recurrence-representation-search.v1') {
    const certified = (report.ranked || []).filter(item => item.label !== 'finite_fit').length;
    return (report.representations || []).length + ' representation(s), ' + certified + ' certified description(s)' +
      (report.status === 'complete' ? '.' : '; budget exhausted, the search is incomplete.');
  }
  return 'Report ready.';
}
