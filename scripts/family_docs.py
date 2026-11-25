"""Render the curated family site from the installed binary's canonical data.

No mathematical generation is implemented here: all rows, draft blocks and
b-files come from polytool. All interpolated metadata is escaped as text.
"""
import html
import json
from pathlib import Path
import subprocess


def call(binary, *arguments):
    return json.loads(subprocess.run([binary, "families", *arguments, "--json"],
                                    check=True, text=True, stdout=subprocess.PIPE).stdout)


def escaped(value):
    return html.escape(str(value), quote=True)


def link(url, label):
    return f'<a href="{escaped(url)}">{escaped(label)}</a>'


def page(title, body, prefix=".."):
    return ('<!doctype html><html lang="en"><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{escaped(title)} — Polytool</title>'
            f'<link rel="stylesheet" href="{prefix}/style.css">'
            f'<nav aria-label="Documentation"><a href="{prefix}/index.html">Polytool</a>'
            '<a href="index.html">Families</a></nav><main>' + body + '</main></html>')


def proof(property, references):
    result = f'<p>Status: <strong>{escaped(property["status"])}</strong>'
    if property.get("as_of"):
        result += f' (as of {escaped(property["as_of"])})'
    result += '</p>'
    if property["statement"]:
        result += f'<p>{escaped(property["statement"])}</p>'
    if property.get("theorem"):
        ref = next(r for r in references if r["key"] == property["theorem"])
        result += '<p>Theorem: ' + link(ref["url"], ref["citation"] + ', ' + ref["locator"]) + '</p>'
    if property.get("lean"):
        lean = property["lean"]
        result += '<p>Public Lean theorem: ' + link(lean["permalink"], lean["declaration"]) + '</p>'
        result += '<p>Manual attestation. The statement above records the scope of the identification.</p>'
    elif property["status"] == 'lean_proved':
        result += '<p>Canonical external Lean attestation; no public source mapping curated here.</p>'
    if property.get("paper_proof"):
        result += f'<p>Paper proof: {escaped(property["paper_proof"])}</p>'
    if property.get("lag"):
        result += f'<p>Lag: {escaped(property["lag"])}</p>'
    if property.get("finite_range"):
        finite = property["finite_range"]
        result += (f'<p>Finite evidence, displayed indices {finite["first"]}–{finite["last"]}: '
                   f'{escaped(finite["evidence"])} Finite success is not an all-index proof.</p>')
    return result


def copy_block(label, content):
    return (f'<h3>{escaped(label)}</h3><button type="button" class="copy">Copy</button>'
            f'<pre><code>{escaped(content)}</code></pre>')


def render_families(binary, output, lean_export=None):
    arguments = ["validate"]
    if lean_export:
        arguments += ["--lean-export", str(lean_export)]
    validation = call(binary, *arguments)
    listing = call(binary, "list")["families"]
    destination = output / "families"
    (destination / "bfiles").mkdir(parents=True)
    (destination / "proposals").mkdir()
    index_rows, records = [], []
    for summary in listing:
        family_id = summary["id"]
        record = call(binary, "show", family_id)["family"]
        rows = call(binary, "generate", family_id, "--count", "6")["rows"]
        draft = call(binary, "export-oeis", family_id, "--count", "12")
        records.append(record)
        bfile = f'bfiles/{family_id}.txt'
        proposal = f'proposals/{family_id}.toml'
        (destination / bfile).write_text(draft["bfile"])
        (destination / proposal).write_text(draft["proposal_toml"])
        status = record["real_rootedness"]["status"]
        lean = record["real_rootedness"].get("lean")
        models = ', '.join(m["family"] + ' (' + m["replay_status"] + ')' for m in record["models"])
        cells = [link(f'{family_id}.html', family_id), escaped(record["title"]),
                 escaped(record["degree_pattern"] or 'Not curated'), escaped(status),
                 escaped(record["interlacing"]["status"]),
                 link(lean["permalink"], 'Lean') if lean else '—', escaped(models or '—')]
        index_rows.append(f'<tr data-status="{escaped(status)}">' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')
        body = f'<h1>{escaped(family_id)}: {escaped(record["title"])}</h1>'
        body += '<h2>Definition</h2><p>' + escaped(record["definition"] or 'Not yet curated; see the canonical generator and OEIS entry.') + '</p>'
        if record["oeis"]["id"]:
            body += '<p>' + link('https://oeis.org/' + record["oeis"]["id"], 'OEIS ' + record["oeis"]["id"]) + '</p>'
        for heading, field in [('Formula', 'formula'), ('Generating function', 'generating_function')]:
            body += f'<h2>{heading}</h2><pre>{escaped(record[field] or "Not curated")}</pre>'
        body += '<h2>First rows</h2><p>Generated exactly by Polytool at this build, coefficients in ascending powers of t.</p><pre>'
        body += escaped('\n'.join(f'{record["first_row"] + i}: ' + ', '.join(row) for i, row in enumerate(rows))) + '</pre>'
        body += '<h2>Canonical generator</h2><pre>' + escaped(json.dumps(record["generator"], indent=2)) + '</pre>'
        if record["generator_status"]:
            body += f'<p>Generator validation status: {escaped(record["generator_status"])}; separate from root/proof status.</p>'
        body += '<h2>Real-rootedness</h2>' + proof(record["real_rootedness"], record["references"])
        body += '<h2>Interlacing</h2>' + proof(record["interlacing"], record["references"])
        body += '<h2>References</h2><ul>'
        for ref in record["references"]:
            identifiers = ', '.join(x for x in [ref.get('doi'), ref.get('arxiv')] if x)
            body += '<li>' + link(ref["url"], ref["key"] + ': ' + ref["citation"] + '; ' + ref["locator"]) + escaped(' ' + identifiers) + '</li>'
        body += '</ul><h2>Models</h2><pre>' + escaped(json.dumps(record["models"], indent=2) if record["models"] else 'Not curated') + '</pre>'
        body += '<h2>Related families</h2><ul>'
        for related in record["related"]:
            body += '<li>' + link(related["id"] + '.html', related["id"]) + ': ' + escaped(related["relation"]) + '</li>'
        body += '</ul><h2>OEIS submission</h2><p>State: ' + escaped(record["submission"]["state"])
        body += '; date: ' + escaped(record["submission"]["date"] or '—') + '; assigned A-number: ' + escaped(record["submission"]["a_number"] or '—') + '.</p>'
        body += '<p>Last OEIS search: ' + escaped(record["oeis"]["checked"] or 'Not recorded') + '.</p><pre>' + escaped('\n'.join(record["oeis"]["queries"])) + '</pre>'
        body += '<p>Draft only: review all fields before submission. A000000 denotes an unassigned A-number. '
        body += 'The TOML skeleton needs author and repository-local generator fields filled in.</p>'
        body += copy_block('OEIS blocks', draft["blocks"]) + copy_block('Proposal TOML skeleton', draft["proposal_toml"])
        body += '<p>' + link(bfile, f'Download b-file ({draft["bfile_terms"]} terms from 12 rows)') + '; ' + link(proposal, 'Download proposal skeleton') + '.</p>'
        body += '<h2>External links</h2><ul>' + ''.join('<li>' + link(url, url) + '</li>' for url in record["links"]) + '</ul>'
        body += '<script src="families.js"></script>'
        (destination / f'{family_id}.html').write_text(page(family_id, body))
    statuses = sorted({r["real_rootedness"]["status"] for r in records})
    body = ('<h1>Polynomial sequence families</h1><p>Curated OEIS entries and submission-bound families. '
            'This catalog records definitions, generators and proof attestations; the atlas remains a discovery census. '
            'Object-indexed graph or matroid classes belong in the external open-problems catalog.</p>'
            '<p>Proof links are manual public commit-pinned attestations. A future proof-owner export will '
            'supply tags for links in both directions from one source of truth.</p>'
            '<label>Search <input id="query" type="search"></label> '
            '<label>Real-rootedness status <select id="status"><option value="">All statuses</option>')
    body += ''.join(f'<option>{escaped(status)}</option>' for status in statuses)
    body += f'</select></label><p id="count">{len(records)} families</p><div class="family-table"><table id="families"><thead><tr>'
    body += ''.join(f'<th scope="col">{name}</th>' for name in ['ID', 'Title', 'Degree pattern', 'Real-rootedness', 'Interlacing', 'Lean', 'Models'])
    body += '</tr></thead><tbody>' + ''.join(index_rows) + '</tbody></table></div>'
    body += '<p>' + link('catalog.json', 'Download validated metadata') + '</p><script src="families.js"></script>'
    (destination / 'index.html').write_text(page('Families', body))
    (destination / 'catalog.json').write_text(json.dumps(records, indent=2) + '\n')
    (destination / 'families.js').write_text('''"use strict";
const table = document.querySelector("#families");
if (table) {
  const query = document.querySelector("#query"), status = document.querySelector("#status");
  function filter() {
    let visible = 0;
    for (const row of table.tBodies[0].rows) {
      row.hidden = !(row.textContent.toLowerCase().includes(query.value.toLowerCase()) &&
        (!status.value || row.dataset.status === status.value));
      if (!row.hidden) visible++;
    }
    document.querySelector("#count").textContent = visible + " families";
  }
  query.addEventListener("input", filter);
  status.addEventListener("change", filter);
}
for (const button of document.querySelectorAll(".copy")) {
  button.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(button.nextElementSibling.textContent);
      button.textContent = "Copied";
    } catch (_) { button.textContent = "Select and copy the text below"; }
  });
}
''')
    print(f'Families: {len(records)} generated pages; {validation["count"]} validated records')
