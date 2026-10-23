// Actual Chromium/Worker/WASM regressions. See README for optional prerequisites.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const {chromium} = require(process.env.POLYTOOL_PLAYWRIGHT || 'playwright');

const wasmDir = process.env.POLYTOOL_WASM_DIR;
assert.ok(wasmDir, 'Set POLYTOOL_WASM_DIR to a built no-modules WASM bundle');
const webDir = path.resolve(__dirname, '..');
const files = {
  '/': [path.join(webDir, 'index.html'), 'text/html'],
  '/favicon.svg': [path.join(webDir, 'favicon.svg'), 'image/svg+xml'],
  '/pkg/polytool_web.js': [path.join(wasmDir, 'polytool_web.js'), 'text/javascript'],
  '/pkg/polytool_web_bg.wasm': [path.join(wasmDir, 'polytool_web_bg.wasm'), 'application/wasm'],
};
let workerDelay = 0;
const server = http.createServer((req, res) => {
  const name = new URL(req.url, 'http://localhost').pathname;
  const file = files[name];
  if (!file) { res.writeHead(404); res.end(); return; }
  const send = () => { res.writeHead(200, {'Content-Type': file[1]}); res.end(fs.readFileSync(file[0])); };
  if (workerDelay && name.endsWith('.wasm')) setTimeout(send, workerDelay);
  else send();
});

async function run(page) {
  await page.click('#btnRec');
  await page.waitForFunction(() => !recurrenceRunning);
  return page.locator('#results-recurrence').innerText();
}
async function request(page, value) {
  await page.locator('#spaceAdvanced').evaluate(el => {el.open = true;});
  await page.fill('#spaceRequest', typeof value === 'string' ? value : JSON.stringify(value));
  return run(page);
}
const source = {coefficients: [['0', '1']], seeds: [['1']]};
const shape = {kind: 'box', length: 1, variable_degree: 1};

(async () => {
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const origin = 'http://127.0.0.1:' + server.address().port;
  let browser;
  try {
    browser = await chromium.launch({headless: true,
      ...(process.env.POLYTOOL_CHROMIUM ? {executablePath: process.env.POLYTOOL_CHROMIUM} : {}),
      args: ['--no-sandbox', '--disable-dev-shm-usage']});
    for (const width of [1366, 375]) {
      const context = await browser.newContext({viewport: {width, height: 900}});
      try {
        await context.grantPermissions(['clipboard-read', 'clipboard-write'], {origin});
        const page = await context.newPage();
        page.setDefaultTimeout(15000);
        const errors = [];
        page.on('pageerror', e => errors.push(String(e)));
        await page.goto(origin + '/#polynomial-sequence');
        await page.waitForFunction(() => wasmReady);
        await page.evaluate(() => document.fonts.ready);
        await page.click('.tab[data-tab="recurrence"]');
        await page.check('input[name="searchMode"][value="spaces"]');
        assert.equal(await page.locator('#tryDenominator').isDisabled(), true);
        assert.equal(await page.locator('#status').getAttribute('aria-live'), 'polite');

        assert.equal(await page.locator('#spaceSource').inputValue(), 'rows');
        assert.equal(await page.locator('#spaceRequest').isVisible(), false);
        assert.equal(await page.locator('#manualOpts').isVisible(), false);
        assert.equal(await page.locator('#legacyRecurrenceOptions').isVisible(), false);
        // Use the main input, including expanded polynomials and exact large integers.
        await page.fill('#input', '1\nx\nx^2\nx^3\nx^4\nx^5\nx^6\nx^7');
        await page.selectOption('#spaceLength', '1');
        await page.selectOption('#spaceDerivative', '0');
        await page.fill('#spaceOrigin', '0');
        let text = await run(page);
        assert.match(text, /fit the supplied rows exactly/);
        assert.match(text, /Recurrence found/);
        assert.equal(await page.evaluate(() => window._recurrenceSpaces.request.row_origin), 0);
        assert.equal(await page.evaluate(() => window._recurrenceSpaces.request.shapes.length), 1);
        await page.fill('#input', '1\ninvalid row');
        assert.equal(await page.evaluate(() => window._recurrenceSpaces), null);
        assert.match(await run(page), /Error:/);
        assert.equal(await page.evaluate(() => recurrenceRunning), false);
        await page.fill('#input', '9007199254740993\n18014398509481986\n36028797018963972\n72057594037927944\n144115188075855888');
        await run(page);
        assert.equal(await page.evaluate(() => window._recurrenceSpaces.request.rows[0][0]), '9007199254740993');
        await page.fill('#spaceOrigin', '-1');
        assert.match(await run(page), /whole number/);
        await page.fill('#spaceOrigin', '1');

        for (const example of ['default', 'hermite', 'automatic']) {
          await page.selectOption('#spaceSource', example);
          assert.equal(await page.locator('button[onclick^="copySpaceFormula"]').count(), 0);
          assert.equal(await page.locator('#spaceOriginField').isVisible(), false);
          const text = await run(page);
          assert.match(text, /Recurrence search results/);
          assert.match(text, /follow from the defining recurrence/);
          assert.doesNotMatch(text, /\bunfinished\b|kernel|projection|certificate|verify-space|frontier/i);
          assert.ok(await page.locator('#results-recurrence .katex').count() > 0, 'Formulas render as mathematics');
          await page.locator('button[onclick^="copySpaceFormula"]').first().click();
          const exportCheck = await page.evaluate(async () => {
            const copied = await navigator.clipboard.readText();
            return {same: copied === window._recurrenceSpaces.results.find(r => r.representatives.length).representatives[0].formula,
              replay: JSON.parse(wasm_bindgen.verify_recurrence_spaces(JSON.stringify(window._recurrenceSpaces)))};
          });
          assert.equal(exportCheck.same, true);
          assert.equal(exportCheck.replay.verified_completed_claims, true);
          assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1),
            'Page must fit the viewport, including mode controls and results');
          if (process.env.POLYTOOL_BROWSER_ARTIFACTS) {
            fs.mkdirSync(process.env.POLYTOOL_BROWSER_ARTIFACTS, {recursive: true});
            await page.screenshot({path: path.join(process.env.POLYTOOL_BROWSER_ARTIFACTS,
              `${width}-${example}.png`), fullPage: true});
            if (example === 'hermite') {
              await page.emulateMedia({colorScheme: 'dark'});
              await page.screenshot({path: path.join(process.env.POLYTOOL_BROWSER_ARTIFACTS,
                `${width}-hermite-dark.png`), fullPage: true});
              await page.emulateMedia({colorScheme: 'light'});
            }
          }
        }

        // Editing invalidates the old result and export, including while a job runs.
        await page.locator('#spaceAdvanced').evaluate(el => {el.open = true;});
        await page.fill('#spaceRequest', '{');
        assert.equal(await page.locator('button[onclick^="copySpaceFormula"]').count(), 0);
        assert.equal(await page.evaluate(() => window._recurrenceSpaces), null);
        assert.match(await run(page), /JSON|expected|EOF/i);
        assert.match(await page.locator('#status').innerText(), /Request rejected/);
        assert.match(await request(page, {rows: [['1'], ['0', '1'], ['0', '0', '1']], shapes: [shape]}),
          /fit the supplied rows exactly/);
        assert.match(await request(page, {source, shapes: [{kind: 'box', length: 1}]}),
          /No recurrence within these bounds/);
        const partial = await request(page, {source, limits: {max_matrix_entries: 2},
          shapes: [{kind: 'automatic_variable_degree', length: 1}]});
        assert.match(partial, /0 of 1 searches finished; 1 unfinished/);
        assert.match(partial, /does not tell us whether a recurrence exists/);
        assert.doesNotMatch(partial, /cover all possible x-degrees/);

        // Delay only the next worker's WASM fetch, making navigation/cancel deterministic.
        workerDelay = 700;
        await page.fill('#spaceRequest', JSON.stringify({source, shapes: [shape]}));
        await page.evaluate(() => { document.getElementById('results-oeis').textContent = 'Previous OEIS result'; });
        await page.click('#btnRec');
        await page.click('.tab[data-tab="oeis"]');
        await page.waitForFunction(() => !recurrenceRunning);
        assert.equal(await page.locator('#results-oeis').innerText(), 'Previous OEIS result');
        assert.match(await page.locator('#results-recurrence').innerText(), /Recurrence search results/);
        await page.click('.tab[data-tab="recurrence"]');
        await page.click('#btnRec');
        await page.click('#btnCancel');
        assert.match(await page.locator('#results-recurrence').innerText(), /cancelled/);
        assert.equal(await page.evaluate(() => recurrenceRunning), false);
        workerDelay = 0;
        assert.match(await run(page), /1 of 1 searches finished/);

        // A future render exception must be visible, with controls usable for retry.
        await page.evaluate(() => {
          window.savedSpaceRenderer = renderRecurrenceSpaces;
          renderRecurrenceSpaces = () => { throw new Error('injected renderer failure'); };
        });
        assert.match(await run(page), /Could not display recurrence result:.*injected renderer failure/);
        assert.equal(await page.locator('#btnRec').isDisabled(), false);
        await page.evaluate(() => { renderRecurrenceSpaces = window.savedSpaceRenderer; });
        assert.match(await run(page), /1 of 1 searches finished/);

        // Existing modes still work after spaces; no numeric data is mocked.
        await page.fill('#input', '1\n2\n4\n8\n16\n32\n64\n128');
        for (const mode of ['representations', 'adaptive', 'manual']) {
          await page.check(`input[name="searchMode"][value="${mode}"]`);
          const text = await run(page);
          assert.ok(text.length > 0);
          assert.doesNotMatch(text, /Could not display|Worker error|Searching for recurrence/);
          assert.match(await page.locator('#status').innerText(), mode === 'representations' ? /Comparison completed/ : /Found in/);
        }
        assert.equal(await page.locator('#tryDenominator').isDisabled(), false);
        assert.deepEqual(errors, []);
        console.log(`${width}px: guided controls, exact rows/origins, examples, copy, backend replay, partial search, errors, stale input, tab routing, cancellation and legacy modes passed`);
      } finally { await context.close(); }
    }
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
