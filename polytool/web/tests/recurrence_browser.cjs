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

        for (const example of ['default', 'hermite', 'automatic']) {
          if (example !== 'default') {
            await page.click(`button[onclick="loadRecurrenceSpaceExample('${example}')"]`);
            assert.equal(await page.locator('button[onclick^="copyRecurrenceSpaces"]').count(), 0);
          }
          const text = await run(page);
          assert.match(text, /Recurrence-space analysis/);
          assert.match(text, /conditional on the supplied source/);
          assert.match(text, /0 incomplete/);
          if (example === 'automatic') assert.match(text, /RHS X-degree cap: 3/);
          if (example === 'hermite') assert.match(text, /Polynomial-index propagation certificate: order 4/);
          await page.click('button[onclick^="copyRecurrenceSpaces"]');
          const exportCheck = await page.evaluate(async () => {
            const copied = await navigator.clipboard.readText();
            return {same: copied === JSON.stringify(window._recurrenceSpaces, null, 2),
              replay: JSON.parse(wasm_bindgen.verify_recurrence_spaces(copied))};
          });
          assert.equal(exportCheck.same, true);
          assert.equal(exportCheck.replay.verified_completed_claims, true);
          assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth + 1),
            'Page must fit the viewport, including mode controls and results');
          if (process.env.POLYTOOL_BROWSER_ARTIFACTS) {
            fs.mkdirSync(process.env.POLYTOOL_BROWSER_ARTIFACTS, {recursive: true});
            await page.screenshot({path: path.join(process.env.POLYTOOL_BROWSER_ARTIFACTS,
              `${width}-${example}.png`), fullPage: true});
          }
        }

        // Editing invalidates the old result and export, including while a job runs.
        await page.fill('#spaceRequest', '{');
        assert.equal(await page.locator('button[onclick^="copyRecurrenceSpaces"]').count(), 0);
        assert.equal(await page.evaluate(() => window._recurrenceSpaces), null);
        assert.match(await run(page), /JSON|expected|EOF/i);
        assert.match(await page.locator('#status').innerText(), /Request rejected/);
        assert.match(await request(page, {rows: [['1'], ['0', '1'], ['0', '0', '1']], shapes: [shape]}),
          /Exact fit to supplied rows; no all-rank conclusion/);
        assert.match(await request(page, {source, shapes: [{kind: 'box', length: 1}]}),
          /rows generated by the supplied source/);
        const partial = await request(page, {source, limits: {max_matrix_entries: 2},
          shapes: [{kind: 'automatic_variable_degree', length: 1}]});
        assert.match(partial, /0 of 1 shapes complete; 1 incomplete/);
        assert.match(partial, /Search at this cap is incomplete/);
        assert.doesNotMatch(partial, /completed search covers all RHS/);

        // Delay only the next worker's WASM fetch, making navigation/cancel deterministic.
        workerDelay = 700;
        await page.fill('#spaceRequest', JSON.stringify({source, shapes: [shape]}));
        await page.evaluate(() => { document.getElementById('results-oeis').textContent = 'Previous OEIS result'; });
        await page.click('#btnRec');
        await page.click('.tab[data-tab="oeis"]');
        await page.waitForFunction(() => !recurrenceRunning);
        assert.equal(await page.locator('#results-oeis').innerText(), 'Previous OEIS result');
        assert.match(await page.locator('#results-recurrence').innerText(), /Recurrence-space analysis/);
        await page.click('.tab[data-tab="recurrence"]');
        await page.click('#btnRec');
        await page.click('#btnCancel');
        assert.match(await page.locator('#results-recurrence').innerText(), /cancelled/);
        assert.equal(await page.evaluate(() => recurrenceRunning), false);
        workerDelay = 0;
        assert.match(await run(page), /1 of 1 shapes complete/);

        // A future render exception must be visible, with controls usable for retry.
        await page.evaluate(() => {
          window.savedSpaceRenderer = renderRecurrenceSpaces;
          renderRecurrenceSpaces = () => { throw new Error('injected renderer failure'); };
        });
        assert.match(await run(page), /Could not display recurrence result:.*injected renderer failure/);
        assert.equal(await page.locator('#btnRec').isDisabled(), false);
        await page.evaluate(() => { renderRecurrenceSpaces = window.savedSpaceRenderer; });
        assert.match(await run(page), /1 of 1 shapes complete/);

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
        console.log(`${width}px: examples, copy/replay, evidence, partial cap, errors, stale input, tab routing, cancellation, render recovery and legacy modes passed`);
      } finally { await context.close(); }
    }
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
