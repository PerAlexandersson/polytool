// Run: node --test web/tests/page-assets.test.cjs
// The page is markup plus separate static assets: no inline style or script, every local
// asset exists, and all cache-busting versions agree with ASSET_VERSION.
const fs = require('node:fs');
const path = require('node:path');
const assert = require('node:assert/strict');
const {test} = require('node:test');
const {loadPage} = require('./page_source.cjs');

const page = loadPage();
const version = /const ASSET_VERSION = '([^']+)';/.exec(page.js);

test('index.html is markup only', () => {
  assert.doesNotMatch(page.markup, /<style\b/i);
  assert.ok(page.scripts.every(item => item.url), 'every script is an external file');
  assert.doesNotMatch(page.markup, /\sstyle=["'][^"']{200,}/);
});

test('local assets exist and ASSET_VERSION is defined once, in the first application script', () => {
  assert.ok(version, 'ASSET_VERSION is defined');
  assert.equal(page.js.match(/const ASSET_VERSION\b/g).length, 1);
  const appScripts = page.scripts.filter(item => !/(^|\/)pkg\//.test(item.url));
  assert.ok(appScripts[0].text.startsWith('const ASSET_VERSION = '), 'the first script defines the version');
  for (const item of [...page.stylesheets, ...page.scripts]) {
    assert.ok(fs.existsSync(item.file) || /(^|\/)pkg\//.test(item.url), `missing ${item.url}`);
  }
});

test('every local asset URL carries the current version query', () => {
  for (const item of [...page.stylesheets, ...page.scripts]) {
    assert.ok(item.url.endsWith('?v=' + version[1]), `${item.url} should end with ?v=${version[1]}`);
  }
});

test('every script under js/ is loaded by the page', () => {
  const loaded = new Set(page.scripts.map(item => item.url && item.url.split('?')[0]));
  for (const name of fs.readdirSync(path.join(page.dir, 'js'))) {
    assert.ok(loaded.has('js/' + name), `js/${name} is not referenced by index.html`);
  }
});

test('worker blobs and wasm loading use the shared version', () => {
  assert.match(page.js, /pkg\/polytool_web_bg\.wasm\?v=' \+ ASSET_VERSION/);
  assert.equal(page.js.match(/importScripts\('\$\{base\}pkg\/polytool_web\.js\?v=\$\{ASSET_VERSION\}'\)/g).length, 2);
});
