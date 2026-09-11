// Run: node --test polytool/web/input-normalization.test.cjs
// Optional POLYTOOL_HTML=/path/to/index.html checks the staged/live page.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const { test } = require('node:test');
const htmlPath = process.env.POLYTOOL_HTML || path.join(__dirname, 'index.html');
const html = fs.readFileSync(htmlPath, 'utf8');
const start = html.indexOf('function stripCommentOnlyLines(');
const end = html.indexOf('function copyRec(', start);
assert.ok(start >= 0 && end > start);
const context = vm.createContext({ document: { getElementById: () => ({ value: '' }) } });
vm.runInContext(html.slice(start, end), context);
const normalize = value => context.normalizePolynomialInput(value);

const oeis = '1;\n\n       1;\n\n       2;\n\n       4,       1;\n\n' +
  '      10,       5;\n\n      26,      26;\n\n      76,     117,      10;\n\n' +
  '     232,     540,     105;\n\n     764,    2445,     931;\n\n' +
  '    2620,   11338,    6909,    280;\n\n    9496,   53033,   48546,   4900;\n\n' +
  '   35696,  253826,  324753,  64295;\n\n  140152, 1235115, 2131855, 691075, 15400;';
const expectedOeis = [
  '1', '1', '2', '4, 1', '10, 5', '26, 26', '76, 117, 10', '232, 540, 105',
  '764, 2445, 931', '2620, 11338, 6909, 280', '9496, 53033, 48546, 4900',
  '35696, 253826, 324753, 64295', '140152, 1235115, 2131855, 691075, 15400',
].join('\n');

for (const [name, input, expected] of [
  ['exact user OEIS paste', oeis, expectedOeis],
  ['exact user Mathematica list', '{{1},{1,2},{1,2,3},{1,2,3,4}}', '1\n1, 2\n1, 2, 3\n1, 2, 3, 4'],
  ['JSON lists', '[[1],[1,2],[1,2,3],]', '1\n1, 2\n1, 2, 3'],
  ['Mathematica assignment', 'rows := {{1}, {1,2}};', '1\n1, 2'],
  ['Python tuple rows', '((1,), (1,2), (1,2,3))', '1\n1, 2\n1, 2, 3'],
  ['one-line OEIS', '1; 1, 2; 1, 2, 3;', '1\n1, 2\n1, 2, 3'],
  ['separate bracket rows', '{1}, {1,2}; [1,2,3]', '1\n1, 2\n1, 2, 3'],
  ['wrapped bracket rows', '[1,\n2];\n[3,\n4]', '1, 2\n3, 4'],
  ['whitespace, commas and signs', '+1  -2, 3;\n4,\t5,', '+1, -2, 3\n4, 5'],
  ['copy artifacts', '\ufeff1;\r\n\u00a0\u22122,\u202f3;\r4\u20285', '1\n-2, 3\n4\n5'],
  ['comment lines', '# copied OEIS\n1;\n  # another row\n2,3;', '1\n2, 3'],
  ['nested comments', '{{1},\n# comment\n{2,3}}', '1\n2, 3'],
  ['arbitrary precision', '{{9007199254740993}, {-1000000000000000000000000000007}}',
    '9007199254740993\n-1000000000000000000000000000007'],
  ['trailing zeros retained', '0; 1, 0, 0;', '0\n1, 0, 0'],
  ['expanded polynomials unchanged', '1 + 2t + t^2\nx**3 - 2*x', '1 + 2t + t^2\nx**3 - 2*x'],
  ['expanded Unicode minus', '1 \u2212 2t + t^2', '1 - 2t + t^2'],
  ['blank input', '', ''],
]) {
  test(name, () => {
    assert.equal(normalize(input), expected);
    assert.equal(normalize(expected), expected, 'normalization must be idempotent');
  });
}

for (const input of [
  '{{1}, {}, {2}}', '[1,,2]', '1,,2;', '1,2,...;', '1,2.5;', '1/2,3;',
  '[1,2}', '{{1},{2}', '{1}garbage', '1; alert(2);', '{{{1}}}', '[1]2',
  '{{1},2}', '[1;2]', ',1,2;', '[[1],[],[2]]',
]) {
  test('invalid input remains invalid: ' + input, () => assert.equal(normalize(input), input));
}

test('excessive nesting cannot overflow the call stack', () => {
  const input = '['.repeat(20000) + '1' + ']'.repeat(20000);
  assert.equal(normalize(input), input);
});

test('large coefficient row stays exact and complete', () => {
  const input = Array(20000).fill('9007199254740993').join(',') + ';';
  assert.equal(normalize(input), Array(20000).fill('9007199254740993').join(', '));
});

test('browser input boundary normalizes without rewriting the pasted text', () => {
  const textarea = { value: oeis };
  context.document.getElementById = id => {
    assert.equal(id, 'input');
    return textarea;
  };
  assert.equal(context.getInput(), expectedOeis);
  assert.equal(textarea.value, oeis);
});

test('page scripts remain syntactically valid', () => {
  for (const match of html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script>/gi)) {
    if (/type=["'](?:application\/ld\+json|importmap)["']/.test(match[1])) continue;
    if (match[2].trim()) new vm.Script(match[2]);
  }
});

// Exercise the actual deployed/staged WASM parser when explicitly supplied.
test('normalized rows reach real WASM without precision loss', {
  skip: !process.env.POLYTOOL_WASM_DIR,
}, async () => {
  const dir = process.env.POLYTOOL_WASM_DIR;
  const sandbox = vm.createContext({ TextDecoder, TextEncoder, WebAssembly });
  vm.runInContext(fs.readFileSync(path.join(dir, 'polytool_web.js'), 'utf8'), sandbox);
  const wasm = vm.runInContext('wasm_bindgen', sandbox);
  await wasm(fs.readFileSync(path.join(dir, 'polytool_web_bg.wasm')));
  const props = JSON.parse(wasm.check_properties(normalize(oeis)));
  assert.equal(props.length, 13);
  for (let i = 0; i < props.length; i++) {
    assert.ok(!props[i].error, JSON.stringify(props[i]));
    assert.deepEqual(props[i].coefficients.map(String), expectedOeis.split('\n')[i].split(', '));
  }
  const huge = JSON.parse(wasm.check_properties(normalize('9007199254740993;')));
  assert.equal(String(huge[0].coefficients[0]), '9007199254740993');
  const nested = JSON.parse(wasm.check_properties(normalize('{{1},{1,2},{1,2,3},{1,2,3,4}}')));
  assert.deepEqual(nested.map(p => p.coefficients.map(String)),
    [['1'], ['1', '2'], ['1', '2', '3'], ['1', '2', '3', '4']]);
});
