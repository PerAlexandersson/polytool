import assert from 'node:assert/strict';
import fs from 'node:fs';
import vm from 'node:vm';

const source = fs.readFileSync(new URL('../index.html', import.meta.url), 'utf8');

function between(start, end) {
  const startIndex = source.indexOf(start);
  const endIndex = source.indexOf(end, startIndex);
  assert.notEqual(startIndex, -1, `missing ${start}`);
  assert.notEqual(endIndex, -1, `missing ${end}`);
  return source.slice(startIndex, endIndex);
}

const helpers = [
  between('function serializeCoeffRow', 'function propsForRows'),
  between('function normalizeIntegerString', 'function gammaPolyText'),
  between('function isZeroValue', 'function copyData'),
].join('\n');

const context = vm.createContext({ BigInt });
vm.runInContext(helpers, context);

const call = expression => vm.runInContext(expression, context);
const huge = '1000000000000000000000000000000';

assert.equal(
  context.normalizeIntegerString('-000' + huge),
  '-' + huge
);
assert.equal(
  call("coeffsToPolyText(['9007199254740993', '-100000000000000000000', '1'], 't')"),
  '9007199254740993 - 100000000000000000000t + t^2'
);
assert.deepEqual(
  Array.from(call("firstDifferences(['9007199254740993', '9007199254740995', '10000000000000000000'])")),
  [
    '2',
    (10000000000000000000n - 9007199254740995n).toString(),
  ]
);
assert.equal(
  call("sumIntegerRow(['100000000000000000000', '2', '-3'], false)"),
  '99999999999999999999'
);
assert.equal(
  call("sumIntegerRow(['100000000000000000000', '2', '-3'], true)"),
  '99999999999999999995'
);
assert.deepEqual(Array.from(call("stripLeadingZeros(['0', '00', '9007199254740993'])")), [
  '9007199254740993',
]);
assert.equal(
  call("serializeCoeffRow(['9007199254740993', '100000000000000000000'])"),
  '9007199254740993, 100000000000000000000'
);

assert.doesNotMatch(source, /Math\.abs\(c\)/);
assert.doesNotMatch(source, /Number\.isFinite\(Number\(value\)\)/);
assert.doesNotMatch(source, /coefficients\.reduce\(/);

console.log('string-safe coefficient helpers passed');
