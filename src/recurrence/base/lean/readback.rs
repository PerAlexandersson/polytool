//! Parser and typed evaluator for the emitted Lean expression subset.

use super::*;

#[derive(Clone, Debug)]
enum Expression {
    Number(BigInt),
    Index,
    Variable,
    Row(Box<Expression>),
    Constant(Box<Expression>),
    Cast(Box<Expression>, Kind),
    Negate(Box<Expression>),
    Binary(char, Box<Expression>, Box<Expression>),
    Derivative(Box<Expression>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Natural,
    Real,
    Polynomial,
}

#[derive(Debug)]
enum Value {
    Natural(BigInt),
    Real(BigRational),
    Polynomial(Vec<BigRational>),
}

struct Parser<'a> {
    tokens: Vec<String>,
    position: usize,
    name: &'a str,
}

impl<'a> Parser<'a> {
    fn new(text: &str, name: &'a str) -> Result<Self, RecurrenceError> {
        let mut tokens = Vec::new();
        let mut chars = text.chars().peekable();
        while let Some(c) = chars.next() {
            if c.is_whitespace() {
                continue;
            }
            if c.is_ascii_alphanumeric() || c == '_' || c == '\'' {
                let mut token = String::from(c);
                while chars
                    .peek()
                    .is_some_and(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '\'')
                {
                    token.push(chars.next().expect("peeked character"));
                }
                tokens.push(token);
            } else if "()+-*/^:.ℕℝ".contains(c) {
                tokens.push(c.to_string());
            } else {
                return Err(invalid(format!("unsupported Lean token {c}")));
            }
        }
        Ok(Self {
            tokens,
            position: 0,
            name,
        })
    }
    fn peek(&self) -> &str {
        self.tokens.get(self.position).map_or("", String::as_str)
    }
    fn take(&mut self) -> Result<String, RecurrenceError> {
        let token = self
            .tokens
            .get(self.position)
            .cloned()
            .ok_or_else(|| invalid("unexpected end of Lean expression"))?;
        self.position += 1;
        Ok(token)
    }
    fn require(&mut self, token: &str) -> Result<(), RecurrenceError> {
        if self.take()? == token {
            Ok(())
        } else {
            Err(invalid(format!("expected Lean token {token}")))
        }
    }
    fn parse(mut self) -> Result<Expression, RecurrenceError> {
        let expression = self.binary(0)?;
        if self.position != self.tokens.len() {
            return Err(invalid("unparsed Lean expression suffix"));
        }
        Ok(expression)
    }
    fn binary(&mut self, minimum: u8) -> Result<Expression, RecurrenceError> {
        let mut left = self.atom()?;
        loop {
            let (operator, precedence) = match self.peek() {
                "+" => ('+', 1),
                "-" => ('-', 1),
                "*" => ('*', 2),
                "/" => ('/', 2),
                "^" => ('^', 3),
                _ => break,
            };
            if precedence < minimum {
                break;
            }
            self.take()?;
            let right = self.binary(if operator == '^' {
                precedence
            } else {
                precedence + 1
            })?;
            left = Expression::Binary(operator, Box::new(left), Box::new(right));
        }
        Ok(left)
    }
    fn atom(&mut self) -> Result<Expression, RecurrenceError> {
        let token = self.take()?;
        let mut expression = match token.as_str() {
            "(" => {
                let value = self.binary(0)?;
                let value = if self.peek() == ":" {
                    self.take()?;
                    let kind = match self.take()?.as_str() {
                        "ℕ" => Kind::Natural,
                        "ℝ" => Kind::Real,
                        _ => return Err(invalid("unsupported Lean type ascription")),
                    };
                    Expression::Cast(Box::new(value), kind)
                } else {
                    value
                };
                self.require(")")?;
                value
            }
            "-" => Expression::Negate(Box::new(self.binary(3)?)),
            "n" => Expression::Index,
            "X" => Expression::Variable,
            "C" => Expression::Constant(Box::new(self.atom()?)),
            row if row == self.name => Expression::Row(Box::new(self.atom()?)),
            number if number.bytes().all(|c| c.is_ascii_digit()) => Expression::Number(
                number
                    .parse()
                    .map_err(|_| invalid("invalid Lean numeral"))?,
            ),
            _ => return Err(invalid(format!("unsupported Lean expression {token}"))),
        };
        while self.peek() == "." {
            self.take()?;
            self.require("derivative")?;
            expression = Expression::Derivative(Box::new(expression));
        }
        Ok(expression)
    }
}

fn coerce(value: Value, kind: Kind) -> Result<Value, RecurrenceError> {
    match (value, kind) {
        (Value::Natural(n), Kind::Natural) => Ok(Value::Natural(n)),
        (Value::Natural(n), Kind::Real) => Ok(Value::Real(BigRational::from_integer(n))),
        (Value::Natural(n), Kind::Polynomial) => {
            Ok(Value::Polynomial(vec![BigRational::from_integer(n)]))
        }
        (Value::Real(n), Kind::Real) => Ok(Value::Real(n)),
        (Value::Polynomial(p), Kind::Polynomial) => Ok(Value::Polynomial(p)),
        _ => Err(invalid("Lean read-back type mismatch")),
    }
}

fn polynomial(value: Value) -> Result<Vec<BigRational>, RecurrenceError> {
    if let Value::Polynomial(p) = value {
        Ok(p)
    } else {
        Err(invalid("expected Lean polynomial"))
    }
}
fn natural(value: Value) -> Result<BigInt, RecurrenceError> {
    if let Value::Natural(n) = value {
        Ok(n)
    } else {
        Err(invalid("expected Lean natural"))
    }
}
fn real(value: Value) -> Result<BigRational, RecurrenceError> {
    if let Value::Real(n) = value {
        Ok(n)
    } else {
        Err(invalid("expected Lean real"))
    }
}

fn evaluate(
    expression: &Expression,
    kind: Kind,
    index: usize,
    rows: &[Vec<BigRational>],
) -> Result<Value, RecurrenceError> {
    use Expression::*;
    let value = match expression {
        Number(n) => coerce(Value::Natural(n.clone()), kind)?,
        Index => coerce(Value::Natural(BigInt::from(index)), kind)?,
        Variable => Value::Polynomial(vec![BigRational::zero(), BigRational::one()]),
        Row(argument) => {
            let row = natural(evaluate(argument, Kind::Natural, index, rows)?)?
                .to_usize()
                .ok_or_else(|| invalid("Lean row index exceeds machine bounds"))?;
            Value::Polynomial(
                rows.get(row)
                    .cloned()
                    .ok_or_else(|| invalid("Lean reference is not strictly earlier"))?,
            )
        }
        Constant(argument) => {
            Value::Polynomial(vec![real(evaluate(argument, Kind::Real, index, rows)?)?])
        }
        Cast(argument, cast_kind) => coerce(evaluate(argument, *cast_kind, index, rows)?, kind)?,
        Derivative(argument) => Value::Polynomial(vec_poly::derivative(&polynomial(evaluate(
            argument,
            Kind::Polynomial,
            index,
            rows,
        )?)?)),
        Negate(argument) => match evaluate(argument, kind, index, rows)? {
            Value::Real(x) => Value::Real(-x),
            Value::Polynomial(p) => Value::Polynomial(vec_poly::neg(&p)),
            Value::Natural(_) => return Err(invalid("negative natural in Lean read-back")),
        },
        Binary('^', left, right) => {
            let exponent = natural(evaluate(right, Kind::Natural, index, rows)?)?
                .to_u32()
                .filter(|&e| e <= 10000)
                .ok_or_else(|| invalid("Lean exponent exceeds read-back limit"))?;
            match evaluate(left, kind, index, rows)? {
                Value::Natural(x) => Value::Natural(x.pow(exponent)),
                Value::Real(x) => Value::Real(x.pow(exponent as i32)),
                Value::Polynomial(p) => Value::Polynomial(vec_poly::pow(&p, exponent as usize)),
            }
        }
        Binary(operator, left, right) => {
            let left = evaluate(left, kind, index, rows)?;
            let right = evaluate(right, kind, index, rows)?;
            match (left, right) {
                (Value::Natural(a), Value::Natural(b)) => Value::Natural(match operator {
                    '+' => a + b,
                    '-' => (a - b).max(BigInt::zero()),
                    '*' => a * b,
                    '/' => {
                        if b.is_zero() {
                            BigInt::zero()
                        } else {
                            a / b
                        }
                    }
                    _ => return Err(invalid("unsupported natural operator")),
                }),
                (Value::Real(a), Value::Real(b)) => Value::Real(match operator {
                    '+' => a + b,
                    '-' => a - b,
                    '*' => a * b,
                    '/' => {
                        if b.is_zero() {
                            BigRational::zero()
                        } else {
                            a / b
                        }
                    }
                    _ => return Err(invalid("unsupported real operator")),
                }),
                (Value::Polynomial(a), Value::Polynomial(b)) => Value::Polynomial(match operator {
                    '+' => vec_poly::add(&a, &b),
                    '-' => vec_poly::sub(&a, &b),
                    '*' => vec_poly::mul(&a, &b),
                    '/' => {
                        return Err(invalid(
                            "polynomial division is outside the emitted Lean subset",
                        ))
                    }
                    _ => return Err(invalid("unsupported polynomial operator")),
                }),
                _ => return Err(invalid("Lean binary expression type mismatch")),
            }
        }
    };
    coerce(value, kind)
}

/// Parse and independently evaluate an emitted Lean definition and examples.
///
/// The accepted grammar consists of a simple named `ℕ → ℝ[X]` definition,
/// consecutive numeric base cases, one `n + B` case, and exactly `B + 2`
/// examples. Expressions support typed natural/real arithmetic, `C`, `X`,
/// structural references `NAME n` or `NAME (n + j)` with literal `0 < j < B`,
/// and polynomial derivatives. Arithmetic expressions in recursive references
/// are refused even when they evaluate to an earlier index. Real division by zero is
/// zero; natural subtraction truncates and natural division floors. Polynomial
/// division, arbitrary declarations and unparsed suffixes are refused.
/// Every supplied row (numbered from zero) and every example must agree with
/// the parsed definition. Base examples must use `rfl`; the two later examples
/// use `sorry`. This is a finite semantic safeguard, not Lean elaboration.
pub fn verify_lean_recurrence_output(
    text: &str,
    name: &str,
    supplied_rows: &[Vec<BigRational>],
) -> Result<(), RecurrenceError> {
    validate_lean_recurrence_name(name)?;
    // Join physical continuation lines; delimiters are top-level line prefixes.
    let mut statements = Vec::<String>::new();
    let lines: Vec<_> = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with("--"))
        .collect();
    let lines = if lines.first() == Some(&"noncomputable section") {
        if lines.get(1) != Some(&"open Polynomial") || lines.last() != Some(&"end") {
            return Err(invalid("invalid Lean section wrapper"));
        }
        &lines[2..lines.len() - 1]
    } else {
        &lines[..]
    };
    for line in lines {
        let line = line.trim();
        if line.is_empty() || line.starts_with("--") {
            continue;
        }
        if line.starts_with("def ") || line.starts_with('|') || line.starts_with("example ") {
            statements.push(line.into());
        } else {
            let last = statements
                .last_mut()
                .ok_or_else(|| invalid("orphan Lean continuation"))?;
            last.push(' ');
            last.push_str(line);
        }
    }
    let expected_header = format!("def {name} : ℕ → ℝ[X]");
    if statements.first() != Some(&expected_header) {
        return Err(invalid("invalid Lean definition header"));
    }
    let mut bases = Vec::new();
    let mut recursive = None;
    let mut examples = Vec::new();
    for statement in statements.iter().skip(1) {
        if let Some(case) = statement.strip_prefix('|') {
            if !examples.is_empty() {
                return Err(invalid("Lean case after examples"));
            }
            let (pattern, body) = case
                .split_once("=>")
                .ok_or_else(|| invalid("missing Lean case arrow"))?;
            if let Some(base) = pattern.trim().strip_prefix("n + ") {
                if recursive.is_some() {
                    return Err(invalid("duplicate Lean recursive case"));
                }
                let base: usize = base
                    .trim()
                    .parse()
                    .map_err(|_| invalid("invalid Lean recursive boundary"))?;
                if base != bases.len() {
                    return Err(invalid("missing Lean base case"));
                }
                recursive = Some((base, Parser::new(body, name)?.parse()?));
            } else {
                if recursive.is_some() {
                    return Err(invalid("base case after Lean recursive case"));
                }
                let base: usize = pattern
                    .trim()
                    .parse()
                    .map_err(|_| invalid("invalid Lean base pattern"))?;
                if base != bases.len() {
                    return Err(invalid("missing or duplicate Lean base case"));
                }
                bases.push(Parser::new(body, name)?.parse()?);
            }
        } else if let Some(example) = statement.strip_prefix(&format!("example : {name} ")) {
            let (index, rest) = example
                .split_once('=')
                .ok_or_else(|| invalid("invalid Lean example equality"))?;
            let index: usize = index
                .trim()
                .parse()
                .map_err(|_| invalid("invalid Lean example index"))?;
            let (body, proof) = rest
                .split_once(":= by")
                .ok_or_else(|| invalid("invalid Lean example proof"))?;
            if index != examples.len() {
                return Err(invalid("missing or duplicate Lean example"));
            }
            let expected = if index < bases.len() { "rfl" } else { "sorry" };
            if proof.trim() != expected {
                return Err(invalid("unexpected Lean example proof"));
            }
            examples.push(Parser::new(body, name)?.parse()?);
        } else {
            return Err(invalid("unsupported Lean statement"));
        }
    }
    let (base, recursion) = recursive.ok_or_else(|| invalid("missing Lean recursive case"))?;
    if examples.len()
        != base
            .checked_add(2)
            .ok_or_else(|| invalid("Lean example count overflow"))?
    {
        return Err(invalid("wrong number of Lean examples"));
    }
    fn uses_index(expression: &Expression) -> bool {
        match expression {
            Expression::Index => true,
            Expression::Number(_) | Expression::Variable => false,
            Expression::Binary(_, left, right) => uses_index(left) || uses_index(right),
            Expression::Row(argument)
            | Expression::Constant(argument)
            | Expression::Cast(argument, _)
            | Expression::Negate(argument)
            | Expression::Derivative(argument) => uses_index(argument),
        }
    }
    if bases.iter().chain(&examples).any(uses_index) {
        return Err(invalid("unbound n in a Lean base case or example"));
    }
    fn validate_recursive_references(
        expression: &Expression,
        base: usize,
    ) -> Result<(), RecurrenceError> {
        match expression {
            Expression::Row(argument) => {
                let allowed = match argument.as_ref() {
                    Expression::Index => base > 0,
                    Expression::Binary('+', left, right)
                        if matches!(left.as_ref(), Expression::Index) =>
                    {
                        matches!(right.as_ref(), Expression::Number(j) if j > &BigInt::zero() && j < &BigInt::from(base))
                    }
                    _ => false,
                };
                if !allowed {
                    return Err(invalid(
                    "Lean recursive reference must be NAME n or NAME (n + j) with literal 0 < j < B"));
                }
            }
            Expression::Binary(_, left, right) => {
                validate_recursive_references(left, base)?;
                validate_recursive_references(right, base)?;
            }
            Expression::Constant(argument)
            | Expression::Cast(argument, _)
            | Expression::Negate(argument)
            | Expression::Derivative(argument) => validate_recursive_references(argument, base)?,
            _ => (),
        }
        Ok(())
    }
    validate_recursive_references(&recursion, base)?;
    let mut computed = Vec::new();
    for k in 0..supplied_rows.len().max(examples.len()) {
        let (expression, index) = if let Some(expression) = bases.get(k) {
            (expression, 0)
        } else {
            (&recursion, k - base)
        };
        let row = vec_poly::trim(&polynomial(evaluate(
            expression,
            Kind::Polynomial,
            index,
            &computed,
        )?)?);
        if let Some(expected) = supplied_rows.get(k) {
            if row != vec_poly::trim(expected) {
                return Err(invalid(format!(
                    "Lean read-back disagrees at supplied row {k}"
                )));
            }
        }
        if let Some(example) = examples.get(k) {
            let expected = polynomial(evaluate(example, Kind::Polynomial, 0, &[])?)?;
            if row != vec_poly::trim(&expected) {
                return Err(invalid(format!("Lean read-back disagrees at example {k}")));
            }
        }
        computed.push(row);
    }
    Ok(())
}
