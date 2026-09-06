# sym-poly

`sym-poly` is the multi-crate symmetric-function workspace used for checked
examples on `symmetricfunctions.com`.

## Tableau Breadcrumbs

- `core/src/tableau.rs` contains the shared `Tableau` and `SkewTableau` types,
  standard and semistandard generators, reading words, descents, RSK, promotion,
  evacuation, crystals, charge, and key-tableau helpers.
- `core/src/ssaf.rs` contains semi-standard augmented fillings, including atom
  fillings, key fillings, lock fillings, permuted basements, Macdonald filling
  statistics, and Mason's SSYT-to-SSAF map.
- `multipoly/src/flagged_schur.rs` contains flagged Schur and flagged skew
  Schur tableau generators, row-interval flags, and polynomial weight
  enumerators.
- `qsym/src/schur_qsym.rs` contains composition-tableau and immaculate-tableau
  generators for quasisymmetric Schur and dual immaculate functions.
- `sym/src/hook_schur.rs` contains hook-tableau generators for ordinary
  supersymmetric, or hook, Schur functions.
- `qsym/README.md` tracks which Mason/Assaf/Searles-style QSym Schur variants
  are implemented, missing, or need definition-level verification.
- `multipoly/src/key_polynomial.rs` computes key polynomials and tests them
  against the SSAF weight enumerator.
- `multipoly/src/nonsymmetric_macdonald.rs` contains the full `q,t`
  permuted-basement filling formula, plus the `q = 0` operator-side
  Hall-Littlewood specialization.
- `multipoly/src/slide_polynomial.rs`, `multipoly/src/lock_polynomial.rs`, and
  `multipoly/src/kohnert.rs` contain slide/glide polynomials, finite lock
  polynomials, Kohnert/Lascoux weight enumerators, connective-K Grothendieck
  polynomials, and Grothendieck-to-Lascoux expansions.

When adding website examples, prefer a small checked Rust example under the
relevant crate's `examples/` directory and cite it from TeX with a short
`% Related Rust:` comment near the example.
