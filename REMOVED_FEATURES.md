# Removed Features

This document records removals of deprecated or obsolete features, newest first.
Dates use `YYYY-MM-DD` and identify when the feature was removed.

## 2026-10-09: Obsolete Mathematica Chain Backends

- Removed the `tri=orth`, `tri=nambu`, `tri=sc`, and `tri=sc2` selectors and their implementations from the Mathematica initializer. These selectors now produce an unknown-backend error.
- Removed the dependent energy-dependent DMFT pairing interface (`dmftscdelta`, `dg`, and `dgminus`).
- For scalar chain reconstruction, use `tri=old` or `tri=rkpw`. Energy-dependent pairing requires an external chain generator and manual coefficient input; scalar reconstruction is not an equivalent replacement.
- Constant `bcsgap` pairing, `manual_nambu`, `manual_nambu_new`, superconducting symmetries, and the separate `band=nambu` and `chain_gauge=nambu` settings are unaffected. C++ runtime behavior and data formats are unchanged.
