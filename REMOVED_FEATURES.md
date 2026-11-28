# Removed Features

This document records removals of deprecated or obsolete features, newest first.
Dates use `YYYY-MM-DD` and identify when the feature was removed.

## 2026-11-28: Doublet Transfer-Mask Environment Variables

- Removed the solver environment variables `NRG_QSTZ_TRANSFER_MASK` and `NRG_QSZTZ_TRANSFER_MASK`. The solver no longer reads them.
- Set `doublet_transfer_mask=true` in the `[param]` block instead. It selects the same restriction for `QSTZ` and `QSZTZ` and is ignored by other symmetry types.

## 2026-11-26: Obsolete Numerical Helpers And Floquet Diagnostic Field

- Removed the old trapezoidal integration headers for `adapt` and `nrgchain`. Both tools use `TabulatedDensity` for integration.
- Removed unused eigenvalue, coefficient, occupation-factor, Matsubara, and quadrature helpers.
- Removed the always-NaN `c_shift` field from the verbose Floquet report. Use the `Clw` line in the console output to read the shift applied to the ranking criterion.

## 2026-11-25: Older Restart Checkpoint Formats

Restart requires unitary v2 and density v3 headers. There is no migration:
use a fresh checkpoint directory if a unitary archive has an unsupported
format. With a compatible unitary chain, `nrg` rebuilds incompatible density
files as needed.

## 2026-10-09: Obsolete Mathematica Chain Backends

- Removed the `tri=orth`, `tri=nambu`, `tri=sc`, and `tri=sc2` selectors and their implementations from the Mathematica initializer. These selectors now produce an unknown-backend error.
- Removed the dependent energy-dependent DMFT pairing interface (`dmftscdelta`, `dg`, and `dgminus`).
- For scalar chain reconstruction, use `tri=old` or `tri=rkpw`. Energy-dependent pairing requires an external chain generator and manual coefficient input; scalar reconstruction is not an equivalent replacement.
- Constant `bcsgap` pairing, `manual_nambu`, `manual_nambu_new`, superconducting symmetries, and the separate `band=nambu` and `chain_gauge=nambu` settings are unaffected. C++ runtime behavior and data formats are unchanged.
