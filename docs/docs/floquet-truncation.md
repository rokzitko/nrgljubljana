# Floquet Truncation

With `floquet=true`, the runtime keeps the eigensolver's extended-zone energies
and assigns a separate quasienergy-aware cost to each state. The generic NRG
truncation controls operate on that cost.

## Shell Units

Let

\[
s_N=\texttt{step.scale()}.
\]

For a normally rescaled calculation, \(s_N\) is the current shell energy scale.
With `absolute=true`, \(s_N=1\). The frequency used alongside the diagonalized
shell energies is

\[
\Omega_N=\frac{\Omega}{s_N}.
\]

Here \(\Omega\) is the physical input value from `[extra]` in bandwidth units.

## Recalculated Mode Operator

Before selecting retained states, the runtime recalculates the singlet operator
\(\hat M\) in the new shell eigenbasis. For state \((I,i)\), it reads

\[
\mu_{Ii}=\langle\psi_{Ii}|\hat M|\psi_{Ii}\rangle.
\]

The seed `m` blocks define the spectral interval used for these expectation
values. Mode mixing can make \(\mu_{Ii}\) nonintegral.

## Mean Energy And Ranking Criterion

Let \(e_{Ii}\) be the raw extended-zone shell energy. Define

\[
e_{\mathrm{mode},Ii}=\mu_{Ii}\Omega_N,
\]

\[
\bar e_{Ii}=e_{Ii}-e_{\mathrm{mode},Ii},
\]

and

\[
c_{Ii}=\bar e_{Ii}+|e_{\mathrm{mode},Ii}|.
\]

Equivalently,

\[
c_{Ii}=e_{Ii}-\mu_{Ii}\Omega_N+|\mu_{Ii}\Omega_N|.
\]

Multiplication by \(s_N\) gives the same expression in physical energy units:

\[
C_{Ii}=E_{Ii}-\mu_{Ii}\Omega+|\mu_{Ii}\Omega|.
\]

## Replica Dependence

For a replica shift by integer \(\ell\),

\[
(e,\mu)\longmapsto(e+\ell\Omega_N,\mu+\ell).
\]

The mean-energy part is invariant:

\[
(e+\ell\Omega_N)-(\mu+\ell)\Omega_N=e-\mu\Omega_N,
\]

while the cost becomes

\[
c_\ell=\bar e+\Omega_N|\mu+\ell|.
\]

The runtime therefore preserves the extended-zone energy and uses the mode
expectation separately when ranking replicas.

## Ordering And Reference Shifts

Across all invariant sectors, define

\[
e_{\min}=\min_{I,i} e_{Ii},
\qquad
c_{\min}=\min_{I,i} c_{Ii}.
\]

Within each sector, states and eigenvector rows are permuted into ascending
criterion order. The shifted values used after preparation are

\[
\widetilde e_{Ii}=e_{Ii}-e_{\min},
\qquad
\widetilde c_{Ii}=c_{Ii}-c_{\min}.
\]

The raw energy array remains available in extended-zone form; corrected-energy
and criterion arrays carry their respective reference shifts.

## Retained-State Selection

`truncate_prepare(...)` collects \(\widetilde c\) from every invariant sector
and applies the existing global controls:

- with `keepenergy<=0`, `keep` supplies the fixed retained-state count;
- with `keepenergy>0`, the shell cutoff is
  \(c_{\mathrm{cut}}=\texttt{keepenergy}\,\texttt{step.unscale()}\), and the
  existing selection includes one state beyond the counted cutoff subject to
  `keepmin` and `keep`;
- `keepall` selects all computed states on its configured iterations;
- `safeguard` and `safeguardmax` can extend the boundary across a nearby
  criterion cluster;
- on the final CFS/FDM shell, all computed states are retained unless
  `lastalloverride=true`.

The resulting global criterion threshold determines the retained count in each
sector except for that final-shell override. `truncate_perform()` then applies
those counts.

## Worked Criterion Example

Consider three states with

\[
e=(10,4,7),\qquad \mu=(2,-1,0.5),\qquad \Omega_N=2.
\]

Their mode energies, mean energies, and costs are

\[
e_{\mathrm{mode}}=(4,-2,1),
\]

\[
\bar e=(6,6,6),
\]

\[
c=(10,8,7).
\]

Thus \(e_{\min}=4\) and \(c_{\min}=7\). Criterion ordering produces shifted
energy and cost arrays

\[
\widetilde e=(3,0,6),\qquad \widetilde c=(0,1,3).
\]

This is the synthetic case covered by the core unit test.

## Runtime Sequence

For each newly diagonalized NRG shell, the sequence is:

1. establish the Floquet energy reference;
2. resolve small eigenspectrum splittings;
3. split eigenvectors into ancestor blocks without discarding them;
4. recalculate `m` in the shell eigenbasis;
5. calculate and assign the criteria;
6. order states by criterion and apply the energy and criterion shifts;
7. call the generic `truncate_prepare(...)` path;
8. later apply the selected counts through `truncate_perform()`.

Operator recalculation, measurements, state archives, and optional later phases
continue through the normal iteration flow after criterion preparation.

## Diagnostics

With `dumpenergies=true`, `energies.nrg` contains raw energies. `dumpcorr=true`
adds corrected energies, and `dumpcrit=true` adds shifted criterion values. The
seed block retains the criterion present in generated input; iterative blocks
contain the criterion prepared by the runtime.

Log letters `0` through `4` expose progressively more Floquet preparation data,
including the scaled frequency, mode diagonal, per-state criterion, physical
scaling, and prepared eigenspectrum report.

## Related Pages

- [Floquet formalism](floquet-formalism.md)
- [Floquet model construction with nrginit](floquet-nrginit.md)
- [Iteration engine](iteration-engine.md)
- [Output format reference](output-formats.md)
