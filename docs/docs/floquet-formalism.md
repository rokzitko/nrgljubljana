# Floquet Formalism

The current Floquet interface represents one periodic drive in a finite Fourier,
or Sambe, space. This page fixes the notation used by the model-construction and
truncation pages.

## Periodic Hamiltonian

For a Hamiltonian with period \(T\),

\[
H(t+T)=H(t), \qquad \Omega=\frac{2\pi}{T}.
\]

Write its Fourier components as

\[
H(t)=\sum_k H_k e^{ik\Omega t}.
\]

The current initializer takes the positive angular frequency \(\Omega\) from
`[extra] Omega`.

## Floquet States

A Floquet state is written as

\[
|\Psi_\alpha(t)\rangle
  =e^{-i\epsilon_\alpha t}|u_\alpha(t)\rangle,
\qquad
|u_\alpha(t+T)\rangle=|u_\alpha(t)\rangle.
\]

The periodic part obeys

\[
\left[H(t)-i\partial_t\right]|u_\alpha(t)\rangle
  =\epsilon_\alpha|u_\alpha(t)\rangle.
\]

Expanding

\[
|u_\alpha(t)\rangle=\sum_m e^{im\Omega t}|u_{\alpha m}\rangle
\]

gives the time-independent Sambe-space equation

\[
\sum_n\left[H_{m-n}+m\Omega\delta_{mn}\right]
  |u_{\alpha n}\rangle
  =\epsilon_\alpha|u_{\alpha m}\rangle.
\]

The corresponding inner product averages over one period:

\[
\langle\!\langle u|v\rangle\!\rangle
  =\frac{1}{T}\int_0^T dt\,\langle u(t)|v(t)\rangle.
\]

## Quasienergy Replicas

Multiplying the periodic part by \(e^{i\ell\Omega t}\) produces another
representation of the same time-dependent state:

\[
|u_{\alpha,\ell}(t)\rangle=e^{i\ell\Omega t}|u_\alpha(t)\rangle,
\qquad
\epsilon_{\alpha,\ell}=\epsilon_\alpha+\ell\Omega.
\]

NRG Ljubljana retains the resulting extended-zone energies. The runtime uses a
separate criterion to order states for truncation; it does not fold the stored
energies into a selected quasienergy zone.

## Finite Mode Window

The initializer restricts the Fourier index to

\[
m=-N_{\mathrm{cut}},\ldots,N_{\mathrm{cut}},
\]

so the Floquet factor has dimension \(2N_{\mathrm{cut}}+1\). The finite-space
shift is the projected operator

\[
S_k=\sum_{\substack{m=-N_{\mathrm{cut}}\\
                    m+k\in[-N_{\mathrm{cut}},N_{\mathrm{cut}}]}}^{N_{\mathrm{cut}}}
  |m+k\rangle\langle m|.
\]

Terms crossing either edge are absent. The two edges are not connected.

## Notation In NRG Ljubljana

| Symbol or name | Meaning |
| --- | --- |
| \(m\) | Integer Fourier-mode label. |
| \(\hat M=\sum_m m|m\rangle\langle m|\) | Floquet number operator. |
| `m` | Serialized singlet operator containing \(\hat M\). |
| `m^2` | Optional observable containing \(\hat M^2\). |
| \(\mu_{Ii}\) | Expectation value of the recalculated \(\hat M\) in runtime state \((I,i)\). |
| `Ncut` | Non-negative mode cutoff supplied through `options=Ncut=...`. |

After mode mixing, \(\mu_{Ii}\) need not be an integer even though the spectrum
of the seed operator \(\hat M\) consists of integer mode labels.

## Related Pages

- [Floquet model construction with nrginit](floquet-nrginit.md)
- [Floquet truncation](floquet-truncation.md)
- [Parameter reference](parameter-reference.md)
