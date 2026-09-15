# `mixchain`

`mixchain` maps a **matrix** hybridization function onto a Wilson chain with matrix coefficients. `adapt` and
`nrgchain` treat a scalar hybridization function; `mixchain` treats the $N\times N$ Hermitian positive semidefinite
case, in which the channels mix. Each eigenvalue branch of $\Gamma(\omega)$ is discretized logarithmically into a star
Hamiltonian, and block Lanczos maps the star onto a chain. The method follows Liu et al. (2016).

It works in two stages with a file between them: the star stage reads the components of $\Gamma$ and writes
`star.dat`; the chain stage reads `star.dat` and writes `chain.dat`. The chain is written in its general matrix form;
writers for the input files of particular symmetry types of `nrg` are not yet provided.

## Usage

```text
mixchain [options] [s|l] [parameter_file]
```

- `parameter_file` defaults to `param`.
- `s` runs the star stage: it reads $\Gamma$ and writes `star.dat`.
- `l` runs the chain stage: it reads `star.dat` and writes `chain.dat`.
- With neither, both stages run in turn. The chain is built from `star.dat` in this case too, so the default mode and `s` followed by `l` give the same chain.
- `--Nz N` runs for $z_i = i/N$, $i = 1,\ldots,N$, with the files for $z_i$ in the directory `i/`. A `z` in the parameter file is then ignored. The star stage does everything that does not depend on $z$ once and shares it between the values of $z$; `l` reads `i/star.dat` and writes `i/chain.dat`, and checks that each star was built for $z_i$.
- `--epsabs VALUE`, `--epsrel VALUE`, `--workspace-limit N` and `--gsl-error-policy ignore|warn|fail` control the CQUAD integration of the representative energies, with the names, defaults and meaning they have for `adapt --integral`.
- `-v` writes the resolved configuration to standard error.
- `-vv` is accepted for uniformity with the other tools and currently reports the same as `-v`.
- `-V` or `--version` prints the project version and exits immediately.
- `-h` or `--help` prints the command synopsis.

Options and positional arguments may be given in either order. Each stage prints its diagnostics and its wall-clock
time to standard output, with the shared setup of the star stage and each value of $z$ timed separately. Input and
output paths are relative to the working directory.

```sh
mixchain                 # star and chain
mixchain s custom.param  # the star only
mixchain -v l            # the chain from an existing star.dat
mixchain --Nz 4          # stars and chains for z = 1/4, 1/2, 3/4, 1 in 1/ .. 4/
```

## Input

One two-column file $(\omega, \text{value})$ per component of $\Gamma$, named `<dos_prefix>_ij-re.dat` and
`<dos_prefix>_ij-im.dat` for $i,j = 1,\ldots,$ `channels`. Each file covers both negative and positive frequencies,
as `Delta.dat` does for the scalar tools.

- All `-re` files are required, including the diagonal.
- The off-diagonal `-im` files are required either all or none. None means that $\Gamma$ is real, and the whole calculation then runs in real arithmetic. A diagonal `-im` file is optional and must vanish.
- Every file must be tabulated on one and the same frequency grid; the input is never reinterpolated.
- $\Gamma_{ji}$ is checked against $\Gamma_{ij}^*$ relative to the largest element of $\Gamma$, with tolerance `hermiticity_tolerance`, and $\Gamma$ is then symmetrized to $(\Gamma+\Gamma^\dagger)/2$. Positive semidefiniteness is checked at every node.

**Normalization.** As for the `dos` file of `adapt` and `nrgchain`, the input is $\pi$ times the spectral function
of the hybridization function. With

$$
\Delta(z) = \sum_k \frac{V_k V_k^\dagger}{z-\epsilon_k}, \qquad
\Gamma(\omega) = -\frac{1}{2\pi i}\left[\Delta(\omega+i0^+) - \Delta(\omega+i0^+)^\dagger\right]
= \sum_k V_k V_k^\dagger\,\delta(\omega-\epsilon_k),
$$

the files hold $\Gamma_{\rm in}(\omega) = \pi\,\Gamma(\omega)$, which for $N=1$ is $-{\rm Im}\,\Delta(\omega+i0^+)$.
An overall factor in $\Gamma$ leaves every chain coefficient unchanged and scales only the hybridization weight
$\Theta = \int\Gamma_{\rm in}\,d\omega = \pi\sum_k V_kV_k^\dagger$, so the convention matters only for the impurity
coupling.

**Band.** `bandrescale` maps the band edge to 1: $\omega\to\omega/$`bandrescale` and $\Gamma\to\Gamma\cdot$`bandrescale`,
which leaves $\int\Gamma\,d\omega$ unchanged. The mesh reaches only $|\omega|\le1$, so weight tabulated beyond the band
edge is discarded, and where the input stops short of the edge the density is continued at its last tabulated value.
Both are reported per diagonal element when they occur.

## Parameters

Read from the `[param]` block of the parameter file, with the names of `adapt` and `nrgchain` wherever the quantity
is the same.

| Parameter | Default | Stage | Meaning |
| --- | --- | --- | --- |
| `channels` | `1` | star | Dimension $N$ of $\Gamma$. |
| `dos_prefix` | `Gamma` | star | Prefix of the input files. |
| `Lambda` | `2` | star | Discretization parameter $\Lambda$. |
| `z` | `1` | star | Twist parameter. |
| `mMAX` | `2*Nmax` | star | Largest interval index; the star has `2*channels*(mMAX+1)` levels. |
| `bandrescale` | `1` | star | Band rescaling. |
| `adapt` | `false` | star | Adaptive mesh, see below. |
| `mesh_weight` | `frobenius` | star | Weight function of the adaptive mesh: `frobenius` or `trace`. |
| `hardgap`, `boundary` | `false`, `0` | star | Accumulation point of the mesh, as a fraction of the rescaled band edge. |
| `density_interpolation` | `linear` | star | `linear` or `steffen`, as for `adapt` and `nrgchain`. |
| `branch_ordering` | `tracked` | star | `tracked` or `sorted`. |
| `allowed_error` | `1e-10` | star | Default relative tolerance of the integral method. |
| `hermiticity_tolerance` | `1e-8` | star | Allowed deviation of the input from a Hermitian matrix. |
| `Nmax` | required | chain | Last site of the chain, which has the sites `0..Nmax`. |
| `preccpp` | `2000` | chain | Precision of the chain stage in bits, as for `nrgchain`. |
| `breakdown_tolerance` | `1e-20` | chain | Relative size below which a matrix to be inverted counts as singular. |

`boundary` is a fraction of the rescaled band edge, as in `adapt`: a gap $\Delta$ in the units of the input with
`bandrescale`$=D$ is `boundary`$=\Delta/D$. The `-v` report prints both values.

The chain stage takes $\Lambda$, $z$ and `bandrescale` from `star.dat`. If the parameter file sets $\Lambda$ or
`bandrescale` to a different value, the stage stops rather than choose between the two. The same holds for $z$: with
`--Nz` it must be $i/N$ for the star in `i/`, and otherwise it must match a `z` given in the parameter file.

## Star stage

### Mesh

The interval of index $m = 0,\ldots,$`mMAX` is $[\epsilon(z+m+2),\,\epsilon(z+m+1)]$, where

$$
\epsilon(x) = \Lambda^{2-x} \quad\text{for } \texttt{adapt=false}, \qquad
\epsilon(x) = W_{\rm mesh}^{-1}\!\left(\Lambda^{2-x}\right) \quad\text{for } \texttt{adapt=true},
$$

for $x>2$, and $\epsilon(x)=1$ otherwise. $W_{\rm mesh}$ is the cumulative of a scalar weight function of $\Gamma$,
normalized to its value at $\omega=1$: $\lVert\Gamma(\omega)\rVert_F$ for `mesh_weight=frobenius`, or
${\rm tr}\,\Gamma(\omega)$ for `mesh_weight=trace`. Both forms are followed by the `hardgap` rescaling
$\epsilon\to(1-b)\,\epsilon+b$, with $b$ = `boundary`.

This adaptive mesh is **not** the adaptive mesh of `adapt`, which solves a differential equation for a guiding
function $g(x)$ and writes it to `GSOL.dat`. `mixchain` reads no such file and solves no differential equation: its
intervals follow the structure of $\Gamma$ directly. Where the weight function vanishes identically below some
frequency, as it does inside a gap, $W_{\rm mesh}$ is flat there and the mesh accumulates at the edge of that region
without `hardgap`.

### Branches

$\Gamma$ is diagonalized at the input nodes, $\Gamma(\omega) = \sum_a \rho_a(\omega)\,u_a(\omega)\,u_a(\omega)^\dagger$.
With `branch_ordering=tracked`, the eigenvectors at neighbouring nodes are matched by overlap going outward from
$\omega=0$, so that a branch keeps its identity where two eigenvalues cross, and degenerate subspaces are rotated to
match the previous node. With `branch_ordering=sorted`, the branches are labelled by descending eigenvalue at every
node, which exposes the crossing artefact: across a crossing the weight of one channel is attached to the eigenvector
of another.

### Star levels

Each branch $\rho_a$ is treated as a scalar density: its interval weight $w_a$ is the integral of the interpolant, and
its representative energy is given by the integral method of `adapt`,

$$
E_a(x) = W_a^{-1}\!\left[\int_x^{x+1} W_a(\epsilon(x'))\,dx'\right],
$$

with $W_a$ the normalized cumulative of $\rho_a$. The eigenvector $u_a$ at $E_a$ comes from $\Gamma$ interpolated
element by element and diagonalized there. Each interval, branch and frequency branch then gives one bath level with
energy $\pm E_a$ and coupling vector $v = \sqrt{w_a}\,u_a$, in the normalization of the input:

$$
H = \sum_k E_k\,c_k^\dagger c_k + \sum_k\sum_i \left(v_{k,i}\,d_i^\dagger c_k + {\rm h.c.}\right).
$$

A branch without weight in some interval gives a level with vanishing coupling, which is kept.

## Chain stage

Block Lanczos maps the star onto

$$
H = \sum_{ij}\left(V_{ij}\,d_i^\dagger f_{0j} + {\rm h.c.}\right)
+ \sum_n\sum_{ij}(E_n)_{ij}\,f_{ni}^\dagger f_{nj}
+ \sum_n\sum_{ij}\left((T_n)_{ij}\,f_{n+1,i}^\dagger f_{nj} + {\rm h.c.}\right),
$$

with $N\times N$ blocks, in the polar gauge. Lanczos fixes each site only up to a unitary rotation of its $N$
orbitals; the polar gauge removes that freedom by taking $V$ and every $T_n$ Hermitian positive semidefinite, the
matrix analogue of choosing $\xi_n>0$. Writing $H_{\rm bath}$ for the diagonal of the star energies and $A$ for
the matrix with $A_{ki} = v_{k,i}^*$, so that $A^\dagger A = \Theta$,

$$
V = \Theta^{1/2}, \qquad Q_0 = A\,\Theta^{-1/2},
$$

$$
E_n = Q_n^\dagger H_{\rm bath} Q_n, \qquad
R = H_{\rm bath}Q_n - Q_nE_n - Q_{n-1}T_{n-1}^\dagger, \qquad
T_n = (R^\dagger R)^{1/2}, \qquad Q_{n+1} = R\,(R^\dagger R)^{-1/2}.
$$

The residual is reorthogonalized against every earlier block. In the polar gauge $V$ and every $T_n$ are Hermitian
positive semidefinite, the gauge involves no preferred basis, and a rotation of the channels rotates every block in
the same way.

**Precision.** The late coefficients fall off as $\Lambda^{-n/2}$, below the resolution of double precision, so the
recursion runs in multiprecision arithmetic. Its precision is chosen at compile time from a ladder of 50, 200 and 800
decimal digits: `preccpp` bits select the smallest rung that covers them, and the `-v` report shows the result. The
default of 2000 bits resolves to 800 digits; requests beyond 800 digits are rejected. The result is written with 18
significant digits, as `nrgchain` writes `xi.dat`.

**Breakdown.** Inverting $\Theta$ fails when $\Gamma$ is rank deficient over the whole band, so that the bath couples
to fewer combinations of the impurity orbitals than there are channels. Inverting $R^\dagger R$ fails part-way down
the chain when the Krylov space of the star is exhausted, for instance when too few levels carry weight. Both stop the
stage with the site and the rank. The star must have at least `channels*(Nmax+1)` levels.

## Outputs

Both files are written to the working directory, or to `i/` for $z_i$ with `--Nz`.

### `star.dat`

One row per bath level. Lines beginning with `#` are comments; the one carrying `channels=` is the header, read back
by the chain stage as whitespace-separated `key=value` pairs.

| Header key | Meaning |
| --- | --- |
| `channels` | Dimension $N$, and the number of components of every coupling vector. |
| `mMAX` | Largest interval index; the file holds `2*channels*(mMAX+1)` rows. |
| `z`, `Lambda` | The discretization. |
| `bandrescale` | The rescaling applied to $\Gamma$; energies are in the rescaled band. |
| `complex` | `1` if the couplings are complex, `0` if real. It fixes the number of columns. |

| Column | Meaning |
| --- | --- |
| `m` | Interval index; a larger `m` lies closer to the Fermi level. |
| `sign` | `+` for the positive, `-` for the negative frequency branch. |
| `a` | Eigenvalue branch of $\Gamma$, `0..channels-1`. |
| `E` | Representative energy, with the sign of its frequency branch. |
| `v1 ... vN` | Coupling vector in the channel basis; for `complex=1` each component is a pair `Re_vi Im_vi`. |

`a` is a label from the branch tracking, not a channel: a branch is an eigenvector direction of $\Gamma$, which in
general points across several channels and rotates with $\omega$. `m`, `sign` and `a` are not used by the chain
stage, which sees only the pairs $(E, v)$; they are written so that the file can be read and checked. The star
diagnostics are written as comments and are not read back.

### `chain.dat`

One row per matrix element. The second line is the header; the third holds the diagnostics of the recursion.

| Header key | Meaning |
| --- | --- |
| `channels` | Dimension of every block. |
| `Nmax` | Last site. |
| `z`, `Lambda`, `bandrescale` | As in `star.dat`. |
| `complex` | `1` if the coefficients are complex, `0` if real. |
| `digits` | Decimal digits of the arithmetic the recursion ran in. |

| Column | Meaning |
| --- | --- |
| `block` | `V`, `E` or `T`. |
| `n` | Site: `0` for `V`, `0..Nmax` for `E`, `0..Nmax-1` for `T`. |
| `i`, `j` | Matrix indices, `1..channels`, as in `Gamma_ij`. |
| `value` | The element; for `complex=1` a pair `Re Im`. |

$V_{ij}$ multiplies $d_i^\dagger f_{0j}$, $(E_n)_{ij}$ multiplies $f_{ni}^\dagger f_{nj}$, and $(T_n)_{ij}$
multiplies $f_{n+1,i}^\dagger f_{nj}$. $V$ is in the normalization of the input, $V^2 = \Theta$; the physical
coupling is $V/\sqrt{\pi}$.

## Diagnostics

The star stage reports:

- per interval, the matrix sum rule $\lVert\sum_a w_a u_a u_a^\dagger - \int\Gamma\rVert / \lVert\int\Gamma\rVert$, with the largest value and where it occurred. It detects a mislabelled branch as well as eigenvectors that rotate too fast for the interval to resolve;
- $\Theta = \sum_k v_k v_k^\dagger$ against $\int\Gamma$ over the covered range. The trace agrees by construction, the off-diagonal elements only approximately;
- the frequencies at which the tracked and the sorted branch orderings diverge;
- weight beyond the band edge that is discarded, and weight added by continuing the input to the edge, per diagonal element;
- intervals that contain no tabulated point of the input, and whether the input ends there or is merely too coarse. In such intervals the star follows the interpolant rather than the data;
- levels lost to double precision near an accumulation point away from zero. Where the distance to that point drops below the spacing of doubles, the bounds of an interval coincide and its levels carry no weight, which in effect truncates the star there;
- the largest CQUAD error estimate of the integral method.

The chain stage reports the ratio of the smallest to the largest eigenvalue of $\Theta$ and, over the chain, of
$R^\dagger R$, together with the largest anti-Hermitian part removed from an on-site block and the largest component
removed by reorthogonalization. The last two measure rounding and loss of orthogonality at the working precision.

## Flat-band benchmark

For a flat band at $z=1$ the star has $E_m = c\,\Lambda^{-m}$ with $c = (1-\Lambda^{-1})/\ln\Lambda$ and weights
$w_m\propto\Lambda^{-m}(1-\Lambda^{-1})$. Wilson's discretization has the same weights and $c_W = (1+\Lambda^{-1})/2$,
so the chain is Wilson's closed form divided by the Campo–Oliveira factor $A_\Lambda = c_W/c$:

$$
\xi_n = \frac{1-\Lambda^{-1}}{\ln\Lambda}\,
\frac{(1-\Lambda^{-n-1})\,\Lambda^{-n/2}}{\sqrt{(1-\Lambda^{-2n-1})(1-\Lambda^{-2n-3})}},
\qquad \zeta_n = 0.
$$

For $\Lambda=2$, `mMAX=80` and $n<20$, `mixchain` reproduces $\xi_n$ to machine precision; $\Gamma=\rho\,\mathbb{1}$
gives $T_n=\xi_n\mathbb{1}$. With the usual `mMAX=2*Nmax`, the truncation of the star shows at the end of the chain:
for `Nmax=20` the deviation starts at $2\times10^{-13}$ and doubles every two sites, reaching $1.7\times10^{-10}$ at
the last site. The unit tests are in `test/unit/mixchain`.

## References

- J.-G. Liu et al., *Physical Review B* **93**, 035102 (2016).
- K. G. Wilson, "The renormalization group: Critical phenomena and the Kondo problem", *Reviews of Modern Physics* **47**, 773 (1975).
- V. L. Campo and L. N. Oliveira, "Alternative discretization in the numerical renormalization-group method", *Physical Review B* **72**, 104432 (2005).
- Rok Zitko, "Adaptive logarithmic discretization for numerical renormalization group methods", *Computer Physics Communications* **180**, 1271-1276 (2009).
- Rok Zitko and Thomas Pruschke, "Energy resolution and discretization artefacts in the numerical renormalization group", *Physical Review B* **79**, 085106 (2009).
