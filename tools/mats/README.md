# `mats`

`mats` computes a thermal Green's function on a fermionic or bosonic Matsubara-frequency mesh from one or more discrete raw spectra. Multiple spectra are z-averaged before the Green's function is evaluated.

## Definition

Let run $q$ contain discrete frequencies $\epsilon_{qj}$ and weights $w_{qj}$. For $N_z$ runs, `mats` merges exactly equal frequencies and computes the averaged weight

$$
\overline{w}(\epsilon)
= \frac{1}{N_z}
  \sum_{q=1}^{N_z}
  \sum_{j:\,\epsilon_{qj}=\epsilon} w_{qj}.
$$

It then evaluates the discrete spectral representation

$$
G(z) = \sum_{\epsilon} \frac{\overline{w}(\epsilon)}{z-\epsilon}
$$

at $z=i\omega_n$. The fermionic frequencies, used by default, are

$$
\omega_n^{\mathrm{F}} = (2n+1)\pi T,
\qquad n=0,1,\ldots,N_{\mathrm{mats}}-1.
$$

With `-b`, the bosonic frequencies are

$$
\omega_n^{\mathrm{B}} = 2n\pi T,
\qquad n=0,1,\ldots,N_{\mathrm{mats}}-1.
$$

Thus the first fermionic frequency is $\pi T$, while the first bosonic frequency is zero. Both meshes have spacing $2\pi T$.

## Usage

```text
mats [options] <name> <Nz> <T> <nrmats>
```

Arguments:

- `<name>` is the raw spectral-data filename.
- `<Nz>` is the number $N_z$ of spectra to average and must be at least 1.
- `<T>` is the temperature, in the same energy units as the input frequencies, and must be positive. The formulas use $k_{\mathrm{B}}=1$.
- `<nrmats>` is the number $N_{\mathrm{mats}}$ of Matsubara points and must be at least 1.

Options:

- `-h` prints command-line help.
- `-v` prints detailed processing information.
- `-o` reads `<name>` directly when `<Nz>` is 1, instead of reading `1/<name>`.
- `-b` uses bosonic frequencies. The default is fermionic frequencies.
- `-O <file>` writes the result to `<file>`. The default is `spec.dat`.
- `-2` reads three-column records and uses the second column as the weight.
- `-3` reads three-column records and uses the third column as the weight.

Lowercase `-o` controls the single-file input layout; uppercase `-O` controls the output filename.

## Input Format

Input files are headerless binary files containing native C++ `double` values.

Without `-2` or `-3`, every record contains two values:

```text
frequency  weight
```

With `-2` or `-3`, every record contains three values:

```text
frequency  second_weight  third_weight
```

`-2` selects `second_weight`, and `-3` selects `third_weight`. Only one weight column is transformed per invocation. Because the files store native binary `double` values, they must have the byte order and floating-point representation expected by the machine running `mats`.

For z-averaging, the default input layout is

```text
1/<name>
2/<name>
...
<Nz>/<name>
```

For a single spectrum, `-o` instead reads `<name>` from the current working directory.

## Output Format

The output is a text file with one row per Matsubara point and three columns:

```text
omega_n  Re[G(i omega_n)]  Im[G(i omega_n)]
```

The first column contains the real number $\omega_n$, not $i\omega_n$. Values are written with 18 significant digits. An existing output file is replaced.

For a bosonic mesh, $G(0)$ is evaluated directly from the spectral sum. If the input contains nonzero weight at $\epsilon=0$, the denominator vanishes; `mats` does not regularize that pole.

## Examples

Evaluate 100 fermionic points from one file and write the default `spec.dat`:

```sh
mats -o spec.bin 1 0.1 100
```

Evaluate 100 bosonic points and choose the output filename:

```sh
mats -o -b -O susceptibility-mats.dat corr.bin 1 0.1 100
```

Average four spectra from `1/spec.bin` through `4/spec.bin`:

```sh
mats spec.bin 4 0.1 100
```

## Building And Testing

From the repository root:

```sh
cmake --build build --target mats
ctest --test-dir build --output-on-failure -R '^mats[12]$'
```

`mats` was written by Rok Zitko, rok.zitko@ijs.si, in November 2012.
