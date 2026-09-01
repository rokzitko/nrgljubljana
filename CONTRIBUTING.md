# Contributing

NRG Ljubljana grows through research use as well as software development. A
contribution can be a published application, a worked calculation, a model, a
numerical comparison, a clearer explanation, or a helpful answer to another
user. You do not need to know C++ or compile the source code to contribute.

Small contributions are welcome. If you are unsure where to begin, describe
the physics problem or proposed improvement in
[GitHub Discussions](https://github.com/rokzitko/nrgljubljana/discussions).

## Ways To Contribute

### Add A Published Application

If a paper explicitly used NRG Ljubljana, please add it to the
[research bibliography](docs/docs/publications.md). You can
[edit the bibliography in your browser](https://github.com/rokzitko/nrgljubljana/edit/master/docs/docs/publications.md)
or open an issue if you would prefer the maintainers to add it.

Please provide:

- the full citation and DOI or another stable link
- one or two sentences saying which model, observable, or physical result used
  NRG Ljubljana
- a link to the part of the paper or supplementary material that identifies
  the software, when this is not clear from the main text

No local build or test is needed for a bibliography contribution.

### Contribute An Example Or Model

A compact worked calculation can be as useful as a new solver feature. Examples
may show how to define an impurity problem, select useful observables, explore a
parameter regime, or recover a known limit. Model contributions may add a new
Hamiltonian or make an existing model easier to use.

Before preparing a large contribution, start a
[Discussion](https://github.com/rokzitko/nrgljubljana/discussions) or comment on
an open issue with the `example` or `physics` label. The maintainers can help
choose the best location and avoid duplicated work.

Please include the parts that are relevant to the calculation:

- the Hamiltonian and the meaning of its parameters
- the symmetry or conserved quantum numbers being used
- the parameter file, units, and normalization conventions
- the observables to calculate and how to interpret the output
- a simple check, such as an exactly known limit, a sum rule, or comparison
  with a published result
- short instructions for reproducing the calculation

Keep examples small enough to understand and run. Generated spectra and raw
output can be large, so include only the compact reference data needed to show
the expected result. If a model requires Mathematica or a particular symmetry,
say so clearly. Maintainers can help identify the appropriate model and test
files.

### Improve Documentation

Corrections to physical conventions, parameter definitions, equations, units,
and explanations are especially valuable. Typographical fixes and reports of
unclear passages are also welcome.

Current documentation is under [`docs/docs/`](docs/docs/). The older `doc/`
tree remains available while material is being moved. For a small correction,
you may use the **Edit** link on a documentation page and let GitHub create the
proposed change; a local source build is not required.

For a larger documentation change, build the current documentation from the
repository root with:

```sh
make -C docs
```

GitHub also checks the documentation when a pull request is opened.

### Report Benchmark Results

Scientific benchmarks help establish where a calculation is reliable. They may
compare NRG Ljubljana with an analytic limit, a published result, another
many-body method, or an earlier software version. Timing and memory comparisons
on different computers are also useful when the numerical setup is stated.

Share an initial result in
[Show and tell](https://github.com/rokzitko/nrgljubljana/discussions/categories/show-and-tell).
Include enough information for another researcher to understand and reproduce
the comparison:

- the physical problem, Hamiltonian, and observable
- units, normalization, and sign conventions
- the NRG Ljubljana version or source revision
- the parameters that control accuracy, such as `Lambda`, z averaging,
  truncation, temperature, and broadening, when relevant
- the reference result and a quantitative measure of agreement
- any convergence check, uncertainty estimate, or known limitation
- a compact script, table, or figure when it helps explain the result

If a benchmark reveals a reproducible discrepancy, it can then be turned into
a focused issue. There is no need to upload a complete calculation directory or
large raw data files.

### Answer A Discussion Question

Experience from your own calculations is valuable. You are welcome to answer a
[Discussion question](https://github.com/rokzitko/nrgljubljana/discussions/categories/q-a),
point to a relevant paper or documentation page, explain a convention, or share
a minimal input that worked for you. You do not need to be a maintainer or have
a complete answer. State any assumptions and software version that matter for
your reply.

### Report A Bug Or Suggest A Feature

Use the [bug report](https://github.com/rokzitko/nrgljubljana/issues/new?template=bug.md)
for reproducible unexpected behavior and the
[feature request](https://github.com/rokzitko/nrgljubljana/issues/new?template=feature.md)
for a proposed capability. Questions about setting up or interpreting a
calculation usually fit better in Discussions.

## Find A Small Task

Several focused tasks are deliberately kept available for new contributors.
These labels describe the expected kind of contribution:

| Label | Meaning |
| --- | --- |
| [`good first issue`](https://github.com/rokzitko/nrgljubljana/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22) | A bounded task suitable for a first contribution. |
| [`documentation`](https://github.com/rokzitko/nrgljubljana/issues?q=is%3Aissue+is%3Aopen+label%3Adocumentation) | Explanations, references, conventions, and corrections. |
| [`example`](https://github.com/rokzitko/nrgljubljana/issues?q=is%3Aissue+is%3Aopen+label%3Aexample) | Worked calculations, model inputs, and reproducible examples. |
| [`physics`](https://github.com/rokzitko/nrgljubljana/issues?q=is%3Aissue+is%3Aopen+label%3Aphysics) | Physical interpretation, validation, and applications. |
| [`help wanted`](https://github.com/rokzitko/nrgljubljana/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22) | A useful task for which outside experience would help. |

Leave a comment if you would like to work on an issue. It is always reasonable
to ask for clarification before starting.

## Sending A Change

For documentation and bibliography edits, the browser-based GitHub editor is
often sufficient. For other changes, work in a fork or branch and open a pull
request. Keep each contribution focused on one question or result.

In the pull request, explain:

- what changed and why it is useful
- the physical assumptions or conventions that a reader needs to know
- how you checked the result
- any remaining limitation or open question

Link the relevant issue, Discussion, or publication. Review is collaborative;
maintainers may suggest a simpler presentation or an additional physical check.

## Appropriate Checks

Different contributions need different checks. Do not install a C++ toolchain
just to improve wording, answer a question, or add a paper.

| Contribution | Useful check |
| --- | --- |
| Published application | Verify the citation, DOI, and description of how the software was used. |
| Documentation | Read the rendered page, check links and equations, and run `make -C docs` when practical. |
| Example or model | Run the smallest representative calculation and compare with a known limit, sum rule, or reference result. |
| Benchmark | Record enough parameters, conventions, and version information for another researcher to repeat it. |
| C++ solver or numerical tool | Build the affected targets and run focused tests before broader checks. |

The full GitHub checks run after a pull request is opened. If a check fails for
reasons unrelated to your contribution, ask a maintainer for help.

## Changing The C++ Solver

The remaining instructions are for changes to the C++ solver, numerical tools,
or build system. They are not prerequisites for the contribution types above.

### Local Build

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
```

For a debug build:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DCMAKE_BUILD_TYPE=Debug
```

### Local Tests

Run the default test suite with:

```sh
ctest --test-dir build --output-on-failure
```

Useful focused runs:

```sh
ctest --test-dir build --output-on-failure -R '^(store|test_clean|test0_clean)$'
ctest --test-dir build --output-on-failure -R '^(adapt|nrgchain)'
ctest --test-dir build --output-on-failure -R '^(test_dmnrg_only|test_fdm_only|test65_algorithms_mats)$'
```

Some tests are only enabled when Mathematica is detected during configuration.
Long-running suites are enabled with `-DTEST_LONG=ON`.

### Optional Core Checks

Sanitizers:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DASAN=ON -DUBSAN=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Static analysis:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DANALYZE_SOURCES=ON
cmake --build build --target nrgljubljana_c --parallel
```

Legacy Sphinx documentation, only when changing `doc/`:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DBuild_Documentation=ON -DSphinx_Only=ON
cmake --build build --target docs_sphinx --parallel
```

Tests and tools are configured to prefer the freshly built `nrgljubljana_c`
library from the build tree. If you change Mathematica-dependent paths, make
sure the relevant `nrginit` and model-based tests are enabled in your local
configuration.
