# Contributing

## Local Build

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
cmake --build build --parallel
```

For a debug build:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/ -DCMAKE_BUILD_TYPE=Debug
```

## Local Tests

Run the default test suite with:

```sh
ctest --test-dir build --output-on-failure
```

Useful focused runs:

```sh
ctest --test-dir build --output-on-failure -R '^(store|test_clean|test0_clean)$'
ctest --test-dir build --output-on-failure -R '^(adapt|nrgchain)$'
ctest --test-dir build --output-on-failure -R '^(test_dmnrg_only|test_fdm_only|test65_algorithms_mats)$'
```

Some tests are only enabled when Mathematica is detected during configuration.
Long-running suites are enabled with `-DTEST_LONG=ON`.

## Optional Developer Checks

Sanitizers:

```sh
cmake -S . -B build -DASAN=ON -DUBSAN=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Static analysis:

```sh
cmake -S . -B build -DANALYZE_SOURCES=ON
cmake --build build --target nrgljubljana_c --parallel
```

Documentation:

```sh
cmake -S . -B build -DBuild_Documentation=ON -DSphinx_Only=ON
cmake --build build --target docs_sphinx --parallel
```

## Notes

- Tests and tools are configured to prefer the freshly built `nrgljubljana_c`
  library from the build tree.
- If you touch Mathematica-dependent paths, make sure the relevant `nrginit`
  and model-based tests are enabled in your local configuration.
