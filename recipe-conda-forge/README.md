# Conda-Forge Recipe

This directory contains the release recipe intended for conda-forge submission.
It is separate from `recipe/`, which is the development and CI recipe that can
build from a Git checkout, branch, or local Git SHA.

For conda-forge, the source must be immutable and reproducible. This recipe uses
the GitHub archive for the release tag and pins it with `sha256`.

## When To Use

Use `recipe/` for local testing and GitHub Actions builds from current source.

Use `recipe-conda-forge/` when preparing a conda-forge submission or updating the
feedstock recipe.

## Release Refresh Steps

1. Update project versions in the source tree.
2. Commit the version changes.
3. Create and push the release tag, for example `2026.05.1`.
4. Download the GitHub tag archive.
5. Compute the archive SHA256.
6. Update `recipe-conda-forge/meta.yaml` with the version, URL, SHA256, and `build/number: 0`.
7. Copy this directory into `staged-recipes/recipes/nrgljubljana/` for first submission, or into the feedstock `recipe/` directory for later updates.

Example commands for the archive hash:

```bash
curl -L -o nrgljubljana-2026.05.1.tar.gz \
  https://github.com/rokzitko/nrgljubljana/archive/refs/tags/2026.05.1.tar.gz
sha256sum nrgljubljana-2026.05.1.tar.gz
```

On macOS, use:

```bash
shasum -a 256 nrgljubljana-2026.05.1.tar.gz
```

## Current Release

The current recipe is prepared for tag `2026.05.1` with:

```yaml
source:
  url: https://github.com/rokzitko/nrgljubljana/archive/refs/tags/{{ version }}.tar.gz
  sha256: 254e136471f69bdcbab48dd7210d59835e0984649aec145151144b59271afa3b
```

Do not move or recreate the `2026.05.1` tag after computing this hash. If the tag
changes, the GitHub archive changes and the SHA256 must be recomputed.

## Local Checks

Render the recipe with conda-forge packages only:

```bash
conda render --override-channels -c conda-forge recipe-conda-forge
```

Build locally without upload:

```bash
conda build --override-channels -c conda-forge --no-anaconda-upload recipe-conda-forge
```

For feedstock/staged-recipes validation, use the conda-forge infrastructure in
the target repository rather than this source repository.

## Maintenance Notes

Keep `recipe-conda-forge/` synchronized with `recipe/` for build flags,
dependencies, tests, and metadata.

Keep `nrgljubljana_enable_mathematica` disabled by default. The package installs
the optional `nrginit` scripts, but package tests must not invoke Mathematica.

Reset `build/number` to `0` for each new upstream release. Increment it only for
recipe-only rebuilds of the same upstream version.
