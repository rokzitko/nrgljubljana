# Installing NRG Ljubljana with Conda

NRG Ljubljana is available from **conda-forge** as the package `nrgljubljana`. Conda installs the compiled NRG executable, auxiliary tools, libraries, and the `nrginit` scripts together with their software dependencies.

The recommended setup is to install NRG Ljubljana into its own Conda environment.

## Option 1: Starting from scratch

If you do not already use Conda, we recommend **Miniforge**, a minimal Conda distribution configured to use conda-forge.

### 1. Install Miniforge

Download the appropriate **Miniforge3** installer for your system from the Miniforge GitHub release page.

Typical installer names are:

```text
Linux x86-64:       Miniforge3-Linux-x86_64.sh
Linux ARM64:        Miniforge3-Linux-aarch64.sh
macOS Apple Silicon: Miniforge3-MacOSX-arm64.sh
macOS Intel:         Miniforge3-MacOSX-x86_64.sh
```

Run the downloaded installer, for example:

```bash
bash Miniforge3-Linux-x86_64.sh
```

or, on an Apple Silicon Mac,

```bash
bash Miniforge3-MacOSX-arm64.sh
```

Accept the default installation location unless you have a reason to change it. When asked whether Conda should initialize your shell, answer **yes**.

Then close and reopen the terminal.

You can check that Conda is available with

```bash
conda --version
```

### 2. Create an NRG Ljubljana environment

Create a dedicated environment:

```bash
conda create -n nrg nrgljubljana
```

Activate it:

```bash
conda activate nrg
```

NRG Ljubljana is now installed.

---

## Option 2: Using an existing Conda installation

If you already have Conda, Miniconda, Miniforge, Mambaforge, Anaconda, etc., there is no need to install another Conda distribution.

To avoid mixing packages from different channels, create a new environment using only conda-forge:

```bash
conda create -n nrg --override-channels -c conda-forge nrgljubljana
```

Then activate it:

```bash
conda activate nrg
```

Alternatively, if your Conda installation is already configured to use conda-forge with strict channel priority, simply use

```bash
conda create -n nrg nrgljubljana
conda activate nrg
```

If you want to configure conda-forge permanently, use

```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```

and then

```bash
conda create -n nrg nrgljubljana
```

---

## Checking the installation

After activating the environment,

```bash
conda activate nrg
```

check the installed package:

```bash
conda list nrgljubljana
```

and verify that the main programs are available:

```bash
command -v nrg
command -v nrginit
```

You can also check the installation prefix with

```bash
echo "$NRGLJUBLJANA_ROOT"
```

This should normally correspond to the active Conda environment.

## Using NRG Ljubljana later

Whenever you open a new terminal, activate the environment before using NRG Ljubljana:

```bash
conda activate nrg
```

When finished, you can leave the environment with

```bash
conda deactivate
```

## Updating NRG Ljubljana

To update to the latest version available from conda-forge:

```bash
conda activate nrg
conda update nrgljubljana
```

For an existing Conda installation that is not otherwise configured for conda-forge, use

```bash
conda update --override-channels -c conda-forge nrgljubljana
```

## Mathematica and `nrginit`

The Conda package includes the `nrginit` scripts. **Wolfram Mathematica itself is not installed by Conda.**

Mathematica is required for the `nrginit` part of the workflow, which generates the initial NRG problem. The numerically intensive C++ NRG executable can be run without Mathematica once the initialization files have been generated.

Thus, on a computing cluster, Mathematica does not need to be installed on every compute node; initialization can be performed on a machine with Mathematica and the resulting files transferred to the compute system.

## Troubleshooting

To see whether NRG Ljubljana is available for your current platform:

```bash
conda search --override-channels -c conda-forge nrgljubljana
```

To see the active Conda environment and platform information:

```bash
conda info
```

If an old Conda installation has complicated or conflicting channel settings, creating the environment with

```bash
conda create -n nrg --override-channels -c conda-forge nrgljubljana
```

is usually the safest approach.

