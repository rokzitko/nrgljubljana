#!/usr/bin/env bash
set -euxo pipefail

cmake_bool() {
  case "${1:-}" in
    1|[Oo][Nn]|[Tt][Rr][Uu][Ee]|[Yy][Ee][Ss]) printf 'ON' ;;
    *) printf 'OFF' ;;
  esac
}

job_count() {
  local value="${1:-0}"
  local fallback="${CPU_COUNT:-1}"

  case "${value}" in
    ''|0) value="${fallback}" ;;
  esac

  case "${value}" in
    *[!0-9]*|'')
      printf 'Invalid job count: %s\n' "${value}" >&2
      return 1
      ;;
  esac

  if [ "${value}" -lt 1 ]; then
    printf 'Invalid job count: %s\n' "${value}" >&2
    return 1
  fi

  printf '%s' "${value}"
}

build_tests="$(cmake_bool "${nrgljubljana_build_tests:-OFF}")"
test_long="$(cmake_bool "${nrgljubljana_test_long:-OFF}")"
build_jobs="$(job_count "${nrgljubljana_build_jobs:-0}")"
test_jobs="$(job_count "${nrgljubljana_test_jobs:-0}")"

cmake -S . -B build -G Ninja \
  ${CMAKE_ARGS:-} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
  -DBuild_Tests="${build_tests}" \
  -DTEST_LONG="${test_long}" \
  -DBuild_Documentation=OFF \
  -DMPI_C_COMPILER="${PREFIX}/bin/mpicc" \
  -DMPI_CXX_COMPILER="${PREFIX}/bin/mpicxx" \
  -DMPIEXEC_EXECUTABLE="${PREFIX}/bin/mpiexec"

cmake --build build --parallel "${build_jobs}"

if [ "${build_tests}" = "ON" ]; then
  ctest --test-dir build --output-on-failure --parallel "${test_jobs}"
fi

cmake --install build

mkdir -p "${PREFIX}/etc/conda/activate.d" "${PREFIX}/etc/conda/deactivate.d"

cat > "${PREFIX}/etc/conda/activate.d/nrgljubljana.sh" <<'ACTIVATE_EOF'
export NRGLJUBLJANA_ROOT="${CONDA_PREFIX}"

if [ "${NRGLJUBLJANA_CONDA_PATH_BACKUP+x}" != "x" ]; then
  export NRGLJUBLJANA_CONDA_PATH_BACKUP="${PATH:-}"
fi

case ":${PATH:-}:" in
  *":${CONDA_PREFIX}/nrginit:"*) ;;
  *) export PATH="${CONDA_PREFIX}/nrginit${PATH:+:${PATH}}" ;;
esac
ACTIVATE_EOF

cat > "${PREFIX}/etc/conda/deactivate.d/nrgljubljana.sh" <<'DEACTIVATE_EOF'
if [ "${NRGLJUBLJANA_CONDA_PATH_BACKUP+x}" = "x" ]; then
  export PATH="${NRGLJUBLJANA_CONDA_PATH_BACKUP}"
  unset NRGLJUBLJANA_CONDA_PATH_BACKUP
fi

unset NRGLJUBLJANA_ROOT
DEACTIVATE_EOF
