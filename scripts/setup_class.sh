#!/usr/bin/env bash
# Clone CLASS at the pinned release into ./class_public (gitignored), then
# build classy into the uv environment.
set -euo pipefail
CLASS_TAG="v3.4.0"
CLASS_COMMIT="64bbab707faf4de4779a9e04edd180fef18d98fa"

cd "$(dirname "$0")/.."
if [ ! -d class_public ]; then
    git clone --branch "$CLASS_TAG" --depth 1 https://github.com/lesgourg/class_public.git class_public
fi
actual=$(git -C class_public rev-parse HEAD)
if [ "$actual" != "$CLASS_COMMIT" ]; then
    echo "class_public is at $actual, expected $CLASS_COMMIT ($CLASS_TAG)" >&2
    exit 1
fi
uv sync

# Experiment B: patched CLASS in a separate tree and install location, so the
# stock build (used by default) is never modified. Select with CLASS_BUILD=patched.
rm -rf class_patched
cp -r class_public class_patched
git -C class_patched apply ../class_patches/idm_g_kernel.patch
(cd class_patched && make clean >/dev/null 2>&1 || true; rm -rf build classy.egg-info python/build)
uv pip install --target build/classy_patched ./class_patched --no-deps --reinstall
