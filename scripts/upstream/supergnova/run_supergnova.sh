#!/usr/bin/env bash
set -euo pipefail

: "${SUPERGNOVA_PY:?Set SUPERGNOVA_PY to the path of SUPERGNOVA/supergnova.py}"
exec python3 "$SUPERGNOVA_PY" "$@"
