#!/usr/bin/env bash
# Build the shinylive site for a target and (optionally) push PROD to GitHub.
#   bash deploy/publish.sh staging        # full build, all tabs -> local preview only
#   bash deploy/publish.sh prod           # client build, only released tabs (local)
#   bash deploy/publish.sh prod --push    # also commit+push prod to the public repo's gh-pages
#
# Staging is LOCAL ONLY (free GitHub plan): preview via the port-forward command printed below.
# Only PROD is published to GitHub.
#
# Gating: edit config/stage_status.yaml first. "Release stage N" = set its
# status: in_progress + released: true, then `bash deploy/publish.sh prod --push`.
#
# Commits are authored as YOU (no co-author trailer). Configure via env:
#   GIT_AUTHOR_NAME / GIT_AUTHOR_EMAIL  (default: Shahrzad's name + email)
#   GITHUB_USER  + GITHUB_TOKEN (a PAT)  -> used for HTTPS push, never printed
#   PROD_REPO=https://github.com/<USER>/tuy35595-rnaseq.git
set -euo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
TARGET="${1:-staging}"
PUSH="${2:-}"
cd "$PROJ"

export R_LIBS_USER="$PROJ/.Rlib"
# libarchive provides the runtime lib for shinylive's `archive` dependency.
module load R/4.5.1 libarchive/3.8.1 2>/dev/null || true

AUTHOR_NAME="${GIT_AUTHOR_NAME:-Shahrzad Ghazisaeidi}"
AUTHOR_EMAIL="${GIT_AUTHOR_EMAIL:-shahrzad67@gmail.com}"
# Public client repo (override with env if needed).
PROD_REPO="${PROD_REPO:-https://github.com/shahrzadg67/RNAseq-DRG.git}"
GITHUB_USER="${GITHUB_USER:-shahrzadg67}"

SITE="$PROJ/deploy/site_${TARGET}"
echo ">> Generating build_config.json for target=$TARGET"
Rscript deploy/make_build_config.R "$TARGET"

echo ">> Exporting shinylive site -> $SITE"
Rscript -e "shinylive::export('app', '${SITE}')"

echo ">> Built: $SITE"
echo "   Local preview:  Rscript -e \"httpuv::runStaticServer('${SITE}', port=8008)\"  (then forward port 8008 in VS Code)"

if [[ "$PUSH" == "--push" ]]; then
  [[ "$TARGET" == "prod" ]] || { echo "!! Only prod is pushed to GitHub. Staging is local-only."; exit 1; }
  [[ -n "$PROD_REPO" ]] || { echo "!! Set PROD_REPO (https://github.com/<USER>/tuy35595-rnaseq.git)"; exit 1; }
  # Load the PAT from the secure file if not already in the environment.
  TOKEN_FILE="${GITHUB_TOKEN_FILE:-/home/sghazis/.config/rnaseq_deploy/github_pat}"
  if [[ -z "${GITHUB_TOKEN:-}" && -f "$TOKEN_FILE" ]]; then
    GITHUB_TOKEN="$(tr -d '[:space:]' < "$TOKEN_FILE")"
  fi
  # Build an authenticated push URL from GITHUB_USER/GITHUB_TOKEN without printing it.
  PUSH_URL="$PROD_REPO"
  if [[ -n "${GITHUB_TOKEN:-}" && -n "${GITHUB_USER:-}" ]]; then
    PUSH_URL="https://${GITHUB_USER}:${GITHUB_TOKEN}@${PROD_REPO#https://}"
  fi
  echo ">> Publishing prod site to ${PROD_REPO} (gh-pages) as ${AUTHOR_NAME}"
  TMP="$(mktemp -d)"
  cp -r "$SITE/." "$TMP/"
  touch "$TMP/.nojekyll"
  ( cd "$TMP"
    git init -q && git checkout -qb gh-pages
    git add -A
    git -c user.email="$AUTHOR_EMAIL" -c user.name="$AUTHOR_NAME" \
        commit -qm "Publish TUY35595 report $(date -u +%FT%TZ)"
    git push -f "$PUSH_URL" gh-pages 2>&1 | sed -E "s/${GITHUB_TOKEN:-__none__}/***/g" )
  rm -rf "$TMP"
  echo ">> Pushed. Enable GitHub Pages (gh-pages branch) once in repo Settings."
fi
