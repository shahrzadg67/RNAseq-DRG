#!/usr/bin/env bash
# Publish the ALREADY-BUILT shinylive site (deploy/site_prod) to GitHub Pages
# (gh-pages branch), authored as YOU — no co-author, no external collaborator.
# Run this yourself:
#   bash deploy/push_ghpages.sh
#
# It uses:
#   - deploy/site_prod            (built by shinylive::export)
#   - your PAT at ~/.config/rnaseq_deploy/github_pat   (already present)
#   - PROD_REPO below             (edit if your repo differs)
set -euo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
SITE="$PROJ/deploy/site_prod"

PROD_REPO="${PROD_REPO:-https://github.com/shahrzadg67/RNAseq-DRG.git}"
GITHUB_USER="${GITHUB_USER:-shahrzadg67}"
AUTHOR_NAME="${GIT_AUTHOR_NAME:-Shahrzad Ghazisaeidi}"
AUTHOR_EMAIL="${GIT_AUTHOR_EMAIL:-shahrzad67@gmail.com}"
TOKEN_FILE="${GITHUB_TOKEN_FILE:-$HOME/.config/rnaseq_deploy/github_pat}"

[[ -d "$SITE" ]] || { echo "!! No built site at $SITE — run the shinylive build first."; exit 1; }
GITHUB_TOKEN="${GITHUB_TOKEN:-$(tr -d '[:space:]' < "$TOKEN_FILE")}"
PUSH_URL="https://${GITHUB_USER}:${GITHUB_TOKEN}@${PROD_REPO#https://}"

echo ">> Publishing $SITE to ${PROD_REPO} (gh-pages) as ${AUTHOR_NAME}"
TMP="$(mktemp -d)"
cp -r "$SITE/." "$TMP/"
touch "$TMP/.nojekyll"                       # GitHub Pages: serve as-is, no Jekyll
( cd "$TMP"
  git init -q && git checkout -qb gh-pages
  git add -A
  git -c user.email="$AUTHOR_EMAIL" -c user.name="$AUTHOR_NAME" \
      commit -qm "Publish TUY35595 RNA-seq app $(date -u +%FT%TZ)"
  git push -f "$PUSH_URL" gh-pages 2>&1 | sed -E "s/${GITHUB_TOKEN}/***/g" )
rm -rf "$TMP"
echo ">> Pushed. In the repo Settings → Pages, set Source = gh-pages branch (once)."
echo ">> Your app will be live at: https://${GITHUB_USER}.github.io/$(basename "${PROD_REPO%.git}")/"
