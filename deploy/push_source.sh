#!/usr/bin/env bash
# Push the SOURCE + data (code, analysis scripts, app/data/*.rds, README) to
# GitHub main, authored as YOU — no co-author, no external collaborator.
# Big raw data / caches / regenerable image trees are excluded (see .gitignore).
# Run this yourself:
#   bash deploy/push_source.sh
#
# Uses your PAT at ~/.config/rnaseq_deploy/github_pat.  NOTE: the previous PAT was
# expired — regenerate a GitHub token (repo scope) and update that file first:
#   printf '%s' 'ghp_YOURNEWTOKEN' > ~/.config/rnaseq_deploy/github_pat
#   chmod 600 ~/.config/rnaseq_deploy/github_pat
set -euo pipefail
PROJ="/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
PROD_REPO="${PROD_REPO:-https://github.com/shahrzadg67/RNAseq-DRG.git}"
GITHUB_USER="${GITHUB_USER:-shahrzadg67}"
AUTHOR_NAME="${GIT_AUTHOR_NAME:-Shahrzad Ghazisaeidi}"
AUTHOR_EMAIL="${GIT_AUTHOR_EMAIL:-shahrzad67@gmail.com}"
TOKEN_FILE="${GITHUB_TOKEN_FILE:-$HOME/.config/rnaseq_deploy/github_pat}"
BRANCH="${BRANCH:-main}"

GITHUB_TOKEN="${GITHUB_TOKEN:-$(tr -d '[:space:]' < "$TOKEN_FILE")}"
PUSH_URL="https://${GITHUB_USER}:${GITHUB_TOKEN}@${PROD_REPO#https://}"

echo ">> Staging source tree (excluding raw data / caches / image trees)…"
TMP="$(mktemp -d)"
rsync -a \
  --exclude='.git' \
  --exclude='results/' --exclude='fastq/' --exclude='work/' --exclude='refs/' \
  --exclude='Annotations/' --exclude='nf_cache/' --exclude='nfcore_rnaseq_*/' \
  --exclude='.nextflow*' --exclude='singularity_cache/' \
  --exclude='.Rlib/' --exclude='.extralib/' --exclude='logs/' --exclude='dev_data/' \
  --exclude='analysis/cache/' --exclude='DE_results/' --exclude='GSEA_results/' \
  --exclude='plots_png/' --exclude='batch_corrected_results/' --exclude='CHAT_SESSION.md' \
  --exclude='.noxy_jobid' --exclude='*.pid' \
  --exclude='app/www/gsea/' --exclude='app/www/gsea_noxy/' --exclude='app/www/gsea_bc/' \
  --exclude='app/www/pathview/' --exclude='app/www/pathview_noxy/' --exclude='app/www/pathview_bc/' \
  --exclude='app/www/figs/' --exclude='app/www/multiqc/' --exclude='app/data/*_bc.rds' \
  --exclude='deploy/site_staging/' --exclude='deploy/site_prod/' --exclude='deploy/app_deploy/' \
  "$PROJ"/ "$TMP"/

echo ">> Staged size: $(du -sh "$TMP" | cut -f1)"
( cd "$TMP"
  git init -q && git checkout -qb "$BRANCH"
  git add -A
  git -c user.email="$AUTHOR_EMAIL" -c user.name="$AUTHOR_NAME" \
      commit -qm "Update TUY35595 RNA-seq app + analysis $(date -u +%FT%TZ)"
  git push -f "$PUSH_URL" "$BRANCH" 2>&1 | sed -E "s/${GITHUB_TOKEN}/***/g" )
rm -rf "$TMP"
echo ">> Pushed source to ${PROD_REPO} (${BRANCH}), authored as ${AUTHOR_NAME}."
