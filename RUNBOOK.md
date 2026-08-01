# nf-core/rnaseq Runbook — TUY35595 (mouse)

Complete reference for the RNA-seq run set up on the SickKids HPC. Everything lives
under `/hpf/projects/msalter/sghazis/rnaseq_TUY35595/` — nothing touches the 10 GB home.

## Goal
Gene-level **and** transcript-level count matrices from 14 paired-end mouse samples,
in a single run. `--aligner star_salmon` makes Salmon quantify both at once.

## Study / inputs
- **Organism:** mouse, Ensembl **GRCm39**, release **116**.
- **Raw data:** `/hpf/projects/msalter/TUY35595.20260429/20260424_LH00403_0176_B23TKL3LT4`
  — 14 paired-end samples (28 `.fastq.ora`), single lane (L001), 151 bp reads.
- **Samples:** `NP_{M,F}1..3`, `Sham_{M,F}1..3`, `SNI_M`, `SNI_F`
  (pain/nerve-injury model: NP = nerve pain, Sham = control, SNI = spared/spinal nerve injury).
- **Key gotcha:** reads are **ORA-compressed** — nf-core can't read them, so they're
  decompressed to `.fastq.gz` first with the `orad/2.7.0` module.

## Environment facts (this cluster)
- Writable lab base: **`/hpf/projects/msalter`** (NOT `/hpf/largeprojects`, which doesn't exist).
- Home = 10 GB → all caches redirected here via `NXF_HOME`, `NXF_SINGULARITY_CACHEDIR`, `NXF_WORK`.
- No persistent personal scratch; per-job node-local `/localhd/$SLURM_JOBID` only.
- Slurm account = `msalter`; default partition `general`.
- **Login AND compute nodes have outbound internet** (Ensembl, quay.io, GitHub, S3) — so
  references/containers were staged directly, no data-mover SSH hop needed.
- Modules: `orad/2.7.0`, `nf-core/4.2.0` (→ nextflow 24.11.0 + Singularity 3.11.3 + python).
- nf-core/rnaseq pinned to **3.19.0**.

---

## Directory layout
```
rnaseq_TUY35595/
├── RUNBOOK.md / README.md       # this file
├── samplesheet.csv              # 14 samples, strandedness=auto
├── sickkids_slurm.config        # Slurm executor + Singularity + resource caps
├── run_rnaseq.sh                # REAL run: sbatch controller (survives logout)
├── refs/                        # GRCm39 fasta + GTF (rel 116) [+ saved indices after run]
├── fastq/                       # 28 decompressed .fastq.gz
├── nf_cache/                    # NXF_HOME
├── singularity_cache/           # NXF_SINGULARITY_CACHEDIR (~50 container images)
├── work/                        # Nextflow work dir (DELETE after success)
├── results/                     # final outputs (KEEP)
├── logs/                        # all job + script logs
└── scripts/
    ├── 01_ora_to_fastq.sh       # ORA -> fastq.gz (Slurm array, 28 tasks)
    ├── 02_get_reference.sh      # download GRCm39 fasta + GTF from Ensembl
    ├── 03_prestage_pipeline.sh  # pre-pull pipeline + Singularity containers
    ├── 04_make_samplesheet.sh   # (re)generate samplesheet.csv
    ├── 05_smoke_test.sh         # -profile test validation run (Slurm)
    ├── 06_cleanup_after_success.sh  # delete work/ once both matrices exist
    └── 07_autochain.sh          # detached launcher: wait staging -> smoke -> real
```

---

## Manual run order (if starting fresh)
```bash
cd /hpf/projects/msalter/sghazis/rnaseq_TUY35595

sbatch scripts/01_ora_to_fastq.sh        # 1. decompress ORA -> fastq.gz (~51GB -> ~180GB)
bash   scripts/02_get_reference.sh        # 2. download GRCm39 fasta + GTF (rel 116)
bash   scripts/03_prestage_pipeline.sh    # 3. pre-pull pipeline + containers (~50 images)
# scripts/04 already produced samplesheet.csv; rerun only if fastq names change
sbatch scripts/05_smoke_test.sh           # 5. validate Slurm + Singularity + caches
sbatch run_rnaseq.sh                       # 6. real run (or chain via afterok, see below)
bash   scripts/06_cleanup_after_success.sh # 7. reclaim space after success
```

## Self-driving / disconnect-proof launch (what was actually used)
A detached launcher chains the rest so it runs with nobody connected:
```bash
cd /hpf/projects/msalter/sghazis/rnaseq_TUY35595
setsid bash scripts/07_autochain.sh >logs/autochain.log 2>&1 </dev/null &
```
It waits for container staging to finish, then submits the smoke test, then submits the
real run with `--dependency=afterok:<smoke>` (real fires only if smoke succeeds).

### What survives a disconnect
- **Slurm jobs** (decompression array, smoke test, real run): survive everything.
- **Login-node background processes** (container staging, the `setsid` launcher):
  survive a normal SSH/VS Code disconnect; if killed, just re-run the relevant script.
- The **real run is an `sbatch` job**, so once submitted it runs to completion regardless.

---

## Monitoring
```bash
cd /hpf/projects/msalter/sghazis/rnaseq_TUY35595
squeue -u sghazis                     # running/queued Slurm jobs
cat  logs/autochain.log               # launcher progress + submitted job IDs
cat  logs/autochain_jobids.txt        # smoke= / real= job numbers (exists once chain fired)
tail -f logs/rnaseq_ctl_*.out         # live pipeline progress (real run)
tail -f logs/rnaseq_test_*.out        # smoke test progress
du -sh singularity_cache              # container staging size (~50 images when done)
df -h /hpf/projects/msalter           # shared free space (started ~1.4 TB)
```
Stage cues: no Slurm pipeline job + staging still downloading = pre-stage phase (normal).
An `rnaseq_test` job in `squeue` = staging done, chain fired. An `rnaseq_ctl` job = real run.

---

## Outputs (the deliverable) — `results/star_salmon/`
- **Gene level:** `salmon.merged.gene_counts.tsv`, `salmon.merged.gene_tpm.tsv`
  (+ `salmon.merged.gene_counts_length_scaled.tsv`, `.rds` SummarizedExperiment).
- **Transcript level:** `salmon.merged.transcript_counts.tsv`,
  `salmon.merged.transcript_tpm.tsv` (+ `.rds`).
- Per-sample Salmon dirs; QC at `results/multiqc/`; run reports in `results/pipeline_info/`.
- Built STAR + Salmon indices saved under `results/genome/` (from `--save_reference`)
  so re-runs skip rebuilding.

## After a successful run — reclaim space (minimize-footprint choice)
```bash
bash scripts/06_cleanup_after_success.sh   # refuses unless both matrices exist; deletes work/
```
Keeps: `results/`, `fastq/`, `refs/`. Removes: `work/`, `work_test/`, `results_test/`.

---

## Troubleshooting
- **Real run stuck `DependencyNeverSatisfied`** → the smoke test failed. Check
  `logs/rnaseq_test_*.err`. If it's a QOS/account rejection on `general`, add
  `--qos=<normal_qos>` to `clusterOptions` in `sickkids_slurm.config`, then resubmit
  the real run manually: `sbatch run_rnaseq.sh`.
- **A pipeline step fails mid-run** → fix the cause, then resubmit `sbatch run_rnaseq.sh`;
  it uses `-resume` and continues from cached steps.
- **Container won't pull on a node** → containers are pre-staged in `singularity_cache`;
  compute nodes also have internet as fallback. Re-run `scripts/03_prestage_pipeline.sh`.
- **Home filling up** → shouldn't happen (all caches redirected); verify with `du -sh ~`.
- **Resume after editing config** → `sbatch run_rnaseq.sh` (the `-resume` flag is built in).

## Resume from scratch on a new login session
All state is on disk. Just re-export nothing manually — the scripts set every
`NXF_*` variable themselves. Re-run whichever step you need from the run order above.
