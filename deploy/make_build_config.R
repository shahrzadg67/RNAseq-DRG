#!/usr/bin/env Rscript
# Turn config/stage_status.yaml into app/data/build_config.json for a given target.
#   Rscript deploy/make_build_config.R <staging|prod>
suppressPackageStartupMessages({ library(yaml); library(jsonlite) })
args   <- commandArgs(trailingOnly = TRUE)
target <- if (length(args)) args[1] else "staging"
stopifnot(target %in% c("staging", "prod"))

PROJ <- "/hpf/projects/msalter/sghazis/rnaseq_TUY35595"
y <- yaml::read_yaml(file.path(PROJ, "config", "stage_status.yaml"))

stages_named   <- list()
stages_ordered <- list()
for (s in y$stages) {
  stages_named[[s$key]] <- list(id = s$id, title = s$title,
                                status = s$status, released = isTRUE(s$released))
  stages_ordered[[length(stages_ordered) + 1]] <- list(
    key = s$key, title = s$title, status = s$status, released = isTRUE(s$released))
}

cfg <- list(target = target, project = y$project,
            stages = stages_named, stages_ordered = stages_ordered)

out <- file.path(PROJ, "app", "data", "build_config.json")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write_json(cfg, out, auto_unbox = TRUE, pretty = TRUE)
cat("Wrote", out, "(target:", target, ")\n")
vis <- vapply(y$stages, function(s)
  target == "staging" || isTRUE(s$released), logical(1))
cat("Visible tabs:", paste(c("methods", vapply(y$stages, `[[`, "", "key")[vis & vapply(y$stages,function(s)s$key,"")!="methods"]), collapse=", "), "\n")
