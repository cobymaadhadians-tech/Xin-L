suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("ANALYSIS_ROOT")
if (!nzchar(root)) {
  root <- normalizePath(getwd(), mustWork = FALSE)
}
out_dir <- file.path(root,"results/d6b_direct_contrast_20260905")
capture_path <- Sys.getenv("D6B_CAPTURE_RDS", file.path(out_dir,"D6B_direct_textinject_capture.rds"))
sem_path <- file.path(root,"results/d6b_published_pi_braincell_stratified_gsem_20260905/D6B_S_LDSCOutput.rds")
enrichment_path <- file.path(root,"results/d6b_published_pi_braincell_stratified_gsem_20260905/D6B_enrichment_long.tsv")
stopifnot(file.exists(capture_path), file.exists(sem_path), file.exists(enrichment_path))
cap <- readRDS(capture_path)
s <- readRDS(sem_path)
en <- fread(enrichment_path)
p <- c("F_COMP~~SCZ_PGC","F_INT~~SCZ_PGC","F_SUD~~SCZ_PGC","F_COMP~~SCZ_FG","F_INT~~SCZ_FG","F_SUD~~SCZ_FG")
en <- en[target_parameter %in% p & target_parameter != "BASELINE"]
ann <- unique(en[,.(Annotation,annotation_key,annotation_basename,annotation_class,published_label)])
setorder(ann, Annotation)
alpha <- 0.05/(nrow(ann)*3)
rows <- vector("list", nrow(ann)*3); k <- 0L
for (i in seq_len(nrow(ann))) {
  a <- ann$Annotation[i]
  ii <- which(vapply(cap,function(x) identical(as.character(names(s$S)[x$n]),a),logical(1))); cc <- if(length(ii)) cap[[ii[1]]] else list()
  er <- en[Annotation==a]
  warn <- any(er$Annotation_model_warning %in% TRUE | er$Warning %in% c(TRUE,"TRUE","1","1.0"))
  err <- any(er$Annotation_model_error %in% TRUE | er$Error %in% c(TRUE,"TRUE","1","1.0"))
  status <- if (length(cc)==0 || is.null(cc$covariance)) "no sandwich capture" else "sandwich captured"
  if (length(cc)>0 && !is.null(cc$covariance) && length(cc$estimate)==6 && length(cc$base_est)==6) {
    idx <- match(p, cc$map)
    if (all(!is.na(idx))) {
      est <- as.numeric(cc$estimate)
      base <- as.numeric(cc$base_est)
      prop <- as.numeric(cc$prop)
      sc <- 1/(base*prop)
      cov0 <- as.matrix(cc$covariance)[idx,idx,drop=FALSE]
      enr <- est*sc
      ce <- diag(sc) %*% cov0 %*% diag(sc)
    } else { enr <- rep(NA_real_,6); ce <- matrix(NA_real_,6,6) }
  } else { enr <- rep(NA_real_,6); ce <- matrix(NA_real_,6,6) }
  for (d in 1:3) {
    k <- k+1L
    ip <- d; ig <- d+3
    delta <- enr[ig]-enr[ip]
    vv <- ce[ig,ig]+ce[ip,ip]-2*ce[ig,ip]
    se <- if (is.finite(vv) && vv>=0) sqrt(vv) else NA_real_
    zz <- if (is.finite(delta) && is.finite(se) && se>0) delta/se else NA_real_
    pp <- if (is.finite(zz)) 2*pnorm(-abs(zz)) else NA_real_
    lr <- er[target_parameter==p[ip]]
    fr <- er[target_parameter==p[ig]]
    rows[[k]] <- data.table(
      Annotation=a, annotation_key=ann$annotation_key[i], annotation_basename=ann$annotation_basename[i],
      annotation_class=ann$annotation_class[i], published_label=ann$published_label[i],
      domain=c("COMP","INT","SUD")[d],
      PGC_enrichment=enr[ip], FinnGen_enrichment=enr[ig],
      delta_FinnGen_minus_PGC=delta, SE=se, z=zz, P_two_sided=pp,
      alpha_492=alpha, Bonferroni_significant_492=isTRUE(is.finite(pp) && pp<=alpha),
      Annotation_model_warning=warn, Annotation_model_error=err,
      QC_warning_free=isTRUE(!warn && !err && is.finite(pp) && is.finite(se)),
      contrast_status=status,
      PGC_enrichment_from_enrich=if(nrow(lr)) lr$Enrichment[1] else NA_real_,
      FinnGen_enrichment_from_enrich=if(nrow(fr)) fr$Enrichment[1] else NA_real_
    )
  }
}
out <- rbindlist(rows,fill=TRUE)
out[,Bonferroni_significant_492_warning_free := Bonferroni_significant_492 & QC_warning_free]
fwrite(out,file.path(out_dir,"D6B_direct_annotation_contrast.tsv"),sep="\t",na="NA",quote=FALSE)
fwrite(data.table(metric=c("n_annotations","n_rows","n_capture_annotations","n_finite_contrasts","n_significant_492","n_significant_492_warning_free"),
                  value=c(nrow(ann),nrow(out),length(cap),sum(is.finite(out$P_two_sided)),sum(out$Bonferroni_significant_492),sum(out$Bonferroni_significant_492_warning_free))),
       file.path(out_dir,"D6B_direct_annotation_contrast_summary.tsv"),sep="\t")
cat("DIRECT_TABLE_DONE annotations=",nrow(ann)," rows=",nrow(out)," captures=",length(cap),
    " finite=",sum(is.finite(out$P_two_sided))," sig=",sum(out$Bonferroni_significant_492),
    " wf_sig=",sum(out$Bonferroni_significant_492_warning_free),"\n",sep="")
writeLines(capture.output(sessionInfo()), file.path(out_dir,"D6B_direct_reproduction_sessionInfo.txt"))
