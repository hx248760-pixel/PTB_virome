# Supplementary fully nested random-forest analysis for the PTB virome study.
# Run from the repository root:
#   Rscript script/06_random_forest_nested.R
# Or specify project root, output directory, and backend:
#   Rscript script/06_random_forest_nested.R . results/nested_rf_v2 randomForest
# Interactive R (source defines the function without starting computation):
#   source("script/06_random_forest_nested.R")
#   res <- run_nested_rf(project_dir=".")
#
# Required packages: Maaslin2 (Bioconductor), pROC, randomForest (CRAN).
# Install once if needed, not automatically during analysis:
#   install.packages(c("BiocManager", "pROC", "randomForest"))
#   BiocManager::install("Maaslin2")
# ranger is optional and NOT needed for the reported randomForest results.
# For exact reproduction, use the package versions in the original sessionInfo.txt.
# Changing backend or software versions may change numerical results.
#
# Required discovery-cohort inputs:
#   data/metadata.tsv: Sample, Group (Health/Patient), Sex, Age; no missing values.
#   data/votu_abundance.tsv: features in rows, samples in columns; first column ID.
#   data/bacterial_species_abundance.tsv: same layout as the viral table.
# Supply full abundance tables, not features preselected using all group labels.
# Threshold 0.01 is in the ORIGINAL abundance-table units (not necessarily 1%).
# Age and sex adjust MaAsLin2 selection; they are not RF predictors.
#
# Reported settings: 5 outer folds x 10 repetitions; 3 inner folds.
# Grid: top 10 / top 30 / all selected features x nodesize 1 / 5.
# Ranking: training-only Gini importance with 300 trees; classification: 1000 trees.
# No additional class weights or resampling; stratified partitions.
# Empty feature sets predict the training-set PTB proportion.
# Explicit seeds are unchanged from the server script used for the manuscript.
# External-cohort validation is performed by a separate script.
#
# Default output: results/nested_rf_v2/ (created automatically).
#   oof_predictions.tsv: one row per sample, repetition, and model; Actual=1 PTB/0 HC.
#   summary.tsv, metrics_per_repeat.tsv: AUC/Brier metrics (AUC in 0-1 units).
#   outer_model_choices.tsv, inner_tuning.tsv, selected_features.tsv, split_audit.tsv.
#   feature_name_mapping.tsv: encoded IDs to original feature names.
#   nested_ROC.pdf/png, calibration.pdf/png, AUC_by_repeat.pdf.
#   sessionInfo.txt, analysis_config.rds, engine.txt; fold checkpoints and Maaslin2 logs.
# SD_AUC describes split-to-split variability, not a population confidence interval.
# Calibration averages participant OOF probabilities across repetitions into 5 bins.
# ROC curves average interpolated per-repetition OOF curves.
#
# Resume: rerun the SAME command with unchanged inputs, settings, code, and packages.
# Completed outer folds are read from checkpoints; an interrupted fold is recomputed.
# A run signature prevents mixing incompatible runs. Do not delete it to bypass checks.
# Use a new output directory for changed settings (including the number of repetitions).
# The serial statistical implementation below matches the original server script.

# Fully nested RF. Source this file to define run_nested_rf; no analysis starts on source.
# Dependencies: Maaslin2, pROC, and ONE RF backend: randomForest OR ranger.
# Ranking uses training-only Gini/impurity importance, not the old permutation ranking.
run_nested_rf <- function(project_dir='.', out_dir=NULL,
                          backend='randomForest', repeats=10L, trees=1000L, threads=4L,
                          min_abundance=0.01,min_prevalence=0.01,q_threshold=0.25) {
  root <- normalizePath(project_dir,mustWork=TRUE)
  out <- if(is.null(out_dir)) file.path(root,'results/nested_rf_v2') else out_dir
  backend <- match.arg(backend,c('randomForest','ranger'))
  for(p in c('Maaslin2','pROC',backend)) {
    if(!requireNamespace(p,quietly=TRUE)) stop('Missing package: ',p,
                                               '. ranger is not required when backend="randomForest".')
  }
  stopifnot(repeats>=1,trees>=10,threads>=1,min_abundance>=0,
            min_prevalence>=0,min_prevalence<1,q_threshold>0,q_threshold<1)
  metadata_file <- file.path(root,'data/metadata.tsv')
  viral_file <- file.path(root,'data','votu_abundance.tsv')
  bacterial_file <- file.path(root,'data','bacterial_species_abundance.tsv')
  for(f in c(metadata_file,viral_file,bacterial_file)) if(!file.exists(f)) stop('Missing input: ',f)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  REPEATS <- as.integer(repeats)
  OUTER <- 5L; INNER <- 3L; TREES <- as.integer(trees)
  THREADS <- as.integer(threads)
  # Thresholds are in ORIGINAL TABLE UNITS; confirm the input scale.
  MIN_ABUND <- min_abundance; MIN_PREV <- min_prevalence; Q <- q_threshold
  meta <- read.delim(metadata_file,check.names=FALSE)
  stopifnot(!anyDuplicated(meta$Sample),all(meta$Group %in% c('Health','Patient')),
            !anyNA(meta[,c('Sample','Group','Sex','Age')]))
  meta$Age <- suppressWarnings(as.numeric(meta$Age))
  stopifnot(all(is.finite(meta$Age)),nrow(meta)>=20)
  rownames(meta)<-meta$Sample
  meta$Group <- factor(meta$Group,levels=c('Patient','Health'))
  meta$Sex <- factor(meta$Sex)
  y <- as.integer(meta$Group=='Patient')
  stopifnot(all(table(factor(y,levels=0:1))>=10))
  load_table <- function(fn,prefix) {
    d<-read.delim(fn,row.names=1,check.names=FALSE)
    stopifnot(setequal(colnames(d),meta$Sample),!anyDuplicated(rownames(d)))
    x<-t(as.matrix(d[,meta$Sample,drop=FALSE])); storage.mode(x)<-'double'
    stopifnot(all(is.finite(x)),all(x>=0))
    colnames(x)<-paste0(prefix,colnames(x)); x
  }
  V<-load_table(viral_file,'V__')
  B<-load_table(bacterial_file,'B__')
  message('Viral column sums: ',paste(range(rowSums(V)),collapse=' to '))
  message('Bacterial column sums: ',paste(range(rowSums(B)),collapse=' to '))
  feature_map <- data.frame(encoded=c(colnames(V),colnames(B)), original=sub('^[VB]__','',c(colnames(V),colnames(B))), type=c(rep('Virus',ncol(V)),rep('Bacteria',ncol(B))))
  # Safe feature IDs prevent MaAsLin2 name sanitization from breaking matching.
  colnames(V) <- sprintf('V%06d',seq_len(ncol(V)))
  colnames(B) <- sprintf('B%06d',seq_len(ncol(B)))
  feature_map$encoded <- c(colnames(V),colnames(B))
  write.table(feature_map,file.path(out,'feature_name_mapping.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  folds <- function(idx,k,seed) {
    set.seed(seed); z<-integer(length(idx))
    for(cl in 0:1) { at<-which(y[idx]==cl); z[at]<-sample(rep(seq_len(k),length.out=length(at))) }
    lapply(seq_len(k),function(j) idx[z==j])
  }
  # Every fitted feature transform stores its TRAINING pseudocount.
  # MaAsLin2 performs its own LOG transform only on this training partition.
  transform_fit <- function(x) {
    eps<-min(x[x>0])/2
    if(!is.finite(eps)) stop('No positive abundances')
    list(x=log2(x+eps),eps=eps)
  }
  rf <- function(x,yy,node=1L,seed=1L,ranking=FALSE) {
    ntree <- if(ranking) min(300L,TREES) else TREES
    mtry <- max(1L,floor(sqrt(ncol(x))))
    set.seed(seed)
    if(backend=='randomForest') {
      randomForest::randomForest(x=as.data.frame(x),y=factor(yy,levels=0:1),
                                 ntree=ntree,mtry=mtry,nodesize=node,importance=FALSE)
    } else {
      ranger::ranger(x=as.data.frame(x),y=factor(yy,levels=0:1),probability=TRUE,
                     num.trees=ntree,mtry=mtry,min.node.size=node,
                     importance=if(ranking) 'impurity' else 'none',seed=seed,num.threads=THREADS)
    }
  }
  rf_importance <- function(fit) {
    if(backend=='randomForest') {
      im <- randomForest::importance(fit,type=2)
      setNames(as.numeric(im[,1]),rownames(im))
    } else fit$variable.importance
  }
  rf_probability <- function(fit,x) {
    if(backend=='randomForest') as.numeric(predict(fit,newdata=as.data.frame(x),type='prob')[,'1'])
    else as.numeric(predict(fit,data=as.data.frame(x))$predictions[,'1'])
  }
  select_one <- function(x,idx,tag) {
    md<-meta[idx,,drop=FALSE]
    # Remove absent and constant features using this training partition only.
    xx<-x[idx,,drop=FALSE]
    keep<-colMeans(xx>MIN_ABUND)>MIN_PREV & apply(xx,2,function(a) length(unique(a))>1)
    xx<-xx[,keep,drop=FALSE]
    if(!ncol(xx)) return(character())
    d<-file.path(out,'maaslin',tag); dir.create(d,recursive=TRUE,showWarnings=FALSE)
    effects<-c('Group',if(length(unique(md$Sex))>1) 'Sex',if(sd(md$Age)>0) 'Age')
    Maaslin2::Maaslin2(input_data=as.data.frame(xx),input_metadata=md,output=d,
                       fixed_effects=effects,reference=c('Group,Patient'),normalization='NONE',
                       transform='LOG',standardize=TRUE,analysis_method='LM',min_abundance=MIN_ABUND,
                       min_prevalence=MIN_PREV,max_significance=Q,correction='BH',
                       plot_heatmap=FALSE,plot_scatter=FALSE,cores=1)
    a<-read.delim(file.path(d,'all_results.tsv'),check.names=FALSE)
    s <- unique(a$feature[a$metadata=='Group' & is.finite(a$qval) & a$qval<Q])
    if(any(!s %in% colnames(xx))) stop('MaAsLin2 feature names changed unexpectedly')
    s
  }
  prepare <- function(idx,tag,seed) {
    sV<-select_one(V,idx,paste0(tag,'_V')); sB<-select_one(B,idx,paste0(tag,'_B'))
    sets<-list(Virus=sV,Bacteria=sB,Combined=c(sV,sB))
    lapply(sets,function(s) {
      if(!length(s)) return(list(features=character(),eps=NA_real_))
      raw<-cbind(V[,intersect(s,colnames(V)),drop=FALSE],B[,intersect(s,colnames(B)),drop=FALSE])
      t<-transform_fit(raw[idx,,drop=FALSE]); fit<-rf(t$x,y[idx],seed=seed,ranking=TRUE)
      list(features=names(sort(rf_importance(fit),decreasing=TRUE)),eps=t$eps)
    })
  }
  X<-cbind(V,B)
  grid<-expand.grid(top_n=c(10L,30L,Inf),node=c(1L,5L))
  model_names<-c('Virus','Bacteria','Combined')
  predict_cfg <- function(pre,train,test,g,seed) {
    s<-head(pre$features,min(length(pre$features),g$top_n))
    if(!length(s)) return(list(prob=rep(mean(y[train]),length(test)),features=s))
    xt<-log2(X[train,s,drop=FALSE]+pre$eps); xv<-log2(X[test,s,drop=FALSE]+pre$eps)
    fit<-rf(xt,y[train],node=g$node,seed=seed)
    list(prob=rf_probability(fit,xv),features=s)
  }
  signature <- list(input_md5=unname(tools::md5sum(c(metadata_file,viral_file,bacterial_file))),
                    backend=backend,repeats=REPEATS,trees=TREES,threads=THREADS,min_abundance=MIN_ABUND,
                    min_prevalence=MIN_PREV,q=Q,outer=OUTER,inner=INNER,grid=grid,
                    versions=sapply(c('Maaslin2',backend,'pROC'),function(p) as.character(utils::packageVersion(p))),
                    code=paste(deparse(body(run_nested_rf)),collapse='\n'))
  guard<-file.path(out,'run_signature.rds')
  if(file.exists(guard)) {
    if(!identical(signature,readRDS(guard))) stop('Output directory contains a different run. Choose a NEW out_dir.')
  } else {
    if(length(list.files(out,pattern='^r[0-9]+_o[0-9]+\\.rds$'))) stop('Unguarded old checkpoints found. Choose a NEW out_dir.')
    saveRDS(signature,guard)
  }
  records<-list(); choices<-list(); selected<-list(); tuning<-list(); audit<-list()
  for(rep in seq_len(REPEATS)) {
    outer<-folds(seq_along(y),OUTER,20261008L+rep)
    for(of in seq_len(OUTER)) {
      test<-outer[[of]]; train<-setdiff(seq_along(y),test)
      stopifnot(!length(intersect(train,test)))
      tag<-sprintf('r%02d_o%d',rep,of)
      # Checkpoint avoids recomputing completed outer folds in this run.
      cp<-file.path(out,paste0(tag,'.rds'))
      if(file.exists(cp)) { z<-readRDS(cp) } else {
        inner<-folds(train,INNER,10000L+rep*100+of)
        inner_pre<-lapply(seq_len(INNER),function(j) prepare(setdiff(train,inner[[j]]),paste0(tag,'_i',j),rep*1000+of*10+j))
        outer_pre<-prepare(train,paste0(tag,'_outer'),rep*1000+of)
        z<-list(records=list(),choices=list(),selected=list(),tuning=list(),audit=list())
        for(m in model_names) {
          scores<-numeric(nrow(grid))
          for(gi in seq_len(nrow(grid))) {
            ip<-rep(NA_real_,length(train)); names(ip)<-train
            for(j in seq_len(INNER)) {
              itr<-setdiff(train,inner[[j]]); iva<-inner[[j]]
              ip[as.character(iva)]<-predict_cfg(inner_pre[[j]][[m]],itr,iva,grid[gi,],rep*1000+of*10+j)$prob
            }
            stopifnot(!anyNA(ip)); scores[gi]<-as.numeric(pROC::auc(pROC::roc(y[train],ip,levels=0:1,direction='<',quiet=TRUE)))
          }
          best<-which.max(scores); p<-predict_cfg(outer_pre[[m]],train,test,grid[best,],rep*1000+of)
          z$records[[m]]<-data.frame(Sample=meta$Sample[test],Actual=y[test],repeat_id=rep,outer_fold=of,model=m,probability=p$prob)
          z$choices[[m]]<-data.frame(repeat_id=rep,outer_fold=of,model=m,top_n=grid$top_n[best],min_node=grid$node[best],inner_auc=scores[best],n_selected=length(p$features),prior_only=!length(p$features))
          z$selected[[m]]<-data.frame(repeat_id=rep,outer_fold=of,model=m,feature=if(length(p$features)) p$features else NA_character_)
          z$tuning[[m]]<-cbind(repeat_id=rep,outer_fold=of,model=m,grid,inner_auc=scores)
        }
        z$audit<-data.frame(repeat_id=rep,outer_fold=of,Sample=meta$Sample,partition=ifelse(seq_along(y)%in%test,'test','train'))
        saveRDS(z,cp)
      }
      records<-c(records,z$records); choices<-c(choices,z$choices); selected<-c(selected,z$selected); tuning<-c(tuning,z$tuning); audit[[tag]]<-z$audit
      message('Completed ',tag)
    }
  }
  oof<-do.call(rbind,records)
  stopifnot(nrow(oof)==nrow(meta)*REPEATS*3,!anyDuplicated(oof[,c('Sample','repeat_id','model')]))
  write_tsv<-function(x,f) write.table(x,file.path(out,f),sep='\t',quote=FALSE,row.names=FALSE)
  write_tsv(oof,'oof_predictions.tsv');write_tsv(do.call(rbind,choices),'outer_model_choices.tsv')
  write_tsv(do.call(rbind,selected),'selected_features.tsv');write_tsv(do.call(rbind,tuning),'inner_tuning.tsv');write_tsv(do.call(rbind,audit),'split_audit.tsv')
  metrics<-do.call(rbind,lapply(split(oof,list(oof$model,oof$repeat_id),drop=TRUE),function(d) data.frame(model=d$model[1],repeat_id=d$repeat_id[1],AUC=as.numeric(pROC::auc(pROC::roc(d$Actual,d$probability,levels=0:1,direction='<',quiet=TRUE))),Brier=mean((d$Actual-d$probability)^2))))
  summary<-do.call(rbind,lapply(split(metrics,metrics$model),function(d) data.frame(model=d$model[1],mean_AUC=mean(d$AUC),SD_AUC=sd(d$AUC),min_AUC=min(d$AUC),max_AUC=max(d$AUC),mean_Brier=mean(d$Brier))))
  write_tsv(metrics,'metrics_per_repeat.tsv');write_tsv(summary,'summary.tsv')
  cols<-c(Virus='#9fcc56',Bacteria='#869bc3',Combined='#bc6363')
  roc_plot<-function() {
    plot(c(0,1),c(0,1),type='n',xlab='False positive rate',ylab='True positive rate',asp=1,main='Fully nested repeated cross-validation');abline(0,1,lty=2,col='grey60')
    fpr<-seq(0,1,length.out=201)
    for(m in model_names) {
      curves<-sapply(seq_len(REPEATS),function(r) {
        d<-oof[oof$model==m & oof$repeat_id==r,];ro<-pROC::roc(d$Actual,d$probability,levels=0:1,direction='<',quiet=TRUE)
        approx(1-ro$specificities,ro$sensitivities,xout=fpr,ties=max,rule=2)$y
      })
      lines(fpr,rowMeans(matrix(curves,nrow=length(fpr))),col=cols[m],lwd=2)
    }
    s<-summary[match(model_names,summary$model),]
    legend('bottomright',legend=sprintf('%s  AUC %.3f (SD %.3f)',s$model,s$mean_AUC,s$SD_AUC),col=cols[model_names],lwd=2,bty='n',cex=.85)
  }
  cal_plot<-function() {
    par(mfrow=c(1,3),mar=c(4,4,3,1))
    for(m in model_names) {
      d<-oof[oof$model==m,];avg<-aggregate(probability~Sample+Actual,d,mean)
      # Display each participant once; probabilities averaged across repetitions.
      bin<-pmin(5,ceiling(rank(avg$probability,ties.method='first')/nrow(avg)*5))
      q<-aggregate(cbind(pred=avg$probability,observed=avg$Actual),list(bin=bin),mean)
      plot(q$pred,q$observed,xlim=0:1,ylim=0:1,xlab='Mean predicted probability',ylab='Observed PTB proportion',main=m,pch=19,col=cols[m]);abline(0,1,lty=2,col='grey60');lines(q$pred,q$observed,col=cols[m])
    }
  }
  pdf(file.path(out,'nested_ROC.pdf'),7,6);roc_plot();dev.off()
  png(file.path(out,'nested_ROC.png'),2100,1800,res=300);roc_plot();dev.off()
  pdf(file.path(out,'calibration.pdf'),12,4);cal_plot();dev.off()
  png(file.path(out,'calibration.png'),3000,1000,res=250);cal_plot();dev.off()
  pdf(file.path(out,'AUC_by_repeat.pdf'),7,5);boxplot(AUC~model,metrics,ylab='OOF AUC per repetition',col=cols[sort(model_names)],ylim=c(0,1));dev.off()
  writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
  saveRDS(list(grid=grid,repeats=REPEATS,outer=OUTER,inner=INNER,trees=TREES,min_abundance=MIN_ABUND,min_prevalence=MIN_PREV,q=Q),file.path(out,'analysis_config.rds'))
  writeLines(c(paste0('backend=',backend)),file.path(out,'engine.txt'))
  print(summary)
  invisible(list(summary=summary,metrics=metrics,oof=oof,output_dir=out))
}
# CLI: Rscript script/06_random_forest_nested.R PROJECT_DIR OUTPUT_DIR BACKEND
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_nested_rf(project_dir=if(length(args)>=1) args[1] else '.',
                out_dir=if(length(args)>=2) args[2] else NULL,
                backend=if(length(args)>=3) args[3] else 'randomForest')
}
