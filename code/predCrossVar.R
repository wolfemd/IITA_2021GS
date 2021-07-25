# functions to ultimately upgrade and supercede those currently
# implemented in "predCrossVar"

# EXPORT THIS TO NAMESPACE
predCrossVars<-function(CrossesToPredict,modelType,
                       AddEffectList,DomEffectList=NULL,
                       predType="VPM",
                       haploMat,recombFreqMat,
                       ncores=1,nBLASthreads=NULL,...){
  starttime<-proc.time()[3]
  # Center posterior distribution of effects
  ## on posterior mean across MCMC samples
  AddEffectList<-purrr::map(AddEffectList,~scale(.,center = T, scale = F));
  if(modelType=="AD"){
    DomEffectList<-purrr::map(DomEffectList,~scale(.,center = T, scale = F)) }

  ## Extract the posterior mean effects vectors
  ## If predType="VPM" this just recovers the original effects
  postMeanAddEffects<-purrr::map(AddEffectList,~attr(.,which = "scaled:center"))
  if(modelType=="AD"){
    postMeanDomEffects<-purrr::map(DomEffectList,~attr(.,which = "scaled:center"))
  } else {
    postMeanDomEffects<-NULL
    }
  parents<-union(CrossesToPredict$sireID,
                 CrossesToPredict$damID)
  haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),]

  if(predType=="VPM"){
    AddEffectList<-NULL;
    if(predType=="VPM" & modelType=="AD"){ DomEffectList<-NULL; } }

  # Set-up a loop over the crosses
  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  crossespredicted<-CrossesToPredict %>%
    mutate(predVars=future_pmap(.,
                                predOneCross,
                                modelType=modelType,
                                haploMat=haploMat,
                                recombFreqMat=recombFreqMat,
                                predType=predType,
                                postMeanAddEffects=postMeanAddEffects,
                                postMeanDomEffects=postMeanDomEffects,
                                AddEffectList=AddEffectList,
                                DomEffectList=DomEffectList,
                                nBLASthreads=nBLASthreads))
  plan(sequential)

  crossespredicted %<>%
    unnest(predVars)

  totcomputetime<-proc.time()[3]-starttime
  print(paste0("Done predicting fam vars. ",
               "Took ",round((totcomputetime)/60,2),
               " mins for ",nrow(crossespredicted)," crosses"))
  return(crossespredicted)
}
# INTERNAL FUNCTION - predict all variance-covariance components for one cross
predOneCross<-function(sireID,damID,modelType,
                       haploMat,recombFreqMat,
                       predType,
                       postMeanAddEffects,
                       postMeanDomEffects=NULL,
                       AddEffectList=NULL,DomEffectList=NULL,
                       nBLASthreads,...){
  starttime<-proc.time()[3]

  if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }

  # Before predicting variances
  # check for and remove SNPs that
  # won't segregate, i.e. are fixed in parents
  ### hopes to save time / mem
  x<-colSums(rbind(haploMat[grep(paste0("^",sireID,"_"),rownames(haploMat)),],
                   haploMat[grep(paste0("^",damID,"_"),rownames(haploMat)),]))
  segsnps2keep<-names(x[x>0 & x<4])

  if(length(segsnps2keep)>0){
    recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep,drop=F]
    haploMat<-haploMat[,segsnps2keep,drop=F]
    postMeanAddEffects<-purrr::map(postMeanAddEffects,~.[segsnps2keep])
    if(modelType=="AD"){
      postMeanDomEffects<-purrr::map(postMeanDomEffects,~.[segsnps2keep]) }

    # calc cross LD matrix
    progenyLD<-calcCrossLD(sireID,damID,recombFreqMat,haploMat)
    rm(recombFreqMat,haploMat); gc()

    # Set-up loop over variance and covarance parameters
    ## Trait variances to-be-predicted
    traits<-names(postMeanAddEffects)
    varcovars<-tibble::tibble(Trait1=traits,
                              Trait2=traits)
    ## If multiple traits
    if(length(traits)>1){
      ## Add covariances to-be-predicted
      varcovars<-dplyr::bind_rows(varcovars, # trait variances
                                  combn(traits,2,simplify = T) %>% # covariances
                                    t(.) %>% #
                                    `colnames<-`(.,c("Trait1","Trait2")) %>%
                                    tibble::as_tibble(.)) }
    varcovars<-varcovars %>%
      dplyr::mutate(predVars=purrr::pmap(.,
                                         predOneCrossVar,
                                         modelType=modelType,
                                         progenyLD=progenyLD,
                                         # haploMat=haploMat,
                                         # recombFreqMat=recombFreqMat,
                                         predType=predType,
                                         postMeanAddEffects=postMeanAddEffects,
                                         postMeanDomEffects=postMeanDomEffects,
                                         AddEffectList=AddEffectList,
                                         DomEffectList=DomEffectList)) %>%
      unnest(predVars)
    computetime<-proc.time()[3]-starttime
    out_thiscross<-tibble(Nsegsnps=length(segsnps2keep),
                          ComputeTime=computetime,
                          predVars=list(varcovars))
  } else {
    computetime<-proc.time()[3]-starttime
    out_thiscross<-tibble(Nsegsnps=length(segsnps2keep),
                          ComputeTime=computetime,
                          predVars=list()) }
  return(out_thiscross)
}
# INTERNAL FUNCTION - predict one variance-covariance component for one cross
predOneCrossVar<-function(Trait1,Trait2,progenyLD,modelType,
                          #haploMat,recombFreqMat,
                          predType,
                          postMeanAddEffects,
                          postMeanDomEffects=NULL,
                          AddEffectList=NULL,
                          DomEffectList=NULL,
                          segsnps2keep=NULL,...){

  if(predType=="PMV"){
    # Posterior Sample Variance-Covariance Matrix of Marker Effects
    postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])
    postVarCovarOfAddEffects<-postVarCovarOfAddEffects[segsnps2keep,
                                                       segsnps2keep,drop=F]
    if(modelType=="AD"){
      postVarCovarOfDomEffects<-(1/(nrow(DomEffectList[[Trait1]])-1))*crossprod(DomEffectList[[Trait1]],DomEffectList[[Trait2]])
      postVarCovarOfDomEffects<-postVarCovarOfDomEffects[segsnps2keep,
                                                         segsnps2keep,drop=F]
      }
  } else {
    postVarCovarOfAddEffects<-NULL;
    postVarCovarOfDomEffects<-NULL;
  }

  ## Predicts variance for one cross and
  ### one variance or covariance paramater (Trait1-Trait2)

  ## Predict cross additive variance
  #### posterior mean (VPM)
  predVarA_vpm<-quadform(D=progenyLD,
                         x=postMeanAddEffects[[Trait1]],
                         y=postMeanAddEffects[[Trait2]])

  #### posterior mean (co)variance (PMV)
  if(predType=="PMV"){ predVarA_pmv<-predVarA_vpm+sum(diag(progenyLD%*%postVarCovarOfAddEffects)) }

  if(modelType=="AD"){
    ## Predict cross dominance variance
    #### VPM
    progenyLDsq<-progenyLD*progenyLD
    predVarD_vpm<-quadform(D=progenyLDsq,
                           x=postMeanDomEffects[[Trait1]],
                           y=postMeanDomEffects[[Trait2]])
    #### PMV
    if(predType=="PMV"){ predVarD_pmv<-predVarD_vpm+sum(diag(progenyLDsq%*%postVarCovarOfDomEffects)) }
  }

  if(modelType=="A"){ rm(progenyLD); gc() }
  if(modelType=="AD"){ rm(progenyLD,progenyLDsq); gc() }

  # Tidy the results
  out<-tibble(predOf="VarA",predVar=predVarA_vpm)
  ### VarA
  if(predType=="PMV"){
    out<-out %>%
      rename(VPM=predVar) %>%
      mutate(PMV=predVarA_pmv) }
  ### If modelType=="AD" - VarD
  if(modelType=="AD"){
    if(predType=="VPM"){
      out %<>% bind_rows(tibble(predOf="VarD",predVar=predVarD_vpm)) }
    if(predType=="PMV"){
      out %<>% bind_rows(tibble(predOf="VarD",
                                 VPM=predVarD_vpm,
                                 PMV=predVarD_pmv)) }
  }
  return(out)
}


calcGameticLD<-function(parentGID,recombFreqMat,haploMat){
  X<-haploMat[paste0(parentGID,c("_HapA","_HapB")),,drop=F]
  p<-colMeans(X);
  D<-recombFreqMat*((0.5*t(X)%*%X)-p%*%t(p));
  return(D) }

calcCrossLD<-function(sireID,damID,recombFreqMat,haploMat){
  return(calcGameticLD(sireID,recombFreqMat,haploMat)+
           calcGameticLD(damID,recombFreqMat,haploMat)) }

# does x %*% D %*% t(x) faster
quadform<-function(D,x,y){ return(as.numeric(colSums(x*(D%*%y)))) }

# Warning: prediction of meanTGV with F-M Eqn. 14.6 appropriate
## only using a+d partition not allele sub. + dom. dev.;
## genotypic NOT classical in terms used by Vitezica et al. 2013.
## For that reason, predCrossMeans has a "predType" not a "modelType" argument
## predType="TGV" uses F-M Eqn. 14.6 and takes add and dom effect
## predType="BV" input should be allele subst. effs, computes mid-parent GEBV
## there is no equivalent to predicting the dominance variance for the mean
## thus the difference from the predCrossVars() function.
predCrossMeans<-function(CrossesToPredict,predType,
                         AddEffectList,DomEffectList=NULL,
                         doseMat,
                         ncores=1,
                         ...){
  #nBLASthreads=NULL,
  means<-tibble(Trait=names(AddEffectList))
  parents<-CrossesToPredict %$% union(sireID,damID)
  doseMat<-doseMat[parents,colnames(AddEffectList[[1]])]

  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

  if(predType=="BV"){
    means<-means %>%
      mutate(predictedMeans=future_map(Trait,function(Trait,
                                                      #nBLASthreads=nBLASthreads,
                                                      ...){

        #if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }

        parentGEBVs<-tcrossprod(doseMat,AddEffectList[[Trait]])

        predmeans<-CrossesToPredict %>%
          left_join(tibble(sireID=rownames(parentGEBVs),
                           sireGEBV=as.numeric(parentGEBVs))) %>%
          left_join(tibble(damID=rownames(parentGEBVs),
                           damGEBV=as.numeric(parentGEBVs))) %>%
          mutate(predOf="MeanBV",
                 predMean=(sireGEBV+damGEBV)/2)
        return(predmeans) }))
    plan(sequential)
  }

  if(predType=="TGV"){
    plan(multisession, workers = ncores)
    means<-means %>%
      mutate(predictedMeans=future_map(Trait,function(Trait,
                                                      #BLASthreads=nBLASthreads,
                                                      ...){

        #if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }

        predmeans<-CrossesToPredict %>%
          mutate(predOf="MeanTGV",
                 predMean=map2_dbl(sireID,damID,function(sireID,damID,...){
                   # Eqn 14.6 from Falconer+MacKay
                   p1<-doseMat[sireID,]/2
                   p2<-doseMat[damID,]/2
                   q<-1-p1
                   y<-p1-p2
                   g<-AddEffectList[[Trait]]*(p1-q-y) + DomEffectList[[Trait]]*((2*p1*q)+y*(p1-q))
                   meanG<-sum(g)
                   return(meanG)}))
        return(predmeans) }))
    plan(sequential)
  }
  means<-means %>%
    unnest(predictedMeans)
  return(means)
}

#' getPropHom
#'
#' Compute the per-individual proportion homozygous.
#' For example, to use as a predictor for a directional dominance model.
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of per-individual proportion homozygous
#' @export
#'
#' @examples
getPropHom<-function(M){
  W<-M; W[which(W==2)]<-0;
  # f = 1 − h/N,
  # where N is the number of SNPs
  f<-1-(rowSums(W)/ncol(W))
  return(f)
}

#' centerDosage
#'
#' Centers dosage matrix, e.g. for use in whole-genome regressions like rrBLUP
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#'
#' @return
#' @export
#'
#' @examples
#' centeredM<-centerDosage(M)
centerDosage<-function(M){
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  Z <- M-2*P
  return(Z)
}

#' dose2domDev
#'
#' Converts a dosage matrix into a matrix of centered dominance deviations.
#' This sets up the "classical" (aka "Statistical") partition of additive dominance in terms of breeding values and dom. deviations.
#' See function dose2domDevGenotypic() to set-up the "genotypic" (aka "biological") partition in terms of genotypic effects.
#' See Vitezica et al. 2013.
#' Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return
#' @export
#'
#' @examples
#' domDev<-dose2domDev(M)
dose2domDev<-function(M){
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  W<-M;
  W[which(W==1)]<-2*P[which(W==1)];
  W[which(W==2)]<-(4*P[which(W==2)]-2);
  W <- W-2*(P^2)
  return(W)
}


#' dose2domDevGenotypic
#'
#' Converts a dosage matrix into a matrix of centered dominance deviations.
#' This sets up the "genotypic" (aka "biological") partition of additive dominance in terms of their genotypic effects instead of in terms of breeding values or dominance deviations.
#' See function dose2domDev() to set-up the "statistical" (aka "classical") partition in terms of genotypic effects.
#' See Vitezica et al. 2013.
#' Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return
#' @export
#'
#' @examples
#' domDev<-dose2domDevGenotypic(M)
dose2domDevGenotypic<-function(M){
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  W<-M; W[which(W==2)]<-0;
  W <- W-(2*P*(1-P))
  return(W)
}

#' getAF
#'
#' get a vector of allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of allele frequencies, names = SNP IDs if in cols of M
#' @export
#'
#' @examples
getAF<-function(M){ colMeans(M,na.rm=T)/2 }

#' getMAF
#'
#' get a vector of _minor_ allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of _minor_ allele frequencies, names = SNP IDs if in cols of M
#' @export
#'
#' @examples
getMAF<-function(M){
  freq<-colMeans(M, na.rm=T)/2; maf<-freq;
  maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
  return(maf) }

#' maf_filter
#'
#' filter a dosage matrix by minor allele frequence.
#' get a vector of allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#' @param thresh threshold value. Columns of M with maf<thresh will be removed
#' @return dosage matrix potentially with columns removed
#' @export
#'
#' @examples
maf_filter<-function(M,thresh){
  freq<-colMeans(M, na.rm=T)/2; maf<-freq;
  maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
  snps1<-M[,which(maf>thresh)];
  return(snps1) }

#' remove_invariant
#'
#' filter a dosage matrix, removing invariant markers. Removes e.g. cases where MAF=0.5 but all dosages == 1 (het).
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#' @param thresh threshold value. Columns of M with maf<thresh will be removed
#' @return dosage matrix potentially with columns removed
#' @export
#'
#' @examples
remove_invariant<-function(M){
  snps1<-M[ ,apply(M, 2, var) != 0]
  return(snps1) }

#' backsolveSNPeff
#'
#' From the GBLUP solutions and a centered SNP matrix backsolve SNP effects
#'
#' @param Z Centered marker matrix (dominance deviations must also be centered)
#' @param g The solutions (blups, i.e. GEBVs) from the GBLUP model
#'
#' @return matrix of SNP effects matching RR-BLUP / SNP-BLUP
#' @export
#'
#' @examples
#' A<-kinship(M,type="add")
#' trainingDF %<>% dplyr::mutate(ga=factor(as.character(id),
#'                                         levels=rownames(A)),
#'                               gd=ga)
#' gblup<-mmer(pheno~1,
#'             random=~vs(ga,Gu = A),
#'             weights=WT,
#'             data=trainingDF,verbose = T)
#' ga<-as.matrix(gblup$U$`u:ga`$pheno,ncol=1)
#' Za<-centerDosage(M)
#' snpeff<-backsolveSNPeff(Za,ga)
backsolveSNPeff<-function(Z,g){
  ZZt<-tcrossprod(Z)
  diag(ZZt)<-diag(ZZt)+1e-8
  bslvEffs<-crossprod(Z,solve(ZZt))%*%g
  return(bslvEffs)
}

#' genoVarCovarMatFunc
#'
#' Compute the *p SNPs*-by-*p SNPs* variance-covariance matrix of SNP dosages.
#' This is an estimator of the LD between loci within a given population.
#'
#' @param Z column-centered matrix of SNP dosages. Assumes SNPs in Z were originally coded 0, 1, 2 were column centered.
#'
#' @return
#'  NOTE: this matrix is going to be big in practice.
#'  The *p SNPs*-by-*p SNPs* variance-covariance matrix of SNP dosages.
#'  may be worth computing in an R session using multi-threaded BLAS
#' @export
#'
#' @examples
#' Z<-centerDosage(M)
#' genoVarCovarMat<-genoVarCovarMatFunc(Z)
genoVarCovarMatFunc<-function(Z){
  SigmaM<-1/nrow(Z)*t(Z)%*%Z
  return(SigmaM)
}


#' genmap2recombfreq
#'
#' Compute the pairwise recombination frequencies between all loci from genetic map positions.
#'
#' @param m vector of centiMorgan-scale genetic positions. names(m) correspond to a SNP_ID. Since m potentially contains all chromosomes, sets recomb. freq. b/t chrom. to 0.5
#' @param nChr number of chromosomes
#'
#' @details names(m) must be formatted as "chr"_"id" with "chr" being integer. For example: 2_QTL1 for a locus on chr. 2.
#' May be worth computing in an R session using multi-threaded BLAS.
#' @return potentially really large matrix of pairwise recombination frequencies between loci
#' @export
#'
#' @examples
genmap2recombfreq<-function(m,nChr){
  d<-as.matrix(dist(m,upper=T,diag = T,method='manhattan'))
  c1<-0.5*(1-exp(-2*d))
  # Since m contains all chromosomes, set recomb. freq. b/t chrom. to 0.5
  for(i in 1:nChr){
    c1[grepl(paste0("^",i,"_"),rownames(c1)),!grepl(paste0("^",i,"_"),colnames(c1))]<-0.5
    c1[!grepl(paste0("^",i,"_"),rownames(c1)),grepl(paste0("^",i,"_"),colnames(c1))]<-0.5
  }
  return(c1)
}

#' crosses2predict
#'
#' Make a data.frame of all pairwise matings given a vector of parent IDs.
#' Include selfs. No reciprocal crosses, i.e. use as male == use as female.
#' Diagonal and upper-triangle of mating matrix.
#'
#' @param parents
#'
#' @return tibble, two columns, sireID and damID, all pairwise crosses (see details).
#' @export
#'
#' @examples
crosses2predict<-function(parents){
  CrossesToPredict<-matrix(NA,nrow=length(parents),ncol=length(parents))
  CrossesToPredict[upper.tri(CrossesToPredict,diag = T)]<-1
  rownames(CrossesToPredict)<-colnames(CrossesToPredict)<-parents
  CrossesToPredict<-CrossesToPredict %>%
    as.data.frame %>%
    tibble::rownames_to_column(var = "sireID") %>%
    tidyr::pivot_longer(cols = (-sireID), names_to = "damID", values_to = "keep") %>%
    dplyr::filter(keep==1) %>%
    dplyr::select(-keep)
  return(CrossesToPredict)
}

#' intensity
#'
#' Compute the standardized selection intensity from the proportion selection.
#'
#' @param propSel proportion selection
#'
#' @return
#' @export
#'
#' @examples
intensity<-function(propSel){ dnorm(qnorm(1-propSel))/propSel } # standardized selection intensity
