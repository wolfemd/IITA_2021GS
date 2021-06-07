# functions to ultimately upgrade and supercede those currently
# implemented in "predCrossVar"

# EXPORT THIS TO NAMESPACE
predCrossVars<-function(CrossesToPredict,modelType,
                       AddEffectList,DomEffectList=NULL,
                       predType="VPM",
                       haploMat,recombFreqMat,
                       ncores=1,...){
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
  require(furrr); options(mc.cores=ncores); plan(multicore)
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
                                DomEffectList=DomEffectList)) %>%
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
                       AddEffectList=NULL,DomEffectList=NULL,...){
  starttime<-proc.time()[3]
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
                                         haploMat=haploMat,
                                         recombFreqMat=recombFreqMat,
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
                          haploMat,recombFreqMat,
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


predCrossMeans<-function(CrossesToPredict,modelType,
                         AddEffectList,DomEffectList=NULL,
                         doseMat,
                         ncores=1,...){

  means<-tibble(Trait=names(AddEffectList))
  parents<-CrossesToPredict %$% union(sireID,damID)
  doseMat<-doseMat[parents,colnames(AddEffectList[[1]])]

  require(furrr); options(mc.cores=ncores); plan(multicore)
  means<-means %>%
    mutate(predictedMeans=future_map(Trait,function(Trait,...){
      parentGEBVs<-tcrossprod(doseMat,AddEffectList[[Trait]])

      predmeans<-CrossesToPredict %>%
        left_join(tibble(sireID=rownames(parentGEBVs),
                         sireGEBV=as.numeric(parentGEBVs))) %>%
        left_join(tibble(damID=rownames(parentGEBVs),
                         damGEBV=as.numeric(parentGEBVs))) %>%
        mutate(predOf="MeanBV",predMean=(sireGEBV+damGEBV)/2)
      if(modelType=="AD"){
        predmeans %<>%
          bind_rows(predmeans %>%
                      mutate(predOf="MeanTGV",
                             predMean=map2_dbl(sireID,damID,function(sireID,damID,...){
                               # Eqn 14.6 from Falconer+MacKay
                               p1<-doseMat[sireID,]/2
                               p2<-doseMat[damID,]/2
                               q<-1-p1
                               y<-p1-p2
                               g<-AddEffectList[[Trait]]*(p1-q-y) + DomEffectList[[Trait]]*((2*p1*q)+y*(p1-q))
                               meanG<-sum(g)
                               return(meanG)})))
      }
      return(predmeans)})) %>%
    unnest(predictedMeans)
  return(means) }

