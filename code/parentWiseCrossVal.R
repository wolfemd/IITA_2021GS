# modelType="A", additive-only, predicts mean and variances for breeding values [BVs])
# modelType="AD", the "classic" add-dom model, predicts family mean and variance
#### for breeding value based on allele sub. effects.
#### predicts family variance-covariances for TGVs (BVs+DDs).
#### doesn't predict family mean TGV
# modelType="DirDom", the "genotypic" add-dom model with prop. homozygous
#### fit as a fixed-effect, to estimate a genome-wide inbreeding effect.
#### obtains add-dom effects, computes allele sub effects (alpha = a + d(q-p))
#### predicts family mean and covars for BVs using allele sub. effects
#### predicts covars AND the means for TGVs using directly the add-dom effects.
#### estimated linear effect of overall homozygosity (b), interpreted as
#### inbreeding depression or heterosis depending on
#### its direction relative to each trait (Xiang et al. 2016).
#### The estimated genome-wide fixed-effect of homozygosity ($b$) can be
#### incorporated into the predicted means and variances by first dividing
#### by the number of effects ($p$) and subtracting that value from
#### the vector of dominance effects ($d*$), to get $d=d*-\frac{b}{p}$

runParentWiseCrossVal<-function(nrepeats,nfolds,seed=NULL,modelType,
                                ncores=1,nBLASthreads=NULL,
                                outName=NULL,
                                ped=ped,gid="GID",blups,
                                dosages,grms,haploMat,recombFreqMat,
                                selInd,SIwts = NULL,...){
  initime<-proc.time()[3]

  ## Make parent-wise folds
  parentfolds<-makeParentFolds(ped=ped,gid="GID",
                               nrepeats=nrepeats,
                               nfolds=nfolds,
                               seed=seed)
  print("Set-up parent-wise cross-validation folds")

  ## Get univariate REML marker effects
  #### modelType can be "A","AD" or "DirDom"
  print("Fitting models to get marker effects")
  starttime<-proc.time()[3]
  markEffs<-getMarkEffs(parentfolds,blups=blups,gid=gid,modelType=modelType,
                        grms=grms,dosages=dosages,
                        ncores=ncores,nBLASthreads=nBLASthreads)
  print(paste0("Marker-effects Computed. Took  ",
               round((proc.time()[3] - starttime)/60/60,5)," hrs"))

  ## Predict cross variances
  print("Predicting cross variances and covariances")
  starttime<-proc.time()[3]
  cvPredVars<-predictCrossVars(modelType=modelType,ncores=ncores,
                               snpeffs=markEffs,
                               parentfolds=parentfolds,
                               haploMat=haploMat,
                               recombFreqMat=recombFreqMat)
  print(paste0("Cross variance parameters predicted. Took  ",
               round((proc.time()[3] - starttime)/60/60,5)," hrs"))

  print("Predicting cross means")
  starttime<-proc.time()[3]
  cvPredMeans<-predictCrossMeans(modelType=modelType,ncores=ncores,
                                 snpeffs=markEffs,
                                 parentfolds=parentfolds,
                                 doseMat=dosages)
  print(paste0("Cross means predicted. Took  ",
               round((proc.time()[3] - starttime)/60/60,5)," hrs"))

  print("Compute prediction accuracies and wrap up.")
  ## Variance prediction accuracies
  starttime<-proc.time()[3]
  varPredAcc<-varPredAccuracy(modelType = modelType,
                                   crossValOut = cvPredVars,
                                   snpeffs = markEffs,
                                   ped = ped,selInd = selInd,SIwts = SIwts)

  ## Mean prediction accuracies
  meanPredAcc<-meanPredAccuracy(modelType = modelType,
                                     crossValOut = cvPredMeans,
                                     snpeffs = markEffs,
                                     ped = ped,selInd = selInd,SIwts = SIwts)

  if(!is.null(outName)){
    print("Saving outputs to disk.")
    saveRDS(parentfolds,file=paste0(outName,"_parentfolds.rds"))
    saveRDS(markEffs,file=paste0(outName,"_markerEffects.rds"))
    saveRDS(cvPredVars,file=paste0(outName,"_predVars.rds"))
    saveRDS(cvPredMeans,file=paste0(outName,"_predMeans.rds"))
    saveRDS(varPredAcc,file=paste0(outName,"_varPredAccuracy.rds"))
    saveRDS(meanPredAcc,file=paste0(outName,"_meanPredAccuracy.rds"))
  }

  accuracy_out<-list(meanPredAccuracy=meanPredAcc,
                     varPredAccuracy=varPredAcc)
  print(paste0("Accuracies predicted. Took  ",
               round((proc.time()[3] - initime)/60/60,5),
               " hrs total.Goodbye!"))
  return(accuracy_out)
}

getMarkEffs<-function(parentfolds,blups,gid,modelType,grms,dosages,ncores,nBLASthreads=NULL){

  traintestdata<-parentfolds %>%
    dplyr::select(Repeat,Fold,trainset,testset) %>%
    pivot_longer(c(trainset,testset),
                 names_to = "Dataset",
                 values_to = "sampleIDs") %>%
    crossing(Trait=blups$Trait) %>%
    left_join(blups) %>%
    rename(blupsMat=blups)

  fitModel<-function(sampleIDs,blupsMat,modelType,gid,grms,dosages,nBLASthreads,...){
    # debug
    # sampleIDs<-traintestdata$sampleIDs[[2]]; blupsMat<-traintestdata$blupsMat[[2]]

    if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
    # workers in plan(multisession) need this call internal to the function, it seems.

    A<-grms[["A"]]
    if(modelType %in% c("AD","DirDom")){ D<-grms[["D"]] }

    trainingdata<-blupsMat %>%
      dplyr::rename(gid=!!sym(gid)) %>%
      filter(gid %in% sampleIDs)

    trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[["gid"]],
                                            levels=rownames(A))
    if(modelType %in% c("AD")){
      trainingdata[[paste0(gid,"d")]]<-trainingdata[[paste0(gid,"a")]] }
    if(modelType %in% c("DirDom")){
      trainingdata[[paste0(gid,"d_star")]]<-trainingdata[[paste0(gid,"a")]] }

    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A)")
    if(modelType %in% c("AD")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D)") }
    if(modelType=="DirDom"){
      randFormula<-paste0(randFormula,"+vs(",gid,"d_star,Gu=D)")
      f<-getPropHom(dosages)
      trainingdata %<>% mutate(f=f[trainingdata$gid]) }

    # Fixed model statements
    fixedFomula<-ifelse(modelType=="DirDom",
                        "drgBLUP ~1+f","drgBLUP ~1")
    # Fit genomic prediction model
    require(sommer)
    fit <- sommer::mmer(fixed = as.formula(fixedFomula),
                        random = as.formula(randFormula),
                        weights = WT,
                        data=trainingdata,
                        date.warning = F)

    # Backsolve SNP effects
    # Compute allele sub effects
    ## Every model has an additive random term
    ga<-as.matrix(fit$U[[paste0("u:",gid,"a")]]$drgBLUP,ncol=1)
    M<-centerDosage(dosages)

    if(modelType %in% c("A","AD")){
      # models A and AD give add effects corresponding to allele sub effects
      allelesubsnpeff<-backsolveSNPeff(Z=M,g=ga) }

    if(modelType %in% c("AD")){
      # model AD the dom effects are dominance deviation effects
      gd<-as.matrix(fit$U[[paste0("u:",gid,"d")]]$drgBLUP,ncol=1)
      domdevsnpeff<-backsolveSNPeff(Z=dose2domDev(dosages),g=gd) }

    if(modelType %in% c("DirDom")){
      # model DirDom is a different add-dom partition,
      ### add effects are not allele sub effects and gblups are not GEBV
      addsnpeff<-backsolveSNPeff(Z=M,g=ga)
      ### dom effects are called d*, gd_star or domstar
      ### because of the genome-wide homoz. term included in model
      gd_star<-as.matrix(fit$U[[paste0("u:",gid,"d_star")]]$drgBLUP,ncol=1)
      domdevMat_genotypic<-dose2domDevGenotypic(dosages)
      domstar_snpeff<-backsolveSNPeff(Z=domdevMat_genotypic,g=gd_star)
      ### b = the estimate (BLUE) for the genome-wide homoz. effect
      b<-fit$Beta[fit$Beta$Effect=="f","Estimate"]
      ### calc. domsnpeff including the genome-wide homoz. effect
      ### divide the b effect up by number of SNPs and _subtract_ from domstar
      domsnpeff<-domstar_snpeff-(b/length(domstar_snpeff))

      ### allele substitution effects using a+d(q-p) where d=d*-b/p
      p<-getAF(dosages)
      q<-1-p
      allelesubsnpeff<-addsnpeff+(domsnpeff*(q-p))
    }

    # Gather the GBLUPs
    if(modelType %in% c("A","AD")){
      gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                     GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)) }
    if(modelType=="AD"){
      gblups %<>% # compute GEDD (genomic-estimated dominance deviation)
        mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP),
               # compute GETGV
               GETGV=rowSums(.[,grepl("GE",colnames(.))])) }
    if(modelType=="DirDom"){
      # re-calc the GBLUP, GEdomval using dom. effects where d=d*-b/p
      ge_domval<-domdevMat_genotypic%*%domsnpeff
      # calc. the GEBV using allele sub. effects where alpha=a+d(p-q), and d=d*-b/p
      gebv<-M%*%allelesubsnpeff
      # Tidy tibble of GBLUPs
      gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                     GEadd=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP),
                     GEdom_star=as.numeric(fit$U[[paste0("u:",gid,"d_star")]]$drgBLUP)) %>%
        left_join(tibble(GID=rownames(ge_domval),GEdomval=as.numeric(ge_domval))) %>%
        left_join(tibble(GID=rownames(gebv),GEBV=as.numeric(gebv))) %>%
        # GETGV from GEadd + GEdomval
        mutate(GETGV=GEadd+GEdomval)
    }

    # Extract variance components
    varcomps<-summary(fit)$varcomp

    # Exract fixed effects
    # for modelType="DirDom", contains estimate of genome-wide homoz. effect
    fixeffs<-summary(fit)$betas

    results<-tibble(gblups=list(gblups),
                    varcomps=list(varcomps),
                    fixeffs=list(fixeffs))
    # Add snpeffects to output
    results %<>% mutate(allelesubsnpeff=list(allelesubsnpeff))
    if(modelType=="AD"){ results %<>% mutate(domdevsnpeff=list(domdevsnpeff)) }
    if(modelType=="DirDom"){
      results %<>% mutate(addsnpeff=list(addsnpeff),
                          domstar_snpeff=list(domstar_snpeff),
                          domsnpeff=list(domsnpeff)) }
    # this is to remove conflicts with dplyr function select() downstream
    detach("package:sommer",unload = T); detach("package:MASS",unload = T)

    # return results
    return(results)
  }

  require(furrr); require(future.callr); plan(callr, workers = ncores, gc=TRUE)
  #require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  traintestdata<-traintestdata %>%
    mutate(modelOut=future_pmap(.,fitModel,
                                modelType=modelType,
                                gid=gid,
                                grms=grms,
                                dosages=dosages,
                                nBLASthreads=nBLASthreads),
           modelType=modelType)

  traintestdata %<>%
    select(-blupsMat,-sampleIDs) %>%
    unnest(modelOut,keep_empty = FALSE) %>% # dropped failed models on unnest
    nest(effects=c(-Repeat,-Fold,-Dataset,-modelType))

  return(traintestdata)
}

predictCrossVars<-function(modelType,snpeffs,parentfolds,
                           haploMat,recombFreqMat,ncores){
  predvars<-snpeffs %>%
    unnest(effects,keep_empty = FALSE) %>%
    filter(Dataset=="trainset") %>%
    dplyr::select(Repeat,Fold,Trait,modelType,contains("snpeff")) %>%
    nest(EffectList=c(Trait,contains("snpeff"))) %>%
    # AlleleSubEffects to predict VarBV for all models
    mutate(AlleleSubEffectList=map(EffectList,
                                   function(EffectList){
                                     allelesubsnpeff<-map(EffectList$allelesubsnpeff,~t(.))
                                     names(allelesubsnpeff)<-EffectList$Trait
                                     return(allelesubsnpeff)}))
  if(modelType=="AD"){
    predvars<-predvars %>%
      # DomDevEffects for model "AD" to predict VarTGV = VarBV + VarDD
      mutate(DomDevEffectList=map(EffectList,
                                  function(EffectList){
                                    domdevsnpeff<-map(EffectList$domdevsnpeff,~t(.))
                                    names(domdevsnpeff)<-EffectList$Trait
                                    return(domdevsnpeff) })) }
  if(modelType=="DirDom"){
    predvars<-predvars %>%
      # AddEffectList + DomEffectList --> VarTGV; AlleleSubEffectList --> VarBV;
      mutate(AddEffectList=map(EffectList,
                               function(EffectList){
                                 addsnpeff<-map(EffectList$addsnpeff,~t(.))
                                 names(addsnpeff)<-EffectList$Trait
                                 return(addsnpeff) }),
             DomEffectList=map(EffectList,
                               function(EffectList){
                                 domsnpeff<-map(EffectList$domsnpeff,~t(.))
                                 names(domsnpeff)<-EffectList$Trait
                                 return(domsnpeff) })) }

  predvars %<>%
    left_join(parentfolds %>%
                dplyr::select(-testparents,-trainset,-testset)) %>%
    dplyr::select(-EffectList)

  require(furrr);
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")


  if(modelType=="A"){
    predvars<-predvars %>%
      mutate(predVars=map2(CrossesToPredict,AlleleSubEffectList,
                           ~predCrossVars(CrossesToPredict=.x,
                                          AddEffectList=.y,
                                          modelType="A",
                                          haploMat=haploMat,
                                          recombFreqMat=recombFreqMat,
                                          ncores=ncores) %>%
                             unnest(predVars) %>%
                             mutate(predOf="VarBV") %>%
                             nest(predVars=c(Trait1,Trait2,predOf,predVar)))) }
  if(modelType=="AD"){
    predvars<-predvars %>%
      mutate(predVars=pmap(.,function(CrossesToPredict,
                                      AlleleSubEffectList,
                                      DomDevEffectList,...){
        out<-predCrossVars(CrossesToPredict=CrossesToPredict,
                           AddEffectList=AlleleSubEffectList,
                           DomEffectList=DomDevEffectList,
                           modelType="AD",
                           haploMat=haploMat,
                           recombFreqMat=recombFreqMat,
                           ncores=ncores)
        out<-out %>%
          unnest(predVars) %>%
          mutate(predOf=ifelse(predOf=="VarA","VarBV","VarDD")) %>%
          nest(predVars=c(Trait1,Trait2,predOf,predVar))
        return(out) })) }
  if(modelType=="DirDom"){
    predvars<-predvars %>%
      mutate(predVars=pmap(.,function(CrossesToPredict,
                                      AlleleSubEffectList,
                                      AddEffectList,
                                      DomEffectList,...){

        predVarTGV<-predCrossVars(CrossesToPredict=CrossesToPredict,
                                  AddEffectList=AddEffectList,
                                  DomEffectList=DomEffectList,
                                  modelType="AD", # no "DirDom" model in predCrossVars() nor is it needed
                                  haploMat=haploMat,
                                  recombFreqMat=recombFreqMat,
                                  ncores=ncores)
        predVarBV<-predCrossVars(CrossesToPredict=CrossesToPredict,
                                 AddEffectList=AlleleSubEffectList,
                                 DomEffectList=NULL,
                                 modelType="A", # no "DirDom" model in predCrossVars() nor is it needed
                                 haploMat=haploMat,
                                 recombFreqMat=recombFreqMat,
                                 ncores=ncores)
        out<-predVarBV %>%
          unnest(predVars) %>%
          mutate(predOf="VarBV") %>%
          bind_rows(predVarTGV %>%
                      unnest(predVars)) %>%
          nest(predVars=c(Trait1,Trait2,predOf,predVar))
        return(out) })) }
  predvars %<>% select(-contains("EffectList"),-CrossesToPredict)
  return(predvars)
}

predictCrossMeans<-function(modelType,snpeffs,parentfolds,
                            doseMat,ncores){
  predmeans<-snpeffs %>%
    unnest(effects,keep_empty = FALSE) %>%
    filter(Dataset=="trainset") %>%
    dplyr::select(Repeat,Fold,Trait,modelType,contains("snpeff")) %>%
    nest(EffectList=c(Trait,contains("snpeff"))) %>%
    # AlleleSubEffects are available from all prediction models
    ## therefore, prediction of MeanBV is available for all models
    mutate(AlleleSubEffectList=map(EffectList,
                                   function(EffectList){
                                     allelesubsnpeff<-map(EffectList$allelesubsnpeff,~t(.))
                                     names(allelesubsnpeff)<-EffectList$Trait
                                     return(allelesubsnpeff)})) %>%
    left_join(parentfolds %>%
                dplyr::select(-testparents,-trainset,-testset))

  ## predict MeanBVs
  predmeans %<>%
    mutate(predMeans=pmap(.,function(AlleleSubEffectList,CrossesToPredict,...){
      return(predCrossMeans(AddEffectList=AlleleSubEffectList,
                            CrossesToPredict=CrossesToPredict,
                            doseMat=doseMat,ncores=ncores,predType="BV")) }))
  ## predict MeanTGVs
  if(modelType=="DirDom"){
    #  Prediction of MeanTGV is only available for the DirDom model
    ### or a model with "genotypic" additive-dominance SNP effects
    ### As implemented, modelType="AD" is the "classical" partition (BVs+ DomDevs)
    predmeans<-predmeans %>%
      # AddEffectList + DomEffectList --> VarTGV; AlleleSubEffectList --> VarBV;
      mutate(AddEffectList=map(EffectList,
                               function(EffectList){
                                 addsnpeff<-map(EffectList$addsnpeff,~t(.))
                                 names(addsnpeff)<-EffectList$Trait
                                 return(addsnpeff) }),
             DomEffectList=map(EffectList,
                               function(EffectList){
                                 domsnpeff<-map(EffectList$domsnpeff,~t(.))
                                 names(domsnpeff)<-EffectList$Trait
                                 return(domsnpeff) }))

    ## prediction of MeanTGVs
    predmeans %<>%
      mutate(predMeans=pmap(.,function(predMeans,AddEffectList,DomEffectList,
                                       CrossesToPredict,...){
        return(predMeans %>%
                 bind_rows(predCrossMeans(AddEffectList=AddEffectList,
                                          DomEffectList=DomEffectList,
                                          CrossesToPredict=CrossesToPredict,
                                          doseMat=doseMat,ncores=ncores,
                                          predType="TGV"))) }))

  }
  predmeans %<>%
    select(-contains("EffectList"),-CrossesToPredict)
  return(predmeans)
}

varPredAccuracy<-function(crossValOut,snpeffs,ped,modelType,
                          selInd=FALSE,SIwts=NULL){

  # Extract and format the GBLUPs from the marker effects object
  gblups<-snpeffs %>%
    unnest(effects,keep_empty = FALSE) %>%
    filter(Dataset=="testset") %>%
    select(Repeat,Fold,modelType,Trait,gblups) %>%
    unnest(gblups) %>%
    nest(testset_gblups=c(-Repeat,-Fold,-modelType))

  # Use the crossValPred object and the pedigree
  # Create a list of the actual members of each family that were predicted
  # in each repeat-fold
  # Join the GBLUPs for each family member for computing
  # cross sample means, variances, covariances
  out<-crossValOut %>%
    unnest(predVars) %>%
    distinct(Repeat,Fold,modelType,sireID,damID) %>%
    left_join(ped) %>%
    nest(CrossesToPredict=c(sireID,damID,GID)) %>%
    left_join(gblups)
  out %<>%
    # remove any gebv/getgv NOT in the crosses-to-be-predicted to save mem
    mutate(testset_gblups=map2(testset_gblups,CrossesToPredict,
                               ~semi_join(.x,.y)))
  # for modelType=="A" remove the GETGV as equiv. to GEBV
  if(modelType=="A"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,Trait,GID,GEBV) %>%
                                  pivot_longer(.,cols = c(GEBV),
                                               names_to = "predOf",
                                               values_to = "GBLUP") %>%
                                  nest(gblups=-predOf))) }
  # for modelType=="AD" remove the GEDD, pivot to long form GEBV/GETGV
  if(modelType=="AD"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,Trait,GID,GEBV,GETGV) %>%
                                  pivot_longer(cols = c(GEBV,GETGV),
                                               names_to = "predOf",
                                               values_to = "GBLUP") %>%
                                  nest(gblups=-predOf)))
  }
  if(modelType=="DirDom"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,Trait,GID,GEBV,GETGV) %>%
                                  pivot_longer(cols = c(GEBV,GETGV),
                                               names_to = "predOf",
                                               values_to = "GBLUP") %>%
                                  nest(gblups=-predOf)))
  }
  out %<>% unnest(testset_gblups)

  # make a matrix of GBLUPs for all traits
  # for each family-to-be-predicted
  # in each rep-fold-predOf combination
  out %<>%
    mutate(famgblups=map2(gblups,CrossesToPredict,
                          ~left_join(.x,.y) %>%
                            pivot_wider(names_from = "Trait",
                                        values_from = "GBLUP") %>%
                            nest(gblupmat=c(-sireID,-damID)) %>%
                            mutate(gblupmat=map(gblupmat,~column_to_rownames(.,var="GID"))))) %>%
    select(-CrossesToPredict,-gblups) %>%
    unnest(famgblups)

  #gblupmat<-out$gblupmat[[1]]
  out %<>%
    # outer loop over rep-fold-predtype
    mutate(obsVars=map(gblupmat,function(gblupmat){
      #gblupmat<-famgblups$gblupmat[[1]]
      covMat<-cov(gblupmat)
      # to match predCrossVar output
      ## keep upper tri + diag of covMat
      obsvars<-covMat
      obsvars[lower.tri(obsvars)]<-NA
      obsvars %<>%
        as.data.frame(.) %>%
        rownames_to_column(var = "Trait1") %>%
        pivot_longer(cols = c(-Trait1),
                     names_to = "Trait2",
                     values_to = "obsVar",
                     values_drop_na = T)
      if(selInd==TRUE){
        covmat<-covMat[names(SIwts),names(SIwts)]
        selIndVar<-SIwts%*%covmat%*%SIwts
        obsvars %<>%
          bind_rows(tibble(Trait1="SELIND",
                           Trait2="SELIND",
                           obsVar=selIndVar),.) }
      obsvars %<>% mutate(obsVar=as.numeric(obsVar))
      return(obsvars) }),
      famSize=map_dbl(gblupmat,nrow)) %>%
    select(-gblupmat) %>%
    unnest(obsVars)

  cvout<-crossValOut %>%
    unnest(predVars) %>%
    unnest(predVars) %>%
    select(Repeat,Fold,modelType,predOf,sireID,damID,Trait1,Trait2,predVar)

  if(modelType=="AD"){
    ## predVarTGV = predVarBV + predVarDD
    cvout %<>%
      filter(predOf=="VarBV") %>%
      bind_rows(cvout %>%
                  group_by(Repeat,Fold,modelType,sireID,damID,Trait1,Trait2) %>%
                  summarize(predVar=sum(predVar),.groups = 'drop') %>%
                  mutate(predOf="VarTGV"))
  }

  if(modelType=="DirDom"){
    ## predVarTGV = predVarA + predVarD
    cvout %<>%
      filter(predOf=="VarBV") %>%
      bind_rows(cvout %>%
                  filter(predOf %in% c("VarA","VarD")) %>%
                  group_by(Repeat,Fold,modelType,sireID,damID,Trait1,Trait2) %>%
                  summarize(predVar=sum(predVar),.groups = 'drop') %>%
                  mutate(predOf="VarTGV"))
  }

  cvout %<>%
    nest(predVars=c(Trait1,Trait2,predVar))

  if(selInd==TRUE){
    # compute predicted selection index variances

    require(furrr); require(future.callr); plan(callr, workers = ncores, gc=TRUE)
    #require(furrr); plan(multisession, workers = ncores)
    options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

    cvout %<>%
      ## loop over each rep-fold-predOf-sireIDxdamID
      mutate(predVars=future_map(predVars,function(predVars){
        gmat<-predVars %>%
          pivot_wider(names_from = "Trait2",
                      values_from = "predVar") %>%
          column_to_rownames(var = "Trait1") %>%
          as.matrix
        gmat[lower.tri(gmat)]<-t(gmat)[lower.tri(gmat)]
        gmat %<>% .[names(SIwts),names(SIwts)]
        predSelIndVar<-SIwts%*%gmat%*%SIwts
        ## add sel index predictions to component trait
        ## var-covar predictions

        predVars<-tibble(Trait1="SELIND",Trait2="SELIND",
                         predVar=as.numeric(predSelIndVar)) %>%
          bind_rows(predVars)
        return(predVars) }))
  }
  out %<>%
    mutate(predOf=ifelse(predOf=="GEBV","VarBV","VarTGV")) %>%
    left_join(cvout %>%
                unnest(predVars)) %>%
    nest(predVSobs=c(sireID,damID,predVar,obsVar,famSize)) %>%
    mutate(AccuracyEst=map_dbl(predVSobs,function(predVSobs){
      out<-psych::cor.wt(predVSobs[,c("predVar","obsVar")],
                         w = predVSobs$famSize) %$% r[1,2] %>%
        round(.,3)
      return(out) }))
  return(out)
}


meanPredAccuracy<-function(crossValOut,snpeffs,ped,modelType,
                           selInd=FALSE,SIwts=NULL){
  # Extract and format the GBLUPs from the marker effects object
  gblups<-snpeffs %>%
    unnest(effects,keep_empty = FALSE) %>%
    filter(Dataset=="testset") %>%
    select(Repeat,Fold,modelType,Trait,gblups) %>%
    unnest(gblups) %>%
    nest(testset_gblups=c(-Repeat,-Fold,-modelType))

  # Use the crossValPred object and the pedigree
  # Create a list of the actual members of each family that were predicted
  # in each repeat-fold
  # Join the GBLUPs for each family member for computing
  # cross sample means
  out<-crossValOut %>%
    unnest(predMeans) %>%
    distinct(Repeat,Fold,modelType,sireID,damID) %>%
    left_join(ped) %>%
    nest(CrossesToPredict=c(sireID,damID,GID)) %>%
    left_join(gblups)

  out %<>%
    # remove any gebv/getgv NOT in the crosses-to-be-predicted to save mem
    mutate(testset_gblups=map2(testset_gblups,CrossesToPredict,
                               ~semi_join(.x,.y)))

  # for modelType=="A" and "AD", MeanBV only (remove GETGV)
  if(modelType %in% c("A","AD")){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,Trait,GID,GEBV) %>%
                                  pivot_longer(.,cols = c(GEBV),
                                               names_to = "predOf",
                                               values_to = "GBLUP") %>%
                                  nest(gblups=-predOf))) }

  if(modelType=="DirDom"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,Trait,GID,GEBV,GETGV) %>%
                                  pivot_longer(cols = c(GEBV,GETGV),
                                               names_to = "predOf",
                                               values_to = "GBLUP") %>%
                                  nest(gblups=-predOf)))
  }
  out %<>% unnest(testset_gblups)

  # make a matrix of GBLUPs for all traits
  # for each family-to-be-predicted
  # in each rep-fold-predOf combination
  out %<>%
    mutate(famgblups=map2(gblups,CrossesToPredict,
                          ~left_join(.x,.y) %>%
                            pivot_wider(names_from = "Trait",
                                        values_from = "GBLUP") %>%
                            nest(gblupmat=c(-sireID,-damID)) %>%
                            mutate(gblupmat=map(gblupmat,~column_to_rownames(.,var="GID"))))) %>%
    select(-CrossesToPredict,-gblups) %>%
    unnest(famgblups)

  out %<>%
    # outer loop over rep-fold-predtype
    mutate(obsMeans=map(gblupmat,function(gblupmat){
      #      gblupmeans<-colMeans(gblupmat) %>% as.list
      gblupmeans<-colMeans(gblupmat) %>% as.data.frame %>% as.matrix

      if(selInd==TRUE){
        selIndMean<-gblupmeans[names(SIwts),]%*%SIwts %>%
          `row.names<-`(.,"SELIND")
        gblupmeans<-rbind(selIndMean,gblupmeans)
      }
      obsmeans<-tibble(Trait=rownames(gblupmeans),
                       obsMean=as.numeric(gblupmeans))
      return(obsmeans) }),
      famSize=map_dbl(gblupmat,nrow)) %>%
    select(-gblupmat) %>%
    unnest(obsMeans)


  cvout<-crossValOut %>%
    unnest(predMeans) %>%
    select(Repeat,Fold,modelType,predOf,sireID,damID,Trait,predMean)

  if(selInd==TRUE){
    # compute predicted selection index variances
    cvout %<>%
      pivot_wider(names_from = "Trait",
                  values_from = "predMean") %>%
      mutate(SELIND=as.numeric(cvout %>%
                                 pivot_wider(names_from = "Trait",
                                             values_from = "predMean") %>%
                                 select(any_of(names(SIwts))) %>%
                                 as.matrix(.)%*%SIwts)) %>%
      pivot_longer(cols = c("SELIND",unique(cvout$Trait)),
                   names_to = "Trait",
                   values_to = "predMean")
  }

  out %<>%
    mutate(predOf=ifelse(predOf=="GEBV","MeanBV","MeanTGV")) %>%
    left_join(cvout) %>%
    nest(predVSobs=c(sireID,damID,predMean,obsMean,famSize)) %>%
    mutate(AccuracyEst=map_dbl(predVSobs,function(predVSobs){
      out<-psych::cor.wt(predVSobs[,c("predMean","obsMean")],
                         w = predVSobs$famSize) %$% r[1,2] %>%
        round(.,3)
      return(out) }))
  return(out)
}


# Prunes out offspring, grandkids, greatgrandkids (up to X4) steps of
# great ancestors.  It is not automatically recursive across any depth of
# pedigree. That depth works for current test pedigree (IITA 2021).
# Must name parent columns in ped "sireID" and "damID".
makeParentFolds<-function(ped,gid,nrepeats=5,nfolds=5,seed=NULL){
  require(rsample)
  set.seed(seed)
  parentfolds<-rsample::vfold_cv(tibble(Parents=union(ped$sireID,
                                                      ped$damID)),
                                 v = nfolds,repeats = nrepeats) %>%
    mutate(folds=map(splits,function(splits){
      #splits<-parentfolds$splits[[1]]
      testparents<-testing(splits)$Parents
      trainparents<-training(splits)$Parents
      ped<-ped %>%
        rename(gid=!!sym(gid))
      offspring<-ped %>%
        filter(sireID %in% testparents | damID %in% testparents) %$%
        unique(gid)
      grandkids<-ped %>%
        filter(sireID %in% offspring | damID %in% offspring) %$%
        unique(gid)
      greatX1grandkids<-ped %>%
        filter(sireID %in% grandkids | damID %in% grandkids) %$%
        unique(gid)
      greatX2grandkids<-ped %>%
        filter(sireID %in% greatX1grandkids |
                 damID %in% greatX1grandkids) %$%
        unique(gid)
      greatX3grandkids<-ped %>%
        filter(sireID %in% greatX2grandkids |
                 damID %in% greatX2grandkids) %$%
        unique(gid)
      greatX4grandkids<-ped %>%
        filter(sireID %in% greatX3grandkids |
                 damID %in% greatX3grandkids) %$%
        unique(gid)

      testset<-unique(c(offspring,
                        grandkids,
                        greatX1grandkids,
                        greatX2grandkids,
                        greatX3grandkids,
                        greatX4grandkids)) %>%
        .[!. %in% c(testparents,trainparents)]

      nontestdescendents<-ped %>%
        filter(!gid %in% testset) %$%
        unique(gid)
      trainset<-union(testparents,trainparents) %>%
        union(.,nontestdescendents)

      out<-tibble(testparents=list(testparents),
                  trainset=list(trainset),
                  testset=list(testset))
      return(out) })) %>%
    unnest(folds)
  if(nrepeats>1){
    parentfolds %<>%
      rename(Repeat=id,Fold=id2) %>%
      select(-splits)
  }
  if(nrepeats==1){
    parentfolds %<>%
      mutate(Repeat="Repeat1") %>%
      rename(Fold=id) %>%
      select(-splits)
    }


  # Crosses To Predict
  parentfolds %<>%
    mutate(CrossesToPredict=map(testparents,
                                ~filter(ped %>%
                                          # only need a list of fams-to-predict
                                          # not the progeny info
                                          distinct(damID,sireID),
                                        sireID %in% . | damID %in% .)))
  return(parentfolds)
}

#nfolds=5; nrepeats=1; gid="GID"; seed=53;
