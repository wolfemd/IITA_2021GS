runParentWiseCrossVal<-function(nrepeats,nfolds,seed=NULL,modelType,
                                     ncores=1,outName=NULL,
                                     ped=ped,gid="GID",blups,
                                     dosages,grms,haploMat,recombFreqMat,
                                     selInd,SIwts = NULL,...){
  # test version will subset to and only analyze one rep-fold combo for speed
  initime<-proc.time()[3]

  ## Make parent-wise folds
  parentfolds<-makeParentFolds(ped=ped,gid="GID",
                               nrepeats=nrepeats,
                               nfolds=nfolds,
                               seed=seed)
  print("Set-up parent-wise cross-validation folds")

  ## Get univariate REML marker effects
  #### modelType=="AD"
  print("Fitting models to get marker effects")
  starttime<-proc.time()[3]
  markEffs<-getMarkEffs(parentfolds,blups=blups,gid=gid,modelType=modelType,
                        grms=grms,dosages=dosages,ncores=ncores)
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
  varPredAccuracy<-varPredAccuracy(modelType = modelType,
                                   crossValOut = cvPredVars,
                                   snpeffs = markEffs,
                                   ped = ped,selInd = selInd,SIwts = SIwts)

  ## Mean prediction accuracies
  meanPredAccuracy<-meanPredAccuracy(modelType = modelType,
                                     crossValOut = cvPredMeans,
                                     snpeffs = markEffs,
                                     ped = ped,selInd = selInd,SIwts = SIwts)

  if(!is.null(outName)){
    print("Saving outputs to disk.")
    saveRDS(parentfolds,file=paste0(outName,"_parentfolds.rds"))
    saveRDS(markEffs,file=paste0(outName,"_markerEffects.rds"))
    saveRDS(cvPredVars,file=paste0(outName,"_predVars.rds"))
    saveRDS(cvPredMeans,file=paste0(outName,"_predMeans.rds"))
    saveRDS(varPredAccuracy,file=paste0(outName,"_varPredAccuracy.rds"))
    saveRDS(meanPredAccuracy,file=paste0(outName,"_meanPredAccuracy.rds"))
  }

  accuracy_out<-list(meanPredAccuracy=meanPredAccuracy,
                     varPredAccuracy=varPredAccuracy)
  print(paste0("Accuracies predicted. Took  ",
               round((proc.time()[3] - initime)/60/60,5),
               " hrs total.\n Goodbye!"))
  return(accuracy_out)
}


getMarkEffs<-function(parentfolds,blups,gid,modelType,grms,dosages,ncores){
  traintestdata<-parentfolds %>%
    dplyr::select(Repeat,Fold,trainset,testset) %>%
    pivot_longer(c(trainset,testset),
                 names_to = "Dataset",
                 values_to = "sampleIDs") %>%
    crossing(Trait=blups$Trait) %>%
    left_join(blups) %>%
    rename(blupsMat=blups)

  ## For each training/testing chunk of sampleIDs and each trait
  ## fit GBLUP model and backsolve SNP-effects

  fitModel<-function(sampleIDs,blupsMat,modelType,gid,grms,dosages,...){
    # debug
    # sampleIDs<-traintestdata$sampleIDs[[2]]; blups<-traintestdata$blups[[2]]
    require(predCrossVar)
    A<-grms[["A"]]
    if(modelType %in% c("AD")){ D<-grms[["D"]] }
    trainingdata<-blupsMat %>%
      dplyr::rename(gid=!!sym(gid)) %>%
      filter(gid %in% sampleIDs)
    trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[["gid"]],
                                            levels=rownames(A))
    if(modelType %in% c("AD")){
      trainingdata[[paste0(gid,"d")]]<-trainingdata[[paste0(gid,"a")]]
    }
    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A)")
    if(modelType %in% c("AD")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D)")
    }

    # Fit genomic prediction model
    require(sommer)
    fit <- sommer::mmer(fixed = drgBLUP ~1,
                        random = as.formula(randFormula),
                        weights = WT,
                        data=trainingdata)

    # Gather the GBLUPs
    gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                   GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP))
    if(modelType %in% c("AD")){
      gblups %<>%
        mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP))
    }
    # Calc GETGVs
    ## Note that for modelType=="A", GEBV==GETGV
    gblups %<>%
      mutate(GETGV=rowSums(.[,grepl("GE",colnames(.))]))

    # Backsolve SNP effects
    ga<-as.matrix(fit$U[[paste0("u:",gid,"a")]]$drgBLUP,ncol=1)
    addsnpeff<-backsolveSNPeff(Z=centerDosage(dosages),g=ga)
    if(modelType %in% c("AD")){
      gd<-as.matrix(fit$U[[paste0("u:",gid,"d")]]$drgBLUP,ncol=1)
      domsnpeff<-backsolveSNPeff(Z=dose2domDev(dosages),g=gd)
    }

    # Extract variance components
    varcomps<-summary(fit)$varcomp

    results<-tibble(gblups=list(gblups),
                    varcomps=list(varcomps),
                    addsnpeff=list(addsnpeff))
    if(modelType %in% c("AD")){
      results %<>%
        mutate(domsnpeff=list(domsnpeff)) }
    # return results
    return(results)
  }

  require(furrr); options(mc.cores=ncores); plan(multicore)
  options(future.globals.maxSize=50000*1024^2)
  options(future.rng.onMisuse="ignore")
  traintestdata<-traintestdata %>%
    mutate(modelOut=future_pmap(.,fitModel,
                                modelType=modelType,
                                gid=gid,
                                grms=grms,
                                dosages=dosages),
           modelType=modelType)

  traintestdata %<>%
    select(-blupsMat,-sampleIDs) %>%
    unnest(modelOut) %>%
    nest(effects=c(Trait,gblups,varcomps,addsnpeff,domsnpeff))

  # this is to remove conflicts with dplyr function select() downstream
  # detach("package:sommer",unload = T); detach("package:MASS",unload = T)

  return(traintestdata)
}
predictCrossVars<-function(modelType,snpeffs,parentfolds,
                           haploMat,recombFreqMat,ncores){
  predvars<-snpeffs %>%
    unnest(effects) %>%
    filter(Dataset=="trainset") %>%
    dplyr::select(Repeat,Fold,Trait,modelType,
                  any_of(c("addsnpeff","domsnpeff"))) %>%
    nest(EffectList=c(Trait,any_of(c("addsnpeff","domsnpeff")))) %>%
    mutate(AddEffectList=map(EffectList,
                             function(EffectList){
                               addsnpeff<-map(EffectList$addsnpeff,~t(.))
                               names(addsnpeff)<-EffectList$Trait
                               return(addsnpeff)}))
  if(modelType %in% c("AD")){
    predvars<-predvars %>%
      mutate(DomEffectList=map(EffectList,
                               function(EffectList){
                                 domsnpeff<-map(EffectList$domsnpeff,~t(.))
                                 names(domsnpeff)<-EffectList$Trait
                                 return(domsnpeff) })) }

  predvars %<>%
    left_join(parentfolds %>%
                dplyr::select(-testparents,-trainset,-testset)) %>%
    dplyr::select(-EffectList)

  require(furrr); options(future.globals.maxSize=40000*1024^2)

  if(modelType=="A"){
    predvars<-predvars %>%
      mutate(predVars=map2(CrossesToPredict,AddEffectList,
                           ~predCrossVars(CrossesToPredict=.x,
                                          AddEffectList=.y,
                                          modelType=modelType,
                                          haploMat=haploMat,
                                          recombFreqMat=recombFreqMat,
                                          ncores=ncores))) }
  if(modelType=="AD"){
    predvars<-predvars %>%
      mutate(predVars=pmap(.,function(CrossesToPredict,
                                      AddEffectList,DomEffectList,...){
        out<-predCrossVars(CrossesToPredict=CrossesToPredict,
                           AddEffectList=AddEffectList,
                           DomEffectList=DomEffectList,
                           modelType=modelType,
                           haploMat=haploMat,
                           recombFreqMat=recombFreqMat,
                           ncores=ncores)
        return(out) })) }
  predvars %<>% select(-AddEffectList,-DomEffectList,-CrossesToPredict)
  return(predvars)
}
predictCrossMeans<-function(modelType,snpeffs,parentfolds,
                            doseMat,ncores){
  predmeans<-snpeffs %>%
    unnest(effects) %>%
    filter(Dataset=="trainset") %>%
    dplyr::select(Repeat,Fold,Trait,modelType,
                  any_of(c("addsnpeff","domsnpeff"))) %>%
    nest(EffectList=c(Trait,any_of(c("addsnpeff","domsnpeff")))) %>%
    mutate(AddEffectList=map(EffectList,
                             function(EffectList){
                               addsnpeff<-map(EffectList$addsnpeff,~t(.))
                               names(addsnpeff)<-EffectList$Trait
                               return(addsnpeff)}))

  if(modelType %in% c("AD")){
    predmeans<-predmeans %>%
      mutate(DomEffectList=map(EffectList,
                               function(EffectList){
                                 domsnpeff<-map(EffectList$domsnpeff,~t(.))
                                 names(domsnpeff)<-EffectList$Trait
                                 return(domsnpeff) })) }

  predmeans %<>%
    left_join(parentfolds %>%
                dplyr::select(-testparents,-trainset,-testset)) %>%
    dplyr::select(-EffectList)

  predmeans %<>%
    mutate(predMeans=pmap(.,predCrossMeans,doseMat=doseMat,ncores=ncores)) %>%
    select(-contains("EffectList"),-CrossesToPredict)
  return(predmeans)
}
varPredAccuracy<-function(crossValOut,snpeffs,ped,modelType,
                          selInd=FALSE,SIwts=NULL){

  # Extract and format the GBLUPs from the marker effects object
  gblups<-snpeffs %>%
    unnest(effects) %>%
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
    select(Repeat,Fold,modelType,sireID,damID) %>%
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
                                ~pivot_longer(.,cols = c(GEBV,GETGV),
                                              names_to = "predOf",
                                              values_to = "GBLUP") %>%
                                  nest(gblups=-predOf) %>%
                                  filter(predOf=="GEBV")))
  }
  # for modelType=="AD" remove the GEDD, pivot to long form GEBV/GETGV
  if(modelType=="AD"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,-GEDD) %>%
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
    select(-CrossesToPredict,-gblups)

  #famgblups<-out$famgblups[[1]]
  out %<>%
    # outer loop over rep-fold-predtype
    mutate(obsVars=map(famgblups,function(famgblups){
      return(famgblups %>%
               # inner loop over families
               mutate(obsvars=map(gblupmat,
                                  function(gblupmat){
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
                                    return(obsvars) }),
                      famSize=map_dbl(gblupmat,nrow)) %>%
               select(-gblupmat) %>%
               unnest(obsvars))})) %>%
    select(-famgblups)

  cvout<-crossValOut %>%
    unnest(predVars) %>%
    unnest(predVars) %>%
    select(Repeat,Fold,modelType,predOf,sireID,damID,Trait1,Trait2,predVar,Nsegsnps)
  if(modelType=="A"){ cvout %<>% mutate(predOf="VarBV") }

  if(modelType=="AD"){
    cvout<-cvout %>%
      filter(predOf=="VarA") %>%
      # Breeding value variance predictions from the predOf=="VarA"
      mutate(predOf="VarBV") %>%
      # for Total Gen Value variance predictions, need to compute:
      ## predVarTot = predVarA + predVarD
      bind_rows(cvout %>%
                  group_by(Repeat,Fold,modelType,sireID,damID,Trait1,Trait2) %>%
                  summarize(predVar=sum(predVar),
                            Nsegsnps=max(Nsegsnps),.groups = 'drop') %>%
                  mutate(predOf="VarTGV"))
  }
  cvout %<>%
    nest(predVars=c(sireID,damID,Trait1,Trait2,predVar,Nsegsnps))

  if(selInd==TRUE){
    # compute predicted selection index variances
    cvout %<>%
      ## loop over each rep-fold-predType
      mutate(predVars=map(predVars,function(predVars){
        predvars<-predVars %>%
          nest(fampredvars=c(-sireID,-damID,-Nsegsnps)) %>%
          ## internal loop over each family
          mutate(predVar=map_dbl(fampredvars,function(fampredvars){
            gmat<-fampredvars %>%
              pivot_wider(names_from = "Trait2",
                          values_from = "predVar") %>%
              column_to_rownames(var = "Trait1") %>%
              as.matrix
            gmat[lower.tri(gmat)]<-t(gmat)[lower.tri(gmat)]
            gmat %<>% .[names(SIwts),names(SIwts)]
            predSelIndVar<-SIwts%*%gmat%*%SIwts
            return(predSelIndVar) }),
            Trait1="SELIND",
            Trait2="SELIND") %>%
          select(-fampredvars)
        ## add sel index predictions to component trait
        ## var-covar predictions
        predvars %<>%
          bind_rows(.,predVars) %>%
          select(sireID,damID,Trait1,Trait2,predVar,Nsegsnps)
        return(predvars) }))
  }

  out %<>%
    mutate(predOf=ifelse(predOf=="GEBV","VarBV","VarTGV")) %>%
    left_join(cvout)
  out %<>%
    mutate(predVSobs=map2(predVars,obsVars,
                          ~left_join(.x,.y) %>%
                            nest(predVSobs=c(sireID,damID,predVar,obsVar,famSize,Nsegsnps)))) %>%
    select(-predVars,-obsVars) %>%
    unnest(predVSobs) %>%
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
  # Extract and format the GBLUPs from the marker effects object
  gblups<-snpeffs %>%
    unnest(effects) %>%
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

  # for modelType=="A" remove the GETGV as equiv. to GEBV
  if(modelType=="A"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~pivot_longer(.,cols = c(GEBV,GETGV),
                                              names_to = "predOf",
                                              values_to = "GBLUP") %>%
                                  nest(gblups=-predOf) %>%
                                  filter(predOf=="GEBV")))
  }
  # for modelType=="AD" remove the GEDD, pivot to long form GEBV/GETGV
  if(modelType=="AD"){
    out %<>%
      mutate(testset_gblups=map(testset_gblups,
                                ~select(.,-GEDD) %>%
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
    select(-CrossesToPredict,-gblups)

  out %<>%
    # outer loop over rep-fold-predtype
    mutate(obsMeans=map(famgblups,function(famgblups){
      return(famgblups %>%
               # inner loop over families
               mutate(obsmeans=map(gblupmat,
                                   function(gblupmat){
                                     gblupmeans<-colMeans(gblupmat) %>% as.list
                                     if(selInd==TRUE){
                                       selIndMean<-list(SELIND=as.numeric(gblupmeans[names(SIwts)])%*%SIwts)
                                       gblupmeans<-c(selIndMean,gblupmeans)
                                     }
                                     obsmeans<-tibble(Trait=names(gblupmeans),
                                                      obsMean=as.numeric(gblupmeans))
                                     return(obsmeans) }),
                      famSize=map_dbl(gblupmat,nrow)) %>%
               select(-gblupmat) %>%
               unnest(obsmeans))})) %>%
    select(-famgblups)

  cvout<-crossValOut %>%
    unnest(predMeans) %>%
    select(Repeat,Fold,modelType,predOf,sireID,damID,Trait,predMean) %>%
    nest(predMeans=c(sireID,damID,Trait,predMean))

  if(selInd==TRUE){
    # compute predicted selection index variances
    cvout %<>%
      ## loop over each rep-fold-predType
      mutate(predMeans=map(predMeans,function(predMeans){
        predmeans<-predMeans %>%
          pivot_wider(names_from = "Trait",
                      values_from = "predMean")
        predmeans %<>%
          select(sireID,damID) %>%
          mutate(Trait="SELIND",
                 predMean=(predmeans %>%
                             select(any_of(names(SIwts))) %>%
                             as.matrix(.)%*%SIwts)) %>%
          ## add sel index predictions to component trait
          ## mean predictions
          bind_rows(predMeans)
        return(predmeans) }))
  }

  out %<>%
    mutate(predOf=ifelse(predOf=="GEBV","MeanBV","MeanTGV")) %>%
    left_join(cvout)

  out %<>%
    mutate(predVSobs=map2(predMeans,obsMeans,
                          ~left_join(.x,.y) %>%
                            nest(predVSobs=c(sireID,damID,predMean,obsMean,famSize)))) %>%
    select(-predMeans,-obsMeans) %>%
    unnest(predVSobs) %>%
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
