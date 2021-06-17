runParentWiseCrossVal<-function(nrepeats,nfolds,seed=NULL,modelType,
                                ncores=1,outName=NULL,
                                ped=ped,gid="GID",blups,
                                dosages,grms,haploMat,recombFreqMat,
                                selInd,SIwts = NULL,...){
  initime<-proc.time()[3]

  ## Make parent-wise folds
  parentfolds<-makeParentFolds(ped=ped,gid="GID",
                               nrepeats=nrepeats,
                               nfolds=nfolds,
                               seed=seed)
  if(!is.null(outName)){ saveRDS(parentfolds,file=paste0(outName,"_parentfolds.rds")) }
  print("Set-up parent-wise cross-validation folds")

  ## Get univariate REML marker effects
  #### modelType=="AD"
  print("Fitting models to get marker effects")
  starttime<-proc.time()[3]
  markEffs<-getMarkEffs(parentfolds,blups=blups,gid=gid,modelType=modelType,
                        grms=grms,dosages=dosages,ncores=ncores)
  if(!is.null(outName)){ saveRDS(markEffs,file=paste0(outName,"_markerEffects.rds")) }
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
  cvPredVars %<>% select(-AddEffectList,-DomEffectList,-CrossesToPredict)
  if(!is.null(outName)){ saveRDS(cvPredVars,file=paste0(outName,"_predVars.rds")) }
  print(paste0("Cross variance parameters predicted. Took  ",
               round((proc.time()[3] - starttime)/60/60,5)," hrs"))

  print("Predicting cross means")
  starttime<-proc.time()[3]
  cvPredMeans<-predictCrossMeans(modelType=modelType,snpeffs=markEffs,
                                 ncores=ncores,
                                 parentfolds=parentfolds,doseMat=dosages)
  if(!is.null(outName)){ saveRDS(cvPredMeans,file=paste0(outName,"_predMeans.rds")) }
  print(paste0("Cross means predicted. Took  ",
               round((proc.time()[3] - starttime)/60/60,5)," hrs"))

  print("Compute prediction accuracies and wrap up.")
  ## Variance prediction accuracies
  starttime<-proc.time()[3]
  varPredAccuracy<-varPredAccuracy(modelType = modelType,
                                   crossValOut = cvPredVars,
                                   markEffs = markEffs,
                                   ped = ped,selInd = selInd,SIwts = SIwts)
  if(!is.null(outName)){ saveRDS(varPredAccuracy,file=paste0(outName,"_varPredAccuracy.rds")) }

  ## Mean prediction accuracies
  meanPredAccuracy<-meanPredAccuracy(modelType = modelType,
                                     crossValOut = cvPredMeans,
                                     markEffs = markEffs,
                                     ped = ped,selInd = selInd,SIwts = SIwts)
  if(!is.null(outName)){ saveRDS(meanPredAccuracy,file=paste0(outName,"_meanPredAccuracy.rds")) }

  accuracy_out<-list(meanPredAccuracy=meanPredAccuracy,
                     varPredAccuracy=varPredAccuracy)
  print(paste0("Accuracies predicted. Took  ",
               round((proc.time()[3] - initime)/60/60,5),
               " hrs total.\n Goodbye!"))
  return(accuracy_out)
}

meanPredAccuracy<-function(crossValOut,markEffs,ped,modelType,
                           selInd=FALSE,SIwts=NULL){

  # Extract and format the GBLUPs from the marker effects object
  gblups<-markEffs %>%
    filter(Dataset=="testset") %>%
    mutate(testset_gblups=map(effects,
                              function(effects){
                                gblups<-effects %>%
                                  mutate(gblups=map(modelOut,
                                                    ~select(.,gblups) %>%
                                                      unnest())) %>%
                                  select(-trainingdata,-modelOut)
                                return(gblups)})) %>%
    select(Repeat,Fold,modelType,testset_gblups)

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
                               ~semi_join(.x %>% unnest(gblups),.y)))
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


varPredAccuracy<-function(crossValOut,markEffs,ped,modelType,
                          selInd=FALSE,SIwts=NULL){

  # Extract and format the GBLUPs from the marker effects object
  gblups<-markEffs %>%
    filter(Dataset=="testset") %>%
    mutate(testset_gblups=map(effects,
                              function(effects){
                                gblups<-effects %>%
                                  mutate(gblups=map(modelOut,
                                                    ~select(.,gblups) %>%
                                                      unnest())) %>%
                                  select(-trainingdata,-modelOut)
                                return(gblups)})) %>%
    select(Repeat,Fold,modelType,testset_gblups)
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
                               ~semi_join(.x %>% unnest(gblups),.y)))
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
    mutate(obsVars=map(famgblups,function(famgblups){
      return(famgblups %>%
               # inner loop over families
               mutate(obsvars=map(gblupmat,
                                  function(gblupmat){
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
                            Nsegsnps=max(Nsegsnps)) %>%
                  ungroup() %>%
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

predictCrossMeans<-function(modelType,snpeffs,parentfolds,
                            doseMat,ncores){
  predmeans<-snpeffs %>%
    unnest(effects) %>%
    unnest(modelOut) %>%
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
    select(-contains("EffectsList"),-CrossesToPredict)
  return(predmeans) }

predictCrossVars<-function(modelType,snpeffs,parentfolds,
                           haploMat,recombFreqMat,ncores){
  predvars<-snpeffs %>%
    unnest(effects) %>%
    unnest(modelOut) %>%
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
                           ~predCrossVar(CrossesToPredict=.x,
                                         AddEffectList=.y,
                                         modelType=modelType,
                                         haploMat=haploMat,
                                         recombFreqMat=recombFreqMat,
                                         ncores=ncores))) }
  if(modelType=="AD"){
    predvars<-predvars %>%
      mutate(predVars=pmap(.,function(CrossesToPredict,
                                      AddEffectList,DomEffectList,...){
        out<-predCrossVar(CrossesToPredict=CrossesToPredict,
                            AddEffectList=AddEffectList,
                            DomEffectList=DomEffectList,
                            modelType=modelType,
                            haploMat=haploMat,
                            recombFreqMat=recombFreqMat,
                            ncores=ncores)
        return(out) })) }
  return(predvars)
}

getMarkEffs<-function(parentfolds,blups,gid,modelType,grms,dosages,ncores){
  traintestdata<-parentfolds %>%
    dplyr::select(Repeat,Fold,trainset,testset) %>%
    pivot_longer(c(trainset,testset),
                 names_to = "Dataset",
                 values_to = "sampleIDs")

  # Internal function
  ## For each training or testing chunk of sampleIDs
  ## fit GBLUP model for each trait
  ## Backsolve SNP-effects from GBLUP

  fitModel<-function(blups,sampleIDs,modelType,grms){
    require(predCrossVar)
    A<-grms[["A"]]
    if(modelType %in% c("AD")){ D<-grms[["D"]] }
    # debug fitModel()
    # sampleIDs<-traintestdata$sampleIDs[[1]]
    out<-blups %>%
      #slice(1:2) %>% # debug (just 2 traits)
      dplyr::mutate(trainingdata=map(blups,function(blups){
        trainingdata<-blups %>%
          rename(gid=!!sym(gid)) %>%
          filter(gid %in% sampleIDs)
        trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[["gid"]],
                                                levels=rownames(A))
        if(modelType %in% c("AD")){
          trainingdata[[paste0(gid,"d")]]<-trainingdata[[paste0(gid,"a")]]
          # factor(trainingdata[["gid"]],
          #                                         levels=rownames(D))
        }
        return(trainingdata) })) %>%
      dplyr::select(-blups)

    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A)")
    if(modelType %in% c("AD")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D)")
    }

    # Fit model across each trait
    out<-out %>%
      mutate(modelOut=map(trainingdata,function(trainingdata){
        # Fit genomic prediction model
        require(sommer)
        fit <- mmer(fixed = drgBLUP ~1,
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
      }))
    return(out)
  }

  require(furrr); options(mc.cores=ncores); plan(multicore)
  options(future.globals.maxSize=10000*1024^2)

  traintestdata<-traintestdata %>%
    #slice(1:2) %>% # debug just 2 chunks
    mutate(effects=future_map(sampleIDs,~fitModel(sampleIDs=.,
                                                  blups=blups,grms=grms,
                                                  modelType=modelType)),
           modelType=modelType)
  # this is to remove conflicts with dplyr function select() downstream
  detach("package:sommer",unload = T); detach("package:MASS",unload = T)

  return(traintestdata)
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
