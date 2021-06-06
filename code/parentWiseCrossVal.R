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
  return(traintestdata)
}

# Prunes out offspring, grandkids, greatgrandkids (up to X4) steps of
# great ancestors.  It is not automatically recursive across any depth of
# pedigree. That depth works for current test pedigree (IITA 2021).
# Must name parent columns in ped "sireID" and "damID".
makeParentFolds<-function(ped,gid,nrepeats=5,nfolds=5,seed=NULL){
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
    unnest(folds) %>%
    rename(Repeat=id,Fold=id2) %>%
    select(-splits)

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
