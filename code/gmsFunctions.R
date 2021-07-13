# Clean TP data
readDBdata<-function(phenotypeFile,metadataFile=NULL){
  indata<-read.csv(phenotypeFile,
                   na.strings = c("#VALUE!",NA,".",""," ","-","\""),
                   stringsAsFactors = F)
  if(!is.null(metadataFile)){
    meta<-read.csv(metadataFile,
                   na.strings = c("#VALUE!",NA,".",""," ","-","\""),
                   stringsAsFactors = F) %>%
      rename(programName=breedingProgramName,
             programDescription=breedingProgramDescription,
             programDbId=breedingProgramDbId)
    indata<-left_join(indata,meta) }
  indata %<>%
    filter(observationLevel=="plot")
  return(indata) }

makeTrialTypeVar<-function(indata){
  # So far, this function is not very general
  # Handles IITA and NRCRI trial names as of September 2020.
  # Can customize this or add lines to grab TrialTypes for each breeding program
  if(indata$programName=="IITA"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("CE|clonal|13NEXTgenC1",studyName,ignore.case = T),"CET",NA),
             TrialType=ifelse(grepl("EC",studyName,ignore.case = T),"ExpCET",TrialType),
             TrialType=ifelse(grepl("PYT",studyName,ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("AYT",studyName,ignore.case = T),"AYT",TrialType),
             TrialType=ifelse(grepl("UYT",studyName,ignore.case = T),"UYT",TrialType),
             TrialType=ifelse(grepl("geneticgain|gg|genetic gain",studyName,ignore.case = T),"GeneticGain",TrialType),
             TrialType=ifelse(grepl("Cassava",studyName,ignore.case = T) & grepl("/",studyName),"GeneticGain",TrialType),
             TrialType=ifelse((grepl("clonal evaluation trial",!grepl("genetic gain",studyDescription,ignore.case = T),
                                     ignore.case = T)),"CET",TrialType),
             TrialType=ifelse(grepl("preliminary yield trial",studyDescription,ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("Crossingblock|\\.CB\\.|cross",studyName) & is.na(TrialType),"CrossingBlock",TrialType),
             TrialType=ifelse(grepl("NCRP",studyName) & is.na(TrialType),"NCRP",TrialType),
             TrialType=ifelse(grepl("conservation",studyName) & is.na(TrialType),"Conservation",TrialType),
             TrialType=ifelse(grepl("seedling|\\.SN",studyName),"SN",TrialType)) }
  if(indata$programName=="NRCRI"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("TP1",studyName,ignore.case = T),"TP1",NA),
             TrialType=ifelse(grepl("TP2",studyName,ignore.case = T),"TP2",TrialType),
             TrialType=ifelse(grepl("C1a",studyName,ignore.case = T),"C1a",TrialType),
             TrialType=ifelse(grepl("C1b",studyName,ignore.case = T),"C1b",TrialType),
             TrialType=ifelse(grepl("C2a",studyName,ignore.case = T),"C2a",TrialType),
             TrialType=ifelse(grepl("C2b",studyName,ignore.case = T),"C2b",TrialType),
             TrialType=ifelse(grepl("NCRP",studyName) & is.na(TrialType),"NCRP",TrialType),
             TrialType=ifelse(grepl("15nextgen60gs-cbUM|crossnblk|crossingblock",studyName,ignore.case = T) &
                                !grepl("CET",studyName),
                              "CrossingBlock",TrialType),
             TrialType=ifelse(grepl("seedling",studyName,ignore.case = T),NA,TrialType)) }

  if(indata$programName=="TARI"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("Advanced Yield|AYT", trialType, ignore.case = T),"AYT",NA),
             TrialType=ifelse(grepl("Clonal|CET", trialType, ignore.case = T),"CET",TrialType),
             TrialType=ifelse(grepl("Preliminary|PYT", trialType, ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("Regional", trialType, ignore.case = T),"RegionalTrial",TrialType),
             TrialType=ifelse(grepl("Uniform|UYT", trialType, ignore.case = T),"UYT",TrialType),
             TrialType=ifelse(grepl("Variety Release", trialType, ignore.case = T),"VarietyTrial",TrialType),
             TrialType=ifelse(grepl("CROSSING", trialType, ignore.case = T),"CrossingBlock",TrialType),
             TrialType=ifelse(grepl("GWAS", trialType, ignore.case = T),"GWAS",TrialType)) }

  return(outdata) }


#' @param  traitabbrevs data.frame with 2 cols (TraitAbbrev and TraitName). TraitName should match exactly to cassava ontology names
#' @param  indata data.frame read from cassavabase download
#' @param  customColsToKeep char. vec. of any custom cols you added and want to keep
renameAndSelectCols<-function(traitabbrevs,indata,
                              customColsToKeep=NULL){
  outdata<-indata %>%
    select(any_of(c("studyYear","programName","locationName","studyName","studyDesign",
                    "plotWidth","plotLength","fieldSize","plantingDate","harvestDate",
                    "germplasmName","observationUnitDbId",
                    "replicate","blockNumber","plotNumber","rowNumber","colNumber","entryType",
                    "trialType","plantsPerPlot","numberBlocks","numberReps")),
           any_of(customColsToKeep),
           any_of(traitabbrevs$TraitName)) %>% ungroup() %>%
    mutate(across(any_of(traitabbrevs$TraitName), as.numeric)) %>% ungroup() %>%
    pivot_longer(cols = any_of(traitabbrevs$TraitName),
                 names_to = "TraitName",
                 values_to = "Value") %>%
    left_join(.,traitabbrevs) %>%
    select(-TraitName) %>%
    pivot_wider(names_from = TraitAbbrev,
                values_from = "Value")
  return(outdata) }

# Curate by trial
nestByTrials<-function(indata){
  nested_indata<-indata %>%
    # Create some explicitly nested variables including loc and year to nest with the trial data
    mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
           trialInLocYr=paste0(yearInLoc,"_",studyName),
           repInTrial=paste0(trialInLocYr,"_",replicate),
           blockInRep=paste0(repInTrial,"_",blockNumber)) %>%
    nest(TrialData=-c(programName,locationName,studyYear,TrialType,studyName))
  return(nested_indata)
}

detectExptDesigns<-function(indata){
  # nestByTrials
  nestedDBdata<-indata %>%
    # Create some explicitly nested variables including loc and year to nest with the trial data
    mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
           trialInLocYr=paste0(yearInLoc,"_",studyName),
           repInTrial=paste0(trialInLocYr,"_",replicate),
           blockInRep=paste0(repInTrial,"_",blockNumber)) %>%
    nest(TrialData=-c(programName,locationName,studyYear,TrialType,studyName))

  # Define complete blocks
  nestedDBdata %>%
    mutate(Nobs=map_dbl(TrialData,~nrow(.)),
           MaxNOHAV=map_dbl(TrialData,~unique(.$MaxNOHAV)),
           Nrep=map_dbl(TrialData,~length(unique(.$replicate))),
           Nblock=map_dbl(TrialData,~length(unique(.$blockInRep))),
           Nclone=map_dbl(TrialData,~length(unique(.$germplasmName))),
           # median number of obs per clone
           medObsPerClone=map_dbl(TrialData,~count(.,germplasmName) %$% round(median(n),1)),
           # median number of obs per replicate
           medObsPerRep=map_dbl(TrialData,~count(.,replicate) %$% round(median(n),1)),
           # Define complete block effects based on the "replicate" variable
           CompleteBlocks=ifelse(Nrep>1 & medObsPerClone==Nrep & Nobs!=Nrep,TRUE,FALSE),
           # Additional trials with imperfect complete blocks
           CompleteBlocks=ifelse(Nrep>1 & medObsPerClone!=Nrep & medObsPerClone>1 & Nobs!=Nrep,TRUE,CompleteBlocks)) -> x
  x %>%
    # Some complete blocks may only be represented by the "blockNumber" column
    mutate(medBlocksPerClone=map_dbl(TrialData,~select(.,blockInRep,germplasmName) %>%
                                       # median number of blockInRep per clone
                                       distinct %>%
                                       count(germplasmName) %$%
                                       round(median(n))),
           # If CompleteBlocks==FALSE (complete blocks not detected based on replicate)
           # and if more than half the clones are represented in more than one block based on the blockInRep variable
           # Copy the blockInRep values into the repInTrial column
           # Recompute Nrep
           # and declare CompleteBlocks==TRUE
           TrialData=ifelse(medBlocksPerClone>1 & CompleteBlocks==FALSE,map(TrialData,~mutate(.,repInTrial=blockInRep)),TrialData),
           Nrep=map_dbl(TrialData,~length(unique(.$repInTrial))),
           CompleteBlocks=ifelse(medBlocksPerClone>1 & CompleteBlocks==FALSE,TRUE,CompleteBlocks)) -> y
  # Define incomplete blocks
  y %>%
    mutate(repsEqualBlocks=map_lgl(TrialData,~all(.$replicate==.$blockNumber)),
           NrepEqualNblock=ifelse(Nrep==Nblock,TRUE,FALSE),
           medObsPerBlockInRep=map_dbl(TrialData,~count(.,blockInRep) %$% round(median(n),1))) -> z
  # Define complete blocked trials with nested sub-blocks
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==TRUE & Nobs!=Nblock & Nblock>1 & medObsPerBlockInRep>1 & NrepEqualNblock==FALSE,TRUE,FALSE))
  # Define clearly unreplicated (CompleteBlocks==FALSE & Nrep==1) trials with nested sub-blocks
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & Nobs!=Nblock & Nblock>1 & medObsPerBlockInRep>1 & Nrep==1,TRUE,IncompleteBlocks))
  # Define additional trials with incomplete blocks (blockInRep) where CompleteBlocks==FALSE but Nrep>1 and Nrep==Block
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & IncompleteBlocks==FALSE &
                                     Nobs!=Nblock & Nblock>1 &  Nobs!=Nrep &
                                     medObsPerBlockInRep>1 & Nrep>1 & NrepEqualNblock==TRUE,TRUE,IncompleteBlocks))
  # Last few cases (2 trials actually) where Nrep>1 and Nblock>1 and Nrep!=Nblock but CompleteBlocks==FALSE
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & IncompleteBlocks==FALSE &
                                     Nobs!=Nblock & Nobs!=Nrep &
                                     medObsPerBlockInRep>1 & Nrep>1,TRUE,IncompleteBlocks))
  z %<>%
    dplyr::select(-MaxNOHAV) %>%
    unnest(TrialData)
  return(z)
}

nestDesignsDetectedByTraits<-function(indata,traits){
  indata %<>%
    select(programName,locationName,studyYear,TrialType,studyName,
           CompleteBlocks,IncompleteBlocks,
           yearInLoc,trialInLocYr,repInTrial,blockInRep,observationUnitDbId,
           germplasmName,FullSampleName,GID,all_of(traits),PropNOHAV) %>%
    mutate(IncompleteBlocks=ifelse(IncompleteBlocks==TRUE,"Yes","No"),
           CompleteBlocks=ifelse(CompleteBlocks==TRUE,"Yes","No")) %>%
    pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") %>%
    filter(!is.na(Value),
           !is.na(GID)) %>%
    nest(MultiTrialTraitData=c(-Trait))
  return(indata)
}

nestTrialsByTrait<-function(indata,traits){
  nested_trialdata<-dbdata %>%
    select(-MaxNOHAV) %>%
    unnest(TrialData) %>%
    pivot_longer(cols = any_of(traits),
                 names_to = "Trait",
                 values_to = "TraitValue") %>%
    nest(TraitByTrialData=-c(Trait,studyYear,programName,locationName,studyName,TrialType))
  return(nested_trialdata)
}

calcPropMissing<-function(TraitValues){ length(which(is.na(TraitValues))) / length(TraitValues) }

curateTrialOneTrait<-function(Trait,TraitByTrialData,GID="GID"){
  require(lme4)

  modelFormula<-paste0("TraitValue ~ (1|",GID,")")
  modelFormula<-ifelse(all(TraitByTrialData$CompleteBlocks),
                       paste0(modelFormula,"+(1|repInTrial)"),modelFormula)
  modelFormula<-ifelse(all(TraitByTrialData$IncompleteBlocks),
                       paste0(modelFormula,"+(1|blockInRep)"),modelFormula)
  modelFormula<-ifelse(grepl("logRTNO",Trait) | grepl("logFYLD",Trait) | grepl("logTOPYLD",Trait),
                       paste0(modelFormula,"+PropNOHAV"),modelFormula)

  propMiss<-calcPropMissing(TraitByTrialData$TraitValue)
  fit_model<-possibly(function(modelFormula,TraitByTrialData){
    model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData)
    if(!is.na(model_out)){
      outliers<-which(abs(rstudent(model_out))>=3.3)
      if(length(outliers)>0){
        model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData,
                        subset=abs(rstudent(model_out))<3.3)
      }
    }
    return(list(model_out=model_out,outliers=outliers)) },
    otherwise = NA)
  model_out<-fit_model(modelFormula,TraitByTrialData)
  if(is.na(model_out)){
    out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula,Noutliers=NA,Outliers=NA,propMiss=propMiss)
  } else {
    varcomps<-as.data.frame(VarCorr(model_out[["model_out"]]))[,c("grp","vcov")] %>%
      spread(grp,vcov)
    Vg<-varcomps$GID
    H2<-Vg/(Vg+varcomps$Residual)
    BLUP<-ranef(model_out[["model_out"]], condVar=TRUE)[[GID]]
    PEV <- c(attr(BLUP, "postVar"))
    blups<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
      mutate(REL=1-(PEV/Vg),
             drgBLUP=BLUP/REL,
             WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
    out <- tibble(H2=H2,
                  VarComps=list(varcomps),
                  BLUPs=list(blups),
                  Model=modelFormula,
                  Noutliers=length(model_out[["outliers"]]),
                  Outliers=list(model_out[["outliers"]]),
                  propMiss=propMiss) }
  return(out)
}

curateTrialsByTrait<-function(nestedTrialData,traits){
  outdata<-nestedTrialData %>%
    mutate(modelOutput=map2(Trait,TraitByTrialData,~curateTrialOneTrait(Trait = .x,TraitByTrialData = .y))) %>%
    dplyr::select(-TraitByTrialData) %>%
    unnest(modelOutput)
  return(outdata)
}

# Get BLUPs
nestForMultiTrialAnalysis<-function(curatedTrialData){
  nested_trialdata<-curatedTrialData %>%
    # remove trait-trial models that failed
    filter(!is.na(H2)) %>%
    # remove some per-trial summaries we don't want at this stage
    select(-H2,-VarComps,-Model,-Noutliers,-propMiss) %>%
    unnest(BLUPs) %>%
    nest(MultiTrialTraitData=c(-Trait))
  return(nested_trialdata)
}

fitMultiTrialModel<-function(curatedTrialData,GID="GID"){
  require(lme4)
  modelFormula<-paste0("drgBLUP ~ (1|",GID,")")
  fit_model<-possibly(function(modelFormula,curatedTrialData){
    model_out<-lmer(as.formula(modelFormula),
                    data=curatedTrialData,
                    weights = WT)
    return(model_out) },
    otherwise = NA)
  model_out<-fit_model(modelFormula,curatedTrialData)
  summary(model_out)
  if(is.na(model_out)){
    out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula)
  } else {
    varcomps<-as.data.frame(VarCorr(model_out))[,c("grp","vcov")] %>%
      spread(grp,vcov)
    Vg<-varcomps$GID
    H2<-Vg/(Vg+varcomps$Residual)
    BLUP<-ranef(model_out, condVar=TRUE)[[GID]]
    PEV <- c(attr(BLUP, "postVar"))
    blups<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
      mutate(REL=1-(PEV/Vg),
             drgBLUP=BLUP/REL,
             WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
    out <- tibble(H2=H2,
                  VarComps=list(varcomps),
                  BLUPs=list(blups),
                  Model=modelFormula) }
  return(out)
}

# Cross-validation

maf_filter<-function(snps,thresh){
  freq<-colMeans(snps, na.rm=T)/2; maf<-freq;
  maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
  snps1<-snps[,which(maf>thresh)];
  return(snps1) }

#' kinship function
#'
#' Function to create additive and dominance genomic relationship matrices from biallelic dosages.
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID.
#' @param type string, "add" or "domGenotypic" and "domClassic". type="add" gives same as rrBLUP::A.mat(), i.e. Van Raden, Method 1. type="domClassic" and type="domGenotypic" give the classical and genotypic parameterization according to Vitezica et al. 2013. Genetics. Both dominance matrices, in combo with the additive matrix, predict total merit identically. Difference in partition of add-dom, meaning of variance components and genomic predictions.
#' @return square symmetric genomic relationship matrix
#' @export
#'
#' @examples
#' K<-kinship(M,"add")
kinship<-function(M,type){
  M<-round(M)
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  if(type=="add"){
    Z <- M-2*P
    varD<-sum(2*freq*(1-freq))
    K <- tcrossprod(Z)/ varD
    return(K)
  }
  if(type=="domGenotypic"){
    W<-M; W[which(W==2)]<-0;
    W <- W-(2*P*(1-P))
    varD<-sum((2*freq*(1-freq))*(1-(2*freq*(1-freq))))
    D <- tcrossprod(W) / varD
    return(D)
  }
  if(type=="domClassic"){
    W<-M;
    W[which(W==1)]<-2*P[which(W==1)];
    W[which(W==2)]<-(4*P[which(W==2)]-2);
    W <- W-2*(P^2)
    varD<-sum((2*freq*(1-freq))^2)
    D <- tcrossprod(W) / varD
    return(D)
  }
}

#' @param byGroup logical, if TRUE, assumes a column named "Group" is present which unique classifies each GID into some genetic grouping.
#' @param modelType string, A, AD or ADE representing model with Additive-only, Add. plus Dominance, and Add. plus Dom. plus. AxD Epistasis (AD), respectively.
#' @param grms list of GRMs where each element is named either A, D, or, AD. Matrices supplied must match required by A, AD and ADE models. For ADE grms=list(A=A,D=D,AD=AD)...
#' @param augmentTP option to supply an additional set of training data, which will be added to each training model but never included in the test set.
#' @param TrainTestData data.frame with de-regressed BLUPs, BLUPs and weights (WT) for training and test. If byGroup==TRUE, a column with Group as the header uniquely classifying GIDs into genetic groups, is expected.
runCrossVal<-function(TrainTestData,modelType,grms,nrepeats,nfolds,ncores=1,
                          byGroup=FALSE,augmentTP=NULL,gid="GID",...){
  require(sommer); require(rsample)
  # Set-up replicated cross-validation folds
  # splitting by clone (if clone in training dataset, it can't be in testing)
  if(byGroup){
    cvsamples<-tibble(GroupName=unique(TrainTestData$Group))
  } else { cvsamples<-tibble(GroupName="None") }
  cvsamples<-cvsamples %>%
    mutate(Splits=map(GroupName,function(GroupName){
      if(GroupName!="None"){
        thisgroup<-TrainTestData %>%
          filter(Group==GroupName) } else { thisgroup<-TrainTestData }
      out<-tibble(repeats=1:nrepeats,
                  splits=rerun(nrepeats,group_vfold_cv(thisgroup, group = gid, v = nfolds))) %>%
        unnest(splits)
      return(out)
    })) %>%
    unnest(Splits)

  ## Internal function
  ## fits prediction model and calcs. accuracy for each train-test split

  fitModel<-possibly(function(splits,modelType,augmentTP,TrainTestData,GroupName,grms){
    starttime<-proc.time()[3]
    # Set-up training set
    trainingdata<-training(splits)
    ## Make sure, if there is an augmentTP, no GIDs in test-sets
    if(!is.null(augmentTP)){
      ## remove any test-set members from augment TP before adding to training data
      training_augment<-augmentTP %>% filter(!(!!sym(gid) %in% testing(splits)[[gid]]))
      trainingdata<-bind_rows(trainingdata,training_augment) }
    if(GroupName!="None"){ trainingdata<-bind_rows(trainingdata,
                                                   TrainTestData %>%
                                                     filter(Group!=GroupName,
                                                            !(!!sym(gid) %in% testing(splits)[[gid]]))) }
    # Subset kinship matrices
    traintestgids<-union(trainingdata[[gid]],testing(splits)[[gid]])
    A1<-grms[["A"]][traintestgids,traintestgids]
    trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[[gid]],levels=rownames(A1))
    if(modelType %in% c("AD","ADE")){
      D1<-grms[["D"]][traintestgids,traintestgids]
      trainingdata[[paste0(gid,"d")]]<-factor(trainingdata[[gid]],levels=rownames(D1))
      if(modelType=="ADE"){
        #AA1<-grms[["AA"]][traintestgids,traintestgids]
        AD1<-grms[["AD"]][traintestgids,traintestgids]
        diag(AD1)<-diag(AD1)+1e-06
        #DD1<-grms[["DD"]][traintestgids,traintestgids]
        #trainingdata[[paste0(gid,"aa")]]<-factor(trainingdata[[gid]],levels=rownames(AA1))
        trainingdata[[paste0(gid,"ad")]]<-factor(trainingdata[[gid]],levels=rownames(AD1))
        #trainingdata[[paste0(gid,"dd")]]<-factor(trainingdata[[gid]],levels=rownames(DD1))
      }
    }
    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A1)")
    if(modelType %in% c("AD","ADE")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D1)")
      if(modelType=="ADE"){
        randFormula<-paste0(randFormula,"+vs(",gid,"ad,Gu=AD1)")
        #"+vs(",gid,"aa,Gu=AA1)",
        #"+vs(",gid,"ad,Gu=AD1)")
        #"+vs(",gid,"dd,Gu=DD1)")
      }
    }
    # Fit genomic prediction model
    fit <- mmer(fixed = drgBLUP ~1,
                random = as.formula(randFormula),
                weights = WT,
                data=trainingdata)
    # Gather the BLUPs
    gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                   GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP))
    if(modelType %in% c("AD","ADE")){
      gblups %<>% mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP))
      if(modelType=="ADE"){
        gblups %<>% mutate(#GEEDaa=as.numeric(fit$U[[paste0("u:",gid,"aa")]]$drgBLUP),
          GEEDad=as.numeric(fit$U[[paste0("u:",gid,"ad")]]$drgBLUP))
        #GEEDdd=as.numeric(fit$U[[paste0("u:",gid,"dd")]]$drgBLUP))
      }
    }
    # Calc GETGVs
    ## Note that for modelType=="A", GEBV==GETGV
    gblups %<>%
      mutate(GETGV=rowSums(.[,grepl("GE",colnames(.))]))
    # Test set validation data
    validationData<-TrainTestData %>%
      dplyr::select(gid,BLUP) %>%
      filter(GID %in% testing(splits)[[gid]])
    # Measure accuracy in test set
    ## cor(GEBV,BLUP)
    ## cor(GETGV,BLUP)
    accuracy<-gblups %>%
      mutate(GETGV=rowSums(.[,grepl("GE",colnames(.))])) %>%
      filter(GID %in% testing(splits)[[gid]]) %>%
      left_join(validationData) %>%
      summarize(accGEBV=cor(GEBV,BLUP, use = 'complete.obs'),
                accGETGV=cor(GETGV,BLUP, use = 'complete.obs'))
    computeTime<-proc.time()[3]-starttime
    accuracy %<>% mutate(computeTime=computeTime)
    return(accuracy)
  },otherwise = NA)
  ## Run models across all train-test splits
  ## Parallelize
  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

  cvsamples<-cvsamples %>%
    mutate(accuracy=future_map2(splits,GroupName,
                                ~fitModel(splits=.x,GroupName=.y,
                                          modelType=modelType,augmentTP=NULL,TrainTestData=TrainTestData,grms=grms),
                                .progress = FALSE)) %>%
    unnest(accuracy)
  return(cvsamples)
}

#' @param blups nested data.frame with list-column "TrainingData" containing BLUPs
#' @param modelType string, values = "A", "AD", "DirDom". "A" = additive-only. "AD" = the classical additive - dominane partition leading to GETGV=GEBV+GEDD. NEW: modelType="DirDom", to fit a genome-wide homozyg. effect and then to subsequently do the work to incorporate a mean dominance effect into the predicted GEBV and GETGV; involves backsolve SNP effects and adding homz. effects.
#' @param grms list of GRMs. Any genotypes in the GRMs get predicted with, or without phenotypes. Each element is named either A or D. Matrices supplied must match required by A, AD and DirDom models. e.g. grms=list(A=A,D=D).
runGenomicPredictions<-function(modelType,
                                selInd,SIwts = NULL,
                                getMarkEffs=FALSE,
                                returnPEV=FALSE,
                                blups,grms,dosages=NULL,gid="GID",
                                ncores=1,
                                nBLASthreads=NULL){

  fitModel<-function(trainingdata,modelType,getMarkEffs,returnPEV,
                     gid="GID",grms,dosages,
                     nBLASthreads,...){
    if(!is.null(nBLASthreads)) { RhpcBLASctl::blas_set_num_threads(nBLASthreads) }
    # workers in plan(multisession) need this call internal to the function, it seems.

    A<-grms[["A"]]
    if(modelType %in% c("AD","DirDom")){ D<-grms[["D"]] }

    trainingdata %<>%
      dplyr::rename(gid=!!sym(gid)) %>%
      filter(gid %in% rownames(A))

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
                        date.warning = F,
                        getPEV = returnPEV)

    if(getMarkEffs==TRUE | modelType=="DirDom"){

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
      } }

    # Gather the GBLUPs
    if(modelType %in% c("A","AD")){
      gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                     GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP))
      if(returnPEV){
        pev_bv<-diag((fit$PevU[[paste0("u:",gid,"a")]]$drgBLUP))
        gblups %<>% left_join(tibble(GID=names(pev_bv),PEVbv=pev_bv))
      }
    }

    if(modelType=="AD"){
      gblups %<>% # compute GEDD (genomic-estimated dominance deviation)
        mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP),
               # compute GETGV
               GETGV=rowSums(.[,grepl("GE",colnames(.))]))
      if(returnPEV){
        pev_dd<-diag((fit$PevU[[paste0("u:",gid,"d")]]$drgBLUP))
        gblups %<>% left_join(tibble(GID=names(pev_dd),PEVdd=pev_dd))
      }
    }
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
      if(returnPEV){
        pev_a<-diag((fit$PevU[[paste0("u:",gid,"a")]]$drgBLUP))
        pev_dstar<-diag((fit$PevU[[paste0("u:",gid,"d_star")]]$drgBLUP))
        gblups %<>%
          left_join(tibble(GID=names(pev_a),PEVd_star=pev_a)) %>%
          left_join(tibble(GID=names(pev_dstar),PEVd_star=pev_dstar))
      }
    }

    # Extract variance components
    varcomps<-summary(fit)$varcomp

    # Exract fixed effects
    # for modelType="DirDom", contains estimate of genome-wide homoz. effect
    fixeffs<-summary(fit)$betas

    results<-tibble(gblups=list(gblups),
                    varcomps=list(varcomps),
                    fixeffs=list(fixeffs))

    if(getMarkEffs==TRUE){
      # Add snpeffects to output
      results %<>% mutate(allelesubsnpeff=list(allelesubsnpeff))
      if(modelType=="AD"){ results %<>% mutate(domdevsnpeff=list(domdevsnpeff)) }
      if(modelType=="DirDom"){
        results %<>% mutate(addsnpeff=list(addsnpeff),
                            domstar_snpeff=list(domstar_snpeff),
                            domsnpeff=list(domsnpeff)) } }
    # this is to remove conflicts with dplyr function select() downstream
    detach("package:sommer",unload = T); detach("package:MASS",unload = T)
    # return results
    return(results)
  }

  require(furrr); plan(multisession, workers = ncores)
  options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")
  predictions<-blups %>%
    mutate(genomicPredOut=future_map(TrainingData,
                                     ~fitModel(trainingdata=.,
                                               modelType=modelType,
                                               getMarkEffs=getMarkEffs,
                                               returnPEV=returnPEV,
                                               gid=gid,
                                               grms=grms,
                                               dosages=dosages,
                                               nBLASthreads=nBLASthreads)))

  predictions %<>%
    select(-TrainingData) %>%
    unnest(genomicPredOut)
  # tidy GBLUP output for e.g. breeders / selections
  gblups<-predictions %>%
    select(Trait,gblups) %>%
    unnest(gblups) %>%
    select(!!sym(gid),Trait,any_of(c("GEBV","GETGV"))) %>%
    pivot_longer(any_of(c("GEBV","GETGV")),
                 values_to = "GBLUP",
                 names_to = "predOf") %>%
    pivot_wider(names_from = "Trait",
                values_from = "GBLUP")
  if(selInd){
    # calc. SELIND and add to tidy output
    gblups %<>%
      mutate(SELIND=as.numeric((gblups %>%
                                  select(names(SIwts)) %>%
                                  as.matrix(.))%*%SIwts)) %>%
      relocate(SELIND, .after = predOf)
  }
  predictions<-tibble(gblups=list(gblups),
                      genomicPredOut=list(predictions))
  return(predictions)
}

predictCrosses<-function(modelType,
                         stdSelInt = 2.67,
                         selInd,SIwts = NULL,
                         CrossesToPredict,
                         snpeffs,dosages,
                         haploMat,recombFreqMat,ncores){
  ## Format SNP effect matrices ~~~~~~~~~~~~~~~~~~~~~~~~

  AlleleSubEffectList<-snpeffs$allelesubsnpeff %>%
    `names<-`(.,snpeffs$Trait) %>%
    map(.,~t(.))

  if(modelType=="AD"){
    # DomDevEffects for model "AD" to predict VarTGV = VarBV + VarDD
    DomDevEffectList=snpeffs$domdevsnpeff %>%
      `names<-`(.,snpeffs$Trait) %>%
      map(.,~t(.)) }

  if(modelType=="DirDom"){
    # AddEffectList + DomEffectList --> VarTGV; AlleleSubEffectList --> VarBV;
    AddEffectList<-snpeffs$addsnpeff %>%
      `names<-`(.,snpeffs$Trait) %>%
      map(.,~t(.))
    DomEffectList<-snpeffs$domsnpeff %>%
      `names<-`(.,snpeffs$Trait) %>%
      map(.,~t(.)) }

  ## Predict cross variances ~~~~~~~~~~~~~~~~~~~~~~~~
  print("Predicting cross variance parameters")
  if(modelType=="A"){
    predvars<-predCrossVars(CrossesToPredict=CrossesToPredict,
                            AddEffectList=AlleleSubEffectList,
                            modelType="A",
                            haploMat=haploMat,
                            recombFreqMat=recombFreqMat,
                            ncores=ncores) %>%
      unnest(predVars) %>%
      mutate(predOf="VarBV") }
  if(modelType=="AD"){
    predvars<-predCrossVars(CrossesToPredict=CrossesToPredict,
                            AddEffectList=AlleleSubEffectList,
                            DomEffectList=DomDevEffectList,
                            modelType="AD",
                            haploMat=haploMat,
                            recombFreqMat=recombFreqMat,
                            ncores=ncores) %>%
      unnest(predVars) %>%
      mutate(predOf=ifelse(predOf=="VarA","VarBV","VarDD"))
  }
  if(modelType=="DirDom"){
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
    predvars<-predVarBV %>%
      unnest(predVars) %>%
      mutate(predOf="VarBV") %>%
      bind_rows(predVarTGV %>%
                  unnest(predVars)) }

  ## Predict cross means ~~~~~~~~~~~~~~~~~~~~~~~~
  print("Predicting cross means")
  ### predict MeanBVs
  predmeans<-predCrossMeans(AddEffectList=AlleleSubEffectList,
                            CrossesToPredict=CrossesToPredict,
                            doseMat=dosages,ncores=ncores,predType="BV")
  if(modelType=="AD"){
    ### DO NOT predict MeanTGV ~but~ duplicate MeanBV as MeanTGV prediction
    ### there IS predVarTGV for this model, output predUC-TGV (i.e. UC_variety)
    predmeans %<>%
      bind_rows(predmeans %>% mutate(predOf="TGV")) }

  predmeans
  if(modelType=="DirDom"){
    ### predict MeanTGVs
    ####  Prediction of MeanTGV is only available for the DirDom model
    #### or a model with "genotypic" additive-dominance SNP effects
    #### As implemented, modelType="AD" is the "classical" partition (BVs+ DomDevs)
    predmeans %<>%
      bind_rows(predCrossMeans(AddEffectList=AddEffectList,
                               DomEffectList=DomEffectList,
                               CrossesToPredict=CrossesToPredict,
                               doseMat=dosages,ncores=ncores,
                               predType="TGV"))
  }

  ## SIMPLIFIED, TIDY, CROSS-WISE OUTPUT ~~~~~~~~~~~~~~~~~~~~~~~~
  # store raw mean and var preds
  rawPreds<-list(predMeans=list(predmeans),
                 predVArs=list(predvars))

  ## tidy pred. means ~~~~~~
  predmeans %<>%
    mutate(predOf=gsub("Mean","",predOf),
           Trait2=Trait) %>% # to match with variance pred. output
    rename(Trait1=Trait) %>% # to match with variance pred. output
    select(sireID,damID,predOf,Trait1,Trait2,predMean)

  ## tidy pred. vars ~~~~~~
  predvars %<>%
    select(sireID,damID,Nsegsnps,predOf,Trait1,Trait2,predVar) %>%
    mutate(predOf=gsub("Var","",predOf))
  if(modelType=="AD"){
    predvars %<>%
      filter(predOf=="BV") %>%
      bind_rows(predvars %>%
                  pivot_wider(names_from = "predOf",
                              values_from = "predVar",
                              names_prefix = "predVar") %>%
                  mutate(predVar=predVarBV+predVarDD,
                         predOf="TGV") %>%
                  select(-predVarBV,-predVarDD))
  }
  if(modelType=="DirDom"){
    predvars %<>%
      filter(predOf=="BV") %>%
      bind_rows(predvars %>%
                  filter(predOf!="BV") %>%
                  pivot_wider(names_from = "predOf",
                              values_from = "predVar",
                              names_prefix = "predVar") %>%
                  mutate(predVar=predVarA+predVarD,
                         predOf="TGV") %>%
                  select(-predVarA,-predVarD))
  }

  ## SELECTION INDEX MEANS AND VARIANCES ~~~~~~~~~~~~~~~~~~~~~~~~
  #### Compute and add to tidy output, if requested
  if(selInd){
    print("Computing SELECTION INDEX means and variances.")
    traits<-unique(predmeans$Trait1)
    ## Compute Mean SELIND
    predmeans %<>%
      select(-Trait2) %>%
      spread(Trait1,predMean) %>%
      select(sireID,damID,predOf,all_of(names(SIwts))) %>%
      mutate(SELIND=as.numeric((predmeans %>%
                                  select(-Trait2) %>%
                                  spread(Trait1,predMean) %>%
                                  select(all_of(names(SIwts))) %>%
                                  as.matrix(.))%*%SIwts)) %>%
      relocate(SELIND, .after = predOf) %>%
      pivot_longer(cols = c(SELIND,all_of(traits)),
                   names_to = "Trait1",
                   values_to = "predMean") %>%
      mutate(Trait2=Trait1) %>%
      select(sireID,damID,predOf,Trait1,Trait2,predMean)
    ## Compute Var SELIND
    require(furrr); plan(multisession, workers = ncores)
    options(future.globals.maxSize=+Inf); options(future.rng.onMisuse="ignore")

    predvars %<>%
      nest(predVars=c(Trait1,Trait2,predVar)) %>%
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
        return(predVars) })) %>%
      unnest(predVars)
  }

  ## USEFULNESS CRITERIA ~~~~~~~~~~~~~~~~~~~~~~~~
  tidyPreds<-predvars %>%
    inner_join(predmeans) %>%
    rename(Trait=Trait1) %>%
    select(sireID,damID,Nsegsnps,predOf,Trait,predMean,predVar) %>%
    mutate(predSD=sqrt(predVar),
           predUsefulness=predMean+(stdSelInt*predSD))

  predictions<-tibble(tidyPreds=list(tidyPreds),
                      rawPreds=list(rawPreds))
  return(predictions)
}
