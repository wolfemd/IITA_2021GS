
#' Estimate selection error based an index
#'
#' for a single trial, compute a selection index using
#' available traits BLUPs
#' Compute the correlation, reg r-squared and MSE
#' \eqn{SELIND GETGV  ~ SELIND BLUPs_trial}
#'
#' @param TrialData data.frame trial data
#' @param CompleteBlocks T/F
#' @param IncompleteBlocks T/F
#' @param SIwts named vector index weights
#' @param getgvs data.frame
#' @param ...
#'
#' @return
estimateSelectionError<-possibly(function(TrialData,CompleteBlocks,IncompleteBlocks,
                                          SIwts,getgvs,...){

  # SET-UP THE DATA TRAIT-BY-TRIAL~~~~~~~~~~~~~~~~~~~~~
  blups<-TrialData %>%
    select(observationUnitDbId,GID,repInTrial,blockInRep,PropNOHAV,
           all_of(names(SIwts))) %>%
    pivot_longer(cols = c(all_of(names(SIwts))),
                 names_to = "Trait",
                 values_to = "TraitValue") %>%
    nest(TraitByTrialData=c(-Trait))

  # FIT MIXED-MODELS TRAIT-BY-TRIAL~~~~~~~~~~~~~~~~~~~~~~

  ## if model fails, by design, returns NULL for a given trait or trait-trial
  ## Output will simply be absent.
  fit_model<-possibly(function(Trait,TraitByTrialData,CompleteBlocks,IncompleteBlocks){
    # debug settings ~~~~~~~~~~
    # TraitByTrialData<-blups$TraitByTrialData[[1]]
    # Trait<-blups$Trait[[1]]
    # TraitByTrialData<-blups$TraitByTrialData[[8]]
    # Trait<-blups$Trait[[8]]
    # rm(TraitByTrialData,Trait)

    # Model formula based on trial design
    modelFormula<-paste0("TraitValue ~ (1|GID)")
    modelFormula<-ifelse(CompleteBlocks=="Yes",
                         paste0(modelFormula,"+(1|repInTrial)"),modelFormula)
    modelFormula<-ifelse(IncompleteBlocks=="Yes",
                         paste0(modelFormula,"+(1|blockInRep)"),modelFormula)
    modelFormula<-ifelse(Trait %in% c("logRTNO","logFYLD","logTOPYLD","logDYLD"),
                         paste0(modelFormula,"+PropNOHAV"),modelFormula)
    require(lme4)
    model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData)
    propMiss<-length(which(is.na(TraitByTrialData$TraitValue))) / length(TraitByTrialData$TraitValue)
    if(!exists("model_out")){
      out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula,propMiss=propMiss)
    } else {
      varcomps<-as.data.frame(VarCorr(model_out))[,c("grp","vcov")] %>%
        spread(grp,vcov)
      Vg<-varcomps$GID
      H2<-Vg/(Vg+varcomps$Residual)
      BLUP<-ranef(model_out, condVar=TRUE)[["GID"]]
      PEV <- c(attr(BLUP, "postVar"))
      blup<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
        mutate(REL=1-(PEV/Vg),
               drgBLUP=BLUP/REL,
               WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
      out <- tibble(H2=H2,
                    VarComps=list(varcomps),
                    BLUPs=list(blup),
                    Model=modelFormula,
                    propMiss=propMiss) }
    return(out) },
    otherwise = NULL)

  blups %<>%
    mutate(modelOut=pmap(.,fit_model,
                         CompleteBlocks=CompleteBlocks,
                         IncompleteBlocks=IncompleteBlocks)) %>%
    select(-TraitByTrialData) %>%
    unnest(modelOut) %>%
    unnest(VarComps)

  # COMPUTE SELIND FROM BLUPs~~~~~~~~~~~~~~~~~~~~~~
  si_blups<-blups %>%
    select(Trait,BLUPs) %>%
    unnest(BLUPs) %>%
    select(GID,Trait,BLUP) %>%
    spread(Trait,BLUP) %>%
    select(GID,any_of(names(SIwts))) %>%
    column_to_rownames(var = "GID") %>%
    as.matrix

  wts<-SIwts[colnames(si_blups)]
  adjWts<-wts*(sum(SIwts)/sum(wts))

  si_blups<-si_blups%*%adjWts %>%
    as.data.frame %>%
    rownames_to_column(var = "GID") %>%
    rename(SI_BLUP=V1) %>%
    left_join(getgvs)

  # Correlation between SELIND (SI GETGV) and the SI of BLUPs for current trial
  cor2si<-si_blups %$% cor(SI_BLUP,SELIND,use = 'complete.obs')

  # Regress SELIND on SI_BLUPs
  regSIonTrialBLUP<-lm(SELIND~SI_BLUP,data = si_blups)

  # TWO MEASURES OF SELECTION ERROR
  ## 1) regression r-squared [or 1 - r2_si actually]
  ## 2) mean squared error,
  ##### not sure if the scaling will make sense across trials
  ##### because of differences in traits included
  r2_si<-regSIonTrialBLUP %>% summary %$% r.squared

  mse<-regSIonTrialBLUP %>% anova %>% as.data.frame %>% .["Residuals","Mean Sq"]

  NcloneForReg<-si_blups %>% filter(!is.na(SI_BLUP),
                                    !is.na(SELIND)) %>% nrow()
  # Collect and return outputs for current trial
  trial_out<-tibble(cor2si=cor2si,
                    r2_si=r2_si,
                    TrialMSE=mse,
                    NcloneForReg=NcloneForReg,
                    SI_BLUPs=list(si_blups),
                    BLUPs=list(blups))
  return(trial_out)
},otherwise = NA)
