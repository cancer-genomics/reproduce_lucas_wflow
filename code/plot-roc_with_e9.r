plotrocs <- function(data, mytheme, textsize=3.5, facet=TRUE) {
    if(facet) {
        lab <- data %>%
            as_tibble() %>%
            mutate(lower=format(round(lower, 2), nsmall=2),
                   auc=format(round(auc, 2), nsmall=2),
                   upper=format(round(upper, 2), nsmall=2),
                   features=features) %>%
            group_by(features) %>%
            summarize(auc=unique(auc),
                      lower=unique(lower),
                      upper=unique(upper)) %>%
            mutate(text=paste0("AUC: ", auc, "(", lower, "-", upper, ")")) %>%
            mutate(x=0.7, y=0.05)
        ## Someone fix
        ##    lab <- data[,.(lower=format(round(unique(lower),2), nsmall=2),
        ##                   auc=format(round(unique(auc),2), nsmall=2),
        ##                   upper=format(round(unique(upper), 2), nsmall=2)),
        ##                by=features]
        ##n <- data[, length(unique(features))]
        ##        lab <- lab[order(features),
        ##                   `:=`(text=paste0("AUC: ", auc, " (", lower, " - ", upper, ")"),
        ##                        x=0.70, y=0.05)]#seq(0.4, 0.15, length.out=n))]
    }
    else {
        lab <- data[,.(lower=format(round(unique(lower),2), nsmall=2),
                       auc=format(round(unique(auc),2), nsmall=2),
                       upper=format(round(unique(upper), 2), nsmall=2))]
        lab <- lab[,
                   `:=`(text=paste0("AUC: ", auc, " (",
                                    lower, " - ", upper, ")"),
                        x=0.70, y=0.05)]#seq(0.4, 0.15, length.out=n))]
    }
    A <- ggplot() +
        geom_vline(xintercept=c(0.80, 0.90, 0.95),
                   color="gray80", size=0.5, linetype="dashed") +
        geom_abline(intercept=1, slope=1, color="gray", size=1) +
        geom_line(data=data, aes(spec, sens), size=1.1,
                  color="steelblue")
    if(facet) A <- A + facet_grid(.~features)
    A <- A + scale_x_reverse(expand=c(0, 0.01),
                             breaks=c(0, 0.25, 0.5, 0.80, 0.95, 1),
                             labels=as.character(
                                 c("0", ".25", ".50", ".80", ".95", ""))) +
        scale_y_continuous(expand=c(0, 0.01),
                           labels=as.character(
                               c("0", ".25", ".50", ".75", "1"))) +
        ## scale_color_manual(values=palette) +
        mytheme +
        geom_text(data=lab, aes(x, y,  label=text),
                  size=textsize, hjust=0) +
        xlab("Specificity") +
        ylab("Sensitivity") +
        ##ggtitle(overall.title) +
        guides(color=guide_legend(title="Approach"))
    A
}

rocstats <- function(obs, score) {
    roc <- pROC::roc
    roc <- suppressMessages(roc(obs, score,
                                levels=c("healthy", "cancer"),
                                ci=TRUE))
    list(sens = rev(roc[["sensitivities"]]),
         spec = rev(roc[["specificities"]]),
         auc = roc$auc,
         lower = roc$ci[1],
         upper = roc$ci[3])
}


sensitivityStats <- function(x, cutoff){
    tp <- sum(x$score.seq >= cutoff)
    fn <- sum(x$score.seq < cutoff)
    ci <- binom.test(tp, tp+fn, conf.level=0.9)$conf.int
    ##ci <- qbeta(c(0.05, 0.9), 0.5+tp, 0.5+fn)
    tibble(tp=tp, fn=fn, sens=tp/(tp+fn),
           n=tp+fn,
           `0.05%`=ci[1],
           `0.95%`=ci[2],
           stage=x$stage[1])
}

sensitivityByStage <- function(cutoff, pred.list){
    pred.list %>%
        map_dfr(sensitivityStats, cutoff=cutoff)
}

specificityStats <- function(cutoff, x){
    tn <- sum(x$score.seq <= cutoff)
    fp <- sum(x$score.seq >= cutoff)
    ci <- binom.test(tn, tn+fp, conf.level=0.9)$conf.int
    ci2 <- binom.test(tn, tn+fp, conf.level=0.95)$conf.int
    tibble(tn=tn,
           fp=fp,
           n=tn+fp,
           spec=tn/(tn+fp),
           ##`0.05%`=qbeta(0.05, 0.5+tn, 0.5+fp),
           `0.05%`=ci[1],
           ##`0.95%`=qbeta(0.95, 0.5+tn, 0.5+fp),
           ##`0.025%`=qbeta(0.025, 0.5+tn, 0.5+fp),
           ##`0.975%`=qbeta(0.975, 0.5+tn, 0.5+fp))
           `0.95%`=ci[2],
           `0.025%`=ci2[1],
           `0.0975%`=ci2[2])
}

performance2 <- function(params, THR){
    prevalence <- params[["prevalence"]]
    sens <- params[["sens"]]
    sens[ sens > 0.999 ] <- 0.999
    spec <- params[["spec"]]
    spec[ spec > 0.999 ] <- 0.999
    compliance <- params[["compliance"]]
    L <- length(compliance)
    screened <- rbinom(L, size=100e3, prob=compliance)
    ## prevalence
    P <- rbinom(L, size=screened, prob=prevalence) ## P = FN + TP
    ##
    N <- screened - P  ## N = TN + FP
    ##
    ## Pr(test + | y=1)
    TP <- rbinom(L, size=P, prob=sens)
    ## Pr(test + | y=0)
    FP <- rbinom(L, size=N, prob=(1-spec))
    ## Pr(test - | y=0)
    TN <- N - FP  ## N = FP + TN
    ## Pr(test - | y=1)
    ##FN <- rbinom(L, size=P, prob=1-sens)
    FN <- P - TP
    ##falseNegatives(study)
    ##acc <- accuracy(study)
    ##err <- errorRate(study)
    fpr <- FP/N
    ##fpr <- falsePositiveRate(study)
    fnr <- FN/P
    ##fnr <- falseNegativeRate(study)
    tnr <- TN/N
    ##tnr <- trueNegativeRate(study)
    ##sens <- tpr <- truePositiveRate(study)
    ##Pr(x = y)
    tpr <- TP/P
    acc <- (TP+TN)/(P+N)
    err <- (FP+FN)/(P+N)
    ppv <- TP/(TP+FP)
    npv <- TN/(TN+FN)
    ##ppv <- positivePredictive(TP, FP)
    ##npv <- negativePredictive(TN, FN)
    stats <- tibble("P"=P, "N"=N, "TP"=TP, "FP"=FP, "TN"=TN,
          "FN"=FN, "acc"=acc, "err"=err,
          "fpr"=fpr, "fnr"=fnr, "tnr"=tnr,
          "tpr"=tpr, "ppv"=ppv, "npv"=npv,
          "number_screened"=screened)
    stats
}

sensValidation <- function(spec_cutoff, x, xlabels) {
    x %>%
        filter(specificity==spec_cutoff, statistic=="sensitivity") %>%
        ggplot(aes(stage, performance)) +
        geom_hline(yintercept=1, linetype="dashed") +
        geom_errorbar(aes(ymin=`0.05%`, ymax=`0.95%`,
                          color=set), width=0.1,
                      position=position_dodge(0.2)) +
        geom_point(size=3, shape=21,
                   aes(color=set,
                       fill=set),
                   position=position_dodge(0.2)) +
        theme_classic(base_size=20) +
        theme(legend.position=c(0.7, 0.25),
              legend.text=element_text(size=15),
              panel.background=element_rect(fill="white",
                                            color="gray30"),
              panel.grid=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        scale_color_manual(values=c("gray30", "steelblue")) +
        scale_fill_manual(values=c("gray30", "steelblue")) +
        ylab("") +
        xlab("") +
        scale_y_continuous(labels=scales::percent_format(accuracy=1),
                           limit=c(0, 1)) +
        scale_x_discrete(labels=xlabels) +
        ##guides(fill=guide_legend(title=""), color=guide_legend(title="")) +
        guides(fill=FALSE, color=FALSE) +
        xlab("")

}

specValidation <- function(spec_cutoff, x, xlabels){
    x %>%
        filter(specificity==spec_cutoff, statistic=="specificity") %>%
        pivot_longer(cols=c(performance, specificity),
                     names_to="perf", values_to="stat") %>%
        mutate(set=c("Valdation set", "Training set (LUCAS)")) %>%
        mutate(`0.05%`=ifelse(set=="Training set (LUCAS)", stat, `0.05%`),
               `0.95%`=ifelse(set=="Training set (LUCAS)", stat, `0.95%`)) %>%
        ggplot(aes(stage, stat))  +
        geom_hline(yintercept=1, linetype="dashed") +
        geom_errorbar(aes(ymin=`0.05%`, ymax=`0.95%`,
                          color=set), width=0.1,
                      position=position_dodge(0.5)) +
        geom_point(size=3, shape=21,
                   aes(color=set,
                       fill=set),
                   position=position_dodge(0.5)) +
        theme_classic(base_size=20) +
        theme(legend.position=c(0.5, 0.25),
              legend.text=element_text(size=15),
              panel.background=element_rect(fill="white",
                                            color="gray30"),
              panel.grid=element_blank()) +
        scale_color_manual(values=c("gray30", "steelblue")) +
        scale_fill_manual(values=c("gray30", "steelblue")) +
        ylab("") +
        xlab("") +
        scale_x_discrete(labels=xlabels) +
        scale_y_continuous(labels=scales::percent_format(accuracy=1),
                           breaks=c(0, 0.5, spec_cutoff, 1),
                           limit=c(0, 1)) +
        ##scale_x_discrete(labels=xlabels) +
        guides(fill=FALSE, color=FALSE)
}



rocStageAndHistology <- function(preds, full_model=TRUE){
    if(!full_model){
        preds$score <- preds$score.seq
    } else preds$score <- preds$score.full.lasso2
    rocI <- preds[grepl("Cancer-free|^I$", Stage)][,rocstats(type, score)]
    rocII <- preds[grepl("Cancer-free|^II$", Stage)][,rocstats(type, score)]
    rocIII <- preds[grepl("Cancer-free|^III$", Stage)][,rocstats(type, score)]
    rocIV <- preds[grepl("Cancer-free|^IV$", Stage)][,rocstats(type, score)]

    rocI[,features:="I"]
    rocII[,features:="II"]
    rocIII[,features:="III"]
    rocIV[,features:="IV"]
    rocstage <- rbind(rocI, rocII, rocIII, rocIV)

    preds[grepl("healthy", type), histology:="Cancer-free"]
    roc.adeno <- preds[grepl("Cancer-free|^Adenocarcinoma$", histology)][, rocstats(type, score)]
    roc.sclc <- preds[grepl("Cancer-free|^SCLC$", histology)][,rocstats(type, score)]
    roc.squamous <- preds[grepl("Cancer-free|^Squamous$", histology)][,rocstats(type, score)]
    roc.met <- preds[grepl("Cancer-free|^Metastasis", histology)][,rocstats(type, score)]

    roc.adeno[,features:="Adeno"]
    roc.sclc[,features:="SCLC"]
    roc.squamous[,features:="Squamous"]
    roc.met[,features:="Met"]
    rochistology <- rbind(roc.adeno, roc.sclc, roc.squamous, roc.met)
    rocdat <- bind_rows(rocstage, rochistology)
    rocdat
}

lucasCutoff <- function(desired_specificity, pred){
    roc2 <- cutpointr::roc
    lucas_cutoff <- pred %>%
        roc2(x = score.seq, class = type,
            pos_class="cancer",
            neg_class="healthy",
            direction = ">=") %>%
        mutate(sens=tp/(tp+fn),
               spec=1-fpr) %>%
        filter(spec >= desired_specificity, is.finite(x.sorted)) %>%
        pull(x.sorted) %>%
        min()
    lucas_cutoff
}

fig3_clindat <- function(se){
    clindat <- colData(se) %>%
        as_tibble() %>%
        mutate(lab_id=colnames(se)) %>%
        mutate(na1=ifelse(is.na(score_seq), 1, 0),
               na2=ifelse(is.na(score_no_priorcancer), 1, 0),
               na3=ifelse(is.na(score_hr_smokers), 1, 0),
               nna=na1+na2+na3) %>%
        filter(nna < 3)  %>%
        mutate(histology=as.character(histology),
               stage=as.character(stage),
               histology=ifelse(histology=="Metastasis from another primary cancer",
                                "Lung metastasis", histology)) %>%
        select(lab_id, histology, patient_type, stage, hist_dx)
    clindat
}


fig3_data <- function(se, dat1, dat2, dat3){
    clindat <- fig3_clindat(se)
    histology1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    histology2 <- clindat %>%
        filter(lab_id %in% dat2$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[2])
    histology3 <- clindat %>%
        filter(lab_id %in% dat3$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[3])
    stage1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    stage2 <- clindat %>%
        filter(lab_id %in% dat2$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[2])
    stage3 <- clindat %>%
        filter(lab_id %in% dat3$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[3])
    nocancer1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    nocancer2 <- clindat %>%
        filter(lab_id %in% dat2$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[2])
    nocancer3 <- clindat %>%
        filter(lab_id %in% dat3$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[3])
    category_levels <- c("No biopsy",
                         "Benign nodule",
                         "I", "II", "III", "IV",
                         "Adenocarcinoma", "Squamous",
                         "SCLC", "Lung metastasis")
    clindat3 <- bind_rows(histology1,
                          histology2,
                          histology3,
                          stage1,
                          stage2,
                          stage3,
                          nocancer1,
                          nocancer2,
                          nocancer3) %>%
        mutate(category=factor(category, category_levels))

    combined <- clindat3 %>%
        left_join(dat, by=c("lab_id", "training_model")) %>%
        mutate(training_model=factor(training_model, model_levels))
    hlevels <- category_levels[7:10]
    nclevels <- category_levels[1:2]
    stlevels <- category_levels[3:6]
    combined2 <- combined %>%
        mutate(groups=case_when(category %in% hlevels~"Histology",
                                category %in% nclevels~"Non-cancer",
                                category %in% stlevels~"Cancer stage"),
               groups=factor(groups, c("Non-cancer", "Cancer stage", "Histology")))
    combined2
}

extfig9_data <- function(se, dat1, model_levels){
    clindat <- fig3_clindat(se)
    histology1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    stage1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    nocancer1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    category_levels <- c("No biopsy",
                         "Benign nodule",
                         "I", "II", "III", "IV",
                         "Adenocarcinoma", "Squamous",
                         "SCLC", "Lung metastasis")
    clindat3 <- bind_rows(histology1,
                          stage1,
                          nocancer1) %>%
        mutate(category=factor(category, category_levels))
    combined <- clindat3 %>%
        left_join(dat1, by=c("lab_id", "training_model")) %>%
        mutate(training_model=factor(training_model, model_levels))
    hlevels <- category_levels[7:10]
    nclevels <- category_levels[1:2]
    stlevels <- category_levels[3:6]
    combined2 <- combined %>%
        mutate(groups=case_when(category %in% hlevels~"Histology",
                                category %in% nclevels~"Non-cancer",
                                category %in% stlevels~"Cancer stage"),
               groups=factor(groups, c("Non-cancer", "Cancer stage", "Histology")))
    combined2
}

fig3_roc <- function(se, dat1, dat2, dat3,
                     model_levels){
    category_levels <- c("No biopsy",
                         "Benign nodule",
                         "I", "II", "III", "IV",
                         "Adenocarcinoma", "Squamous",
                         "SCLC", "Lung metastasis")
    clindat <- fig3_clindat(se)
    histology1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    histology2 <- clindat %>%
        filter(lab_id %in% dat2$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[2])
    histology3 <- clindat %>%
        filter(lab_id %in% dat3$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[3])
    stage1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    stage2 <- clindat %>%
        filter(lab_id %in% dat2$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[2])
    stage3 <- clindat %>%
        filter(lab_id %in% dat3$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[3])
    nocancer1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    nocancer2 <- clindat %>%
        filter(lab_id %in% dat2$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[2])
    nocancer3 <- clindat %>%
        filter(lab_id %in% dat3$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[3])
    stage_model1 <- performance_wrapper(stage1, dat1, nocancer1)
    stage_model2 <- performance_wrapper(stage2, dat2, nocancer2)
    stage_model3 <- performance_wrapper(stage3, dat3, nocancer3)
    hist_model1 <- performance_wrapper(histology1, dat1, nocancer1)
    hist_model2 <- performance_wrapper(histology2, dat2, nocancer2)
    hist_model3 <- performance_wrapper(histology3, dat3, nocancer3)
    roc_stage <- select(stage_model1, category, performance)
    roc_stage$performance <- roc_stage$performance %>%
        map2(stage_model2$performance, bind_rows) %>%
        map2(stage_model3$performance, bind_rows)
    roc_hist <- select(hist_model1, category, performance)
    roc_hist$performance <- roc_hist$performance %>%
        map2(hist_model2$performance, bind_rows) %>%
        map2(hist_model3$performance, bind_rows)
    roc_stage2 <- unnest(roc_stage, "performance")
    roc_hist2 <- unnest(roc_hist, "performance")
    roc_strata <- bind_rows(roc_stage2, roc_hist2) %>%
        mutate(category=factor(category, category_levels[3:10]))
    roc_strata
}

extfig9_roc <- function(se, dat1, model_levels){
    category_levels <- c("No biopsy",
                         "Benign nodule",
                         "I", "II", "III", "IV",
                         "Adenocarcinoma", "Squamous",
                         "SCLC", "Lung metastasis")
    clindat <- fig3_clindat(se)
    histology1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(histology)) %>%
        select(lab_id, histology) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    stage1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(!is.na(stage)) %>%
        select(lab_id, stage) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    nocancer1 <- clindat %>%
        filter(lab_id %in% dat1$lab_id) %>%
        filter(is.na(stage)) %>%
        select(lab_id, hist_dx) %>%
        mutate(hist_dx=ifelse(hist_dx=="No baseline cancer", "No biopsy",
                              "Benign nodule")) %>%
        set_colnames(c("lab_id", "category")) %>%
        mutate(training_model=model_levels[1])
    stage_model1 <- performance_wrapper(stage1, dat1, nocancer1)
    hist_model1 <- performance_wrapper(histology1, dat1, nocancer1)
    roc_stage <- select(stage_model1, category, performance)
    ##    roc_stage$performance <- roc_stage$performance %>%
    ##        map2(stage_model2$performance, bind_rows) %>%
    ##        map2(stage_model3$performance, bind_rows)
    roc_hist <- select(hist_model1, category, performance)
##    roc_hist$performance <- roc_hist$performance %>%
##        map2(hist_model2$performance, bind_rows) %>%
##        map2(hist_model3$performance, bind_rows)
    roc_stage2 <- unnest(roc_stage, "performance")
    roc_hist2 <- unnest(roc_hist, "performance")
    roc_strata <- bind_rows(roc_stage2, roc_hist2) %>%
        mutate(category=factor(category, category_levels[3:10]))
    roc_strata
}


add_nocancer <- function(x, nocancer){
    xx <- bind_rows(x, nocancer) %>%
        mutate(class=rep(c("Cancer", "No cancer"),
                         c(nrow(x), nrow(nocancer))),
               class=factor(class, c("No cancer", "Cancer")))
    xx
}
join_scores <- function(x, scores){
    xx <- left_join(x, scores, by=c("lab_id",
                                    "training_model"))
    xx
}
performance <- function(x){
    p <- pROC::roc(x$class, x$score, ci=TRUE)
               ##levels=c("healthy", "cancer"),
               ##ci=TRUE))
    dat <- tibble(sens = rev(p[["sensitivities"]]),
                  spec = rev(p[["specificities"]])) %>%
        mutate(auc=as.numeric(p$auc),
               lower=p$ci[1],
               upper=p$ci[3],
               training_model=x$training_model[1])
    dat
}
performance_wrapper <- function(model, scores, nocancer){
    model <- model %>%
        group_by(category) %>%
        nest() %>%
        arrange(category)
    model$performance <- model$data %>%
        map(add_nocancer, nocancer=nocancer) %>%
        map(join_scores, scores) %>%
        map(performance)
    model
}
