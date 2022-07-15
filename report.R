## ----init1, echo=T, warning=F, message=FALSE, cache=F-------------------------
#load the R packages we need. 
library("pander")
library("coin")
library("knitr")
library("rmarkdown")
library("summarytools")
library("parallel")
library("DescTools")
library("ggplot2")
library("rcompanion")
library("ggpubr")

library("ggforce")
source("statfuns.r")

options(stringsAsFactors=F)
`%nin%` = Negate(`%in%`)

#file = "BoosterVaccineStudy_Database_20220510.csv"
#file = "BoosterVaccineStudy_Database_20220517.csv"
#file = "BoosterVaccineStudy_Database_20220518.csv"
#file = "BoosterVaccineStudy_Database_20220520.csv"
#file = "BoosterVaccineStudy_Database_20220531.csv"
file = "BoosterVaccineStudy_Database_20220604.csv"

#load the booster dataset
#x = read.csv("booster_master_20211030.csv")
#x = read.csv("booster_master_20211102.csv")
x = read.csv(file)

#clean up the header names
#h1 = read.csv("booster_master_20211030.csv", header=F, nrow=1)
#h1 = read.csv("booster_master_20211102.csv", header=F, nrow=1)
h1 = read.csv(file, header=F, nrow=1)
#h1 = make.names(gsub(h1[1,], pat="\\(.*\\)", rep=""))
h1 = make.names(gsub(h1, pat="\\.\\.", rep="."))
h1 = make.names(gsub(h1, pat="\\.$", rep=""))
colnames(x) = h1

#trim whitespace and set to lowercase from every character
for (v in colnames(x)){
    if (is.character(x[,v])){
        x[,v] = tolower(trimws(x[,v]))
        #get rid of wierd white spaces
        x[,v] = trimws(x[, v], whitespace = "[\\h\\v]")
    }
}

#exclude patients from first column
table(x[,"Patients.to.exclude"] %in% c("1", "duplicate"))
#x = x[is.na(x[,"Patients.to.exclude"]),]
x = x[x[,"Patients.to.exclude"]=="",]

#exclude any row without an MRN
x = x[!is.na(x$MRN),]

#fix the spike.ab numbers
fix.num.vars = c(
    "Post.vacination.Spike.Ab.Result..AU.mL",
    "Baseline.Study.Spike.Ab.Result..AU.mL",
    "X4.week.Spike.Ab.Result..AU.mL",
    "Baseline.T.cell.Assay.Result..mIU.mL",
    "X4.week.T.cell.Assay.Result..mIU.mL",
    "IgA",
    "IgM",
    "IgG"
)


for (v in fix.num.vars){
    x[,v] = gsub(x[,v], pat="[<,>]", rep="")
    x[,v] = gsub(x[,v], pat="au/ml", rep="")
    x[,v] = as.numeric(x[,v])
}


#fix the test result vars
fix.res.vars = c(
    "X4.week.Spike.Antibody.Result..Positive.or.Negative"
)
for (v in fix.res.vars){
    x[x[,v] %nin% c("positive", "negative"),v] = NA
}

#fix the sinai variables
x[,"Sinai.neutralization.assay.4.week.WT.ID50"] = as.numeric(x[,"Sinai.neutralization.assay.4.week.WT.ID50"])
x[,"Sinai.neutralization.assay.4.week.WT.omicron.ID50"] = as.numeric(x[,"Sinai.neutralization.assay.4.week.WT.omicron.ID50" ])



## ----results="asis"-----------------------------------------------------------
#run a quick summary for each variable
for (v in colnames(x)){
    cat("##", v, "\n")
    cat("\n")
    if (length(unique(x[,v])) < 40 & length(unique(x[,v])) > 0){
        cat(pander(table(x[,v], useNA="ifany")))
        cat("\n")
    }else{
        cat(pander(summary(x[,v])))
        cat("\n")
    }
}



## ---- warning=F---------------------------------------------------------------
#is the time between less than median of cohort?
t1 = x[,"Time.from.vaccination.completion.to.Booster..months"]
median(t1)
t2 = t1 <= median(t1)
x[,"time.vac.to.booster.lessThanMedian"] = t2



#Highly immune suppressed (AZ + BF + BN + BP + BJ) vs. all others (AZ thru CH excluding AZ, BF, BN, BP, BJ)
his.vars = c(
    "chemo",
    "SCT",
    "CD20",
    "CAR.T",
    "BTK"
)
#is the patient highly immune suppresed?
#do they have a 1 in any of the his.vars?
x$is.highly.immune.suppressed = rowSums(x[,his.vars]) > 0

#stratify by 30 and 90
x$time.from.last.cancer.therapy.to.booster.lessthan30days = x$Time.from.last.cancer.therapy.to.booster..days <= 30
x$time.from.last.cancer.therapy.to.booster.lessthan90days = x$Time.from.last.cancer.therapy.to.booster..days <= 90

#split by 6month
x$time.since.last.cd20.prior.booster.lessthan6month = x$Time.since.last.CD20.prior.to.booster..days <= 6*30

#strat by 1 year
x$time.since.last.transplant.prior.booster.lessthan1y = x$Time.since.last.transplant.prior.to.booster..months <= 12


#make a variable describing bcell directed therapy vs others
#If a patient has anything in BL+BP+BR they get a 1.
#If a patient has anything in BB through CJ (and nothing in BL+BP+BR), they get a 0
#If a patient has nothing in BB through CJ they are ignored?
bcell.therapies = c("BTK", "CD20", "CAR.T")
other.therapies = colnames(x)[which(colnames(x) == "chemo"):which(colnames(x) == "PARP")]
other.therapies = other.therapies[other.therapies %nin% bcell.therapies]

therapy = NA
therapy[rowSums(x[,other.therapies]) > 0] = "other"
therapy[rowSums(x[,bcell.therapies]) > 0] = "bcell"
x$therapy = therapy


x$BTK = as.logical(x$BTK)
x$SCT = as.logical(x$SCT)
x$CD20 = as.logical(x$CD20)
x$BTK_CD20 = x$BTK | x$CD20
#stratifying variables for analysis A and B
vars.AB = c(
    "Cancer.Diagnosis",
    "Malignancy.category..Solid.or.Liquid",
    "SOLID.LYMPHOID.MYELOID",
    "Previous.Type.Vaccine.Given",
    "time.vac.to.booster.lessThanMedian",
    "Type.of.Booster.Given",
    "is.highly.immune.suppressed",
    "time.from.last.cancer.therapy.to.booster.lessthan30days",
    "time.from.last.cancer.therapy.to.booster.lessthan90days",
    "time.since.last.cd20.prior.booster.lessthan6month",
    "time.since.last.transplant.prior.booster.lessthan1y",
    "Prior.COVID.infection...yes.no",
    "On.cancer.therapy.at.the.time.of.booster...yes.no",
    "SCT",
    "BTK",
    "CD20",
    "BTK_CD20",
    "therapy"
)

vars.CD = c(
    "SCT",
    "BTK",
    "CD20",
    "BTK_CD20"
)

vars.AB %in% colnames(x)
#x[,vars.AB]



#summarize a variable include # of NAs
mysummary <- function(x){
    c(summary(x), n = sum(!is.na(x)), na=sum(is.na(x)))
}

#For paired categorical variables
#split 2x2 table(var1, var2) by levels of var3
#and run statistical tests
var1 = "Baseline.Study.Spike.Ab.Result..Positive.or.Negative"
var2 = "X4.week.Spike.Antibody.Result..Positive.or.Negative"
var3 = "Malignancy.category..Solid.or.Liquid"
runStratifiedPairedTables <- function(x, var1, var2, var3){
    tab0 = table(x[,var1], x[,var2], x[,var3], dnn=c(var1, var2, var3))

    #run the homogeneity of stratum effects test
    tab2kx2 = c()
    for(i in 1:dim(tab0)[3]){
        tab2kx2 = rbind(tab2kx2, tab0[,,i])    
    }
    res = HSE(tab2kx2)
    cat(paste0("HSE test p-value: ", res$p, "\n\n"))

    #run mcnemar on every level of stratifying variable
    levs = unique(x[,var3])
    levs = levs[!is.na(levs)]
    for (v3 in levs){
        cat(paste0("####", v3, "\n"))
        #make sure vars are factors so levels aren't dropped
        x2 = x[,c(var1, var2, var3)]
        x2[,var1] = factor(x2[,var1])
        x2[,var2] = factor(x2[,var2])
        x2[,var3] = factor(x2[,var3])
        x3 = x2[x2[,var3] == v3,]

        #split into 2x2 table
        ctab1 = ctable(x3[,var1], x3[,var2], dnn=c(var1, var2), useNA="no")
        cat(pander(ctab1[[1]]))

        tab1 = table(x3[,var1], x3[,var2], dnn=c(var1, var2))
        #cat(pander(tab1))

        #run mcnemars test 
        cat(pander(mcnemar.test(tab1)))
        cat(paste0("\n"))
    }
}
#runStratifiedPairedTables(x, var1, var2, var3)



runPropTest2 <- function(x, var1, var2){
    #skip everything if its the same variable
    if (var1 == var2){
        return (F)
    }


    cat(paste0("##", var1, " by ", var2, "\n"))
    #make sure vars are factors so levels aren't dropped
    x2 = x[,c(var1, var2)]
    x2[,var1] = factor(x2[,var1])
    x2[,var2] = factor(x2[,var2])

    #split into 2x2 table
    ctab1 = ctable(x2[,var1], x2[,var2], dnn=c(var1, var2), useNA="no")
    cat(pander(ctab1[[1]]))

    tab1 = table(x2[,var1], x2[,var2], dnn=c(var1, var2))

    try({
        res = fisher.test(tab1)
        if (res$p.value < .1){

            cat(paste0("### * \n"))
        }
        #run fisher's exact test 
        cat(pander(fisher.test(tab1)))
        cat(paste0("\n"))
        knitr::normal_print(fisher.test(tab1))
        cat(paste0("\n"))
    }, silent=T)
}



#simple fisher test of proportions
runPropTest <- function(x, var1, var2){
    cat(paste0("##", var1, " by ", var2, "\n"))
    #make sure vars are factors so levels aren't dropped
    x2 = x[,c(var1, var2)]
    x2[,var1] = factor(x2[,var1])
    x2[,var2] = factor(x2[,var2])

    #split into 2x2 table
    ctab1 = ctable(x2[,var1], x2[,var2], dnn=c(var1, var2), useNA="no")
    cat(pander(ctab1[[1]]))

    tab1 = table(x2[,var1], x2[,var2], dnn=c(var1, var2))

    #run fisher test 
    cat(pander(fisher.test(tab1)))
    cat(paste0("\n"))
}


#run kruskal wallis on a numeric vs categorical variable
#var1 the name of the numeric var
#var2 the name of the categorical var
runDiffAssociationTest <- function(x, var1, var2){
    cat(paste0("##", var1, " by ", var2, "\n"))
    #make sure vars are factors so levels aren't dropped
    x2 = x[,c(var1, var2)]
    x2[,var2] = factor(x2[,var2])

    vals = x2[,var1]
    group = x2[,var2]
    cat(pander(kruskal.test(vals~group)))
    cat(pander(by(vals, group, FUN=summary)))
}


#run kruskal wallis on a numeric vs categorical variable
#var1 the name of the numeric var
#var2 the name of the categorical var
runDiffAssociationTest2 <- function(x, var1, var2){
    #skip everything if its the same variable
    if (var1 == var2){
        return(F)
    }

    cat(paste0("##", var1, " by ", var2, "\n"))
    #make sure vars are factors so levels aren't dropped
    x2 = x[,c(var1, var2)]
    x2[,var2] = factor(x2[,var2])

    vals = x2[,var1]
    group = x2[,var2]
    try({
        res = kruskal.test(vals~group)
        if (res$p.value < .1){

            cat(paste0("### * \n"))
        }
        cat(pander(kruskal.test(vals~group)))
        cat(pander(by(vals, group, FUN=mysummary)))
    }, silent=T)
}


#look at the counts of patients that hvae higher/lower titers at 2nd test(t2) compared to 1st(t1)
#tabulate with respect to a stratifying variable v
#run fisher test to check for association
runTiterTables <- function(t1, t2, v){
    tab1 = table(t1.gte.t2=t1>=t2, v)    
    cat(pander(ftable(tab1)))
    cat(pander(fisher.test(tab1)))
}


#look for correclation between two variables
runCorTest <- function(x, var1, var2){

    cat(paste0("##", var1, " by ", var2, "\n"))
    
    vals1 = x[,var1]
    vals2 = x[,var2]

    #simple scatter plot
    df1 = data.frame(vals1, vals2)
    #p1 = ggplot(df1, aes(x = vals1, y = vals2)) + geom_point() +
    #    xlab(var1) + ylab(var2)
    #print(p1)
    #cat(paste0("\n"))

    #cat(pander(cor.test(vals1, vals2, method="kendall")))
    cat(pander(cor.test(vals1, vals2, method="spearman")))

    #if significant show a star
    try({
        res = cor.test(vals1, vals2, method="spearman")
            if (res$p.value < .10){
                cat(paste0("### *\n"))
    }}, silent=T)
    cat(paste0("\n"))
}



#runStratifiedPairedVals
#check for association between the 
#difference of paired values(var1, var2)
#at various levels of var3
runStratifiedPairedVals <- function(x, var1, var2, var3){ 
    t1 = x[,var1]
    t2 = x[,var2]
    diff = t2 - t1
    group = x[,var3]   

    cat(pander(kruskal.test(diff~group)))
    cat(pander(by(diff, group, FUN=summary)))

    
}

plotTitersByVar <- function(t1, t2, group, ylim=NA, timepoint = "4w", ylab="titers"){
    df1 = data.frame(rbind(
        data.frame(time=0, samp = 1:length(t1), titer=t1, group=group),
        data.frame(time=timepoint, samp = 1:length(t1), titer=t2, group=group)
    ))
    df1$time = as.factor(df1$time) 

    #take out NA valueso
    #df1 = df1[!is.na(df1$group),]

    p1 = ggplot(df1, aes(x=time, y=titer, group=samp, color=group)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.25) + ylab(ylab) + theme_classic()
    if (!is.na(ylim)){
        p1 = p1 + ylim(ylim)
    }

    #calculate CIs
    #https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
    ##add mean trends for each level
    df2 = lapply(unique(group), function(g){
        alpha=0.05
        t=qt((1-alpha)/2 + .5, sum(group==g)-1)   # tend to 1.96 if sample size is big enough
        rbind(
            data.frame(time=0, samp = g, titer=mean(t1[group==g], na.rm=T), titer.ci = sd(t1[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g),
            data.frame(time=timepoint, samp = g, titer=mean(t2[group==g], na.rm=T), titer.ci = sd(t2[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g)
        )
    })
    df2 = do.call(rbind, df2)
    df2$time = as.factor(df2$time) 
    p2 = p1 + geom_line(alpha=1, data=df2, size=1.5) 
    p3 = p2 + geom_errorbar(aes(ymin = titer - titer.ci, ymax = titer + titer.ci), data=df2, size=.75, width=.15)
    print(p2)
    print(p3)


    df3 = data.frame(rbind(
        data.frame(time=0, titer=t1, group=group),
        data.frame(time=timepoint, titer=t2, group=group)
    ))

    df3 = df3[!is.na(df3$group),]
    #b1 = ggplot(df3, aes(x=time, y=titer, fill=group)) +  geom_boxplot(outlier.shape=NA) +  geom_point(position=position_jitterdodge(jitter.width=.1), color="black") + ylab(ylab) + theme_classic()

    b2 = ggplot(df3, aes(x=time, y=titer, fill=group)) +
        stat_summary(position=position_dodge(), fun.data ="mean_cl_boot", geom="crossbar")  + 
        geom_point(position=position_jitterdodge(jitter.width=.1), color="black") + 
        ylab(ylab) + 
        theme_classic()

    print(b2)
    print(b2+scale_y_log10())


}





## -----------------------------------------------------------------------------
#prep variables

#Column W- assess how many pts had T-cell results evaluable (positive/neg/borderline). Exclude invalid and #NA
table(x[,"Baseline.T.cell.Assay.Result..pos..bor..neg..or.inv"])
tcell.res = x[,"Baseline.T.cell.Assay.Result..pos..bor..neg..or.inv"]
tcell.res[grepl(tcell.res, pat="bor")] = "bor"
tcell.res[tcell.res %nin% c("bor", "pos", "neg")] = NA
table(tcell.res)
#treat borderline cases as negative
tcell.res[tcell.res == "bor"] = "neg"
x$tcell.res = tcell.res



#make a new variable for age > 65
x$age.over65 = x$Age > 65

#Among the above we'd like to get their 
#median age (col E), 
#male vs female (column F), 
#race (col E) 
#solid vs liquid (col AO) , 
#cancer diagnosis (col H), 
#previous vaccine type (col I) 
#baseline spike ab positive vs negative (col U) 
#and prior treatments received (cols BA thru CI)
tcell.vars = c(
    "Age",
    "age.over65",
    "Sex..Male.Female",
    "Race..White..Black..Hispanic..Asian..Other",
    "Malignancy.category..Solid.or.Liquid",
    "SOLID.LYMPHOID.MYELOID",
    "Cancer.Diagnosis",
    "Previous.Type.Vaccine.Given",
    "Baseline.Study.Spike.Ab.Result..Positive.or.Negative",
    "chemo",
    "surgery",
    "ADC",
    "Trial",
    "Imid",
    "Proteasome",
    "SCT",
    "TGF",
    "XRT",
    "BCL2",
    "BTK",
    "TKI",
    "Immunotherapy",
    "hormone",
    "CD20",
    "CD38",
    "CAR.T",
    "CD19",
    "CD30",
    "Steroids",
    "EGFR",
    "VEGF",
    "CDK4.6",
    "HER2",
    "mTOR",
    "IL.6",
    "Cytokine",
    "Anti.viral",
    "IVIG",
    "AR.targeted",
    "supportive.care",
    "none",
    "PI3Ki",
    "SLAMF7",
    "PARP",
    "therapy"
)




#Assess if booster vaccine led to change in T-cell activity from negative to positive at 4 weeks (col AH)

table(x[,"X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv"])
tcell.4w.res = x[,"X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv"]
tcell.4w.res[grepl(tcell.4w.res, pat="bor")] = "bor"
table(tcell.4w.res)
tcell.4w.res[tcell.4w.res %nin% c("bor", "pos", "neg")] = NA
table(tcell.4w.res)

#treat borderline cases as neg
tcell.4w.res[tcell.4w.res == "bor"] = "neg"
x$tcell.4w.res = tcell.4w.res



## ---- results="asis"----------------------------------------------------------

for (v in tcell.vars){
    if (length(unique(x[,v])) < 30){
        runPropTest2(x, "tcell.res", v)
    }else{
        runDiffAssociationTest2(x, v, "tcell.res")
    }
}



## -----------------------------------------------------------------------------


#might want to treat borderline patients as either negative or positive
#x$tcell.res[x$tcell.res == "bor"] = "neg"
#x$tcell.4w.res[x$tcell.4w.res == "bor"] = "neg"

tab1 = table(x[,c("tcell.res", "tcell.4w.res")], useNA="no")
print(tab1)
pander(tab1)

pander(mcnemar.test(tab1))


#Assess if this change was influenced by 
#age (E), 
#solid vs liquid (AO), 
#cancer diagnosis (H) 
#and prior therapies (BA thru CI)


## ---- results="asis"----------------------------------------------------------

for (v in tcell.vars[-c(1,3, 4,7, 8)]){
    cat(paste0("###", v, "\n"))
    runStratifiedPairedTables(x,"tcell.res", "tcell.4w.res", v)
    cat("\n\n")

}


## -----------------------------------------------------------------------------

x.6mfu = x[which(x[,"include.in.6.month.follow.up.cohort"] == 1),]
dim(x.6mfu)


#there are 3 times points baseline, 4week, 6month

#columns V, AE, DG and DM
#Baseline.Study.Spike.Ab.Result..AU.mL
#X4.week.Spike.Ab.Result..AU.mL
#Pre.4th.dose.spike.ab.value
#X5.6.Month.Study.Spike.Ab.Result..AU.mL
#Evusheld.Regencov.sotrivimab
#X4th.dose

tps = c(1, 30, 180)

t1 = x.6mfu[, "Baseline.Study.Spike.Ab.Result..AU.mL"]
t2 = x.6mfu[, "X4.week.Spike.Ab.Result..AU.mL"]
t3 = x.6mfu[, "X5.6.Month.Study.Spike.Ab.Result..AU.mL"]
#fix values
t3 = as.numeric(gsub(t3, pat="<|>|,", rep=""))
evu.reg.sot = x.6mfu[,"Evusheld.Regencov.sotrivimab"]
group = as.factor(evu.reg.sot)

#dose
dose4 = as.factor(x.6mfu[, "X4th.dose"])

df1 = data.frame(rbind(
    data.frame(time=tps[1], samp = 1:length(t1), titer=t1, group=group, dose4=dose4),
    data.frame(time=tps[2], samp = 1:length(t1), titer=t2, group=group, dose4=dose4),
    data.frame(time=tps[3], samp = 1:length(t1), titer=t3, group=group, dose4=dose4)
))
#df1$time = factor(df1$time, levels=tps) 
df1$treatment = NA
df1$treatment[df1$group != ""] = "evu.reg.sot"
df1$treatment[df1$group == ""] = "none"
treatment = df1$treatment


asis_output("##Spike.ab Evusheld.Regencov.sotrivimab\n")
p1 = ggplot(df1, aes(x=time, y=titer, group=samp, color=group)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.5) + theme_classic()
p1

asis_output("##Spike.ab lumped(Evusheld.Regencov.sotrivimab)\n")
p1a = ggplot(df1, aes(x=time, y=titer, group=samp, color=treatment)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.5) + theme_classic()

#medians by none vs evu.reg.sot
by(t3-t2, group, median)
#medians by none vs lumped(evu.reg.sot)
by(t3-t2, evu.reg.sot=="", median)

##medians by none percentages
#by(((t2-t3)/t2)*100, evu.reg.sot=="", median)

#medians of baseline
by(t1, evu.reg.sot=="", median)
#4week
by(t2, evu.reg.sot=="", median)
#6mo
by(t3, evu.reg.sot=="", median)

#counts of treatment types
table(evu.reg.sot)
p1a

    ##########################
    #add a plot with error bars
    df2 = lapply(unique(treatment), function(g){
        alpha=0.05
        t=qt((1-alpha)/2 + .5, sum(treatment==g)-1)   # tend to 1.96 if sample size is big enough
        rbind(
            data.frame(time=1, samp = g, titer=mean(t1[treatment==g], na.rm=T), titer.ci = sd(t1[treatment==g], na.rm=T)/sqrt(sum(treatment==g))*t, treatment=g),
            data.frame(time=30, samp = g, titer=mean(t2[treatment==g], na.rm=T), titer.ci = sd(t2[treatment==g], na.rm=T)/sqrt(sum(treatment==g))*t, treatment=g),
            data.frame(time=180, samp = g, titer=mean(t3[treatment==g], na.rm=T), titer.ci = sd(t3[treatment==g], na.rm=T)/sqrt(sum(treatment==g))*t, treatment=g)
        )
    })
    df2 = do.call(rbind, df2)
    df3 = df1
    levs = unique(c(df3$samp, df2$samp))
    df3$samp = factor(df3$samp, levels = levs)
    df2$samp = factor(df2$samp, levels = levs)
    p1b = ggplot(df3, aes(x=time, y=titer, group=samp, color=treatment)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.5) + theme_classic()
    p2 = p1b + geom_line(alpha=1, data=df2, size=1.5)
    p3 = p2 + geom_errorbar(aes(ymin = titer - titer.ci, ymax = titer + titer.ci), data=df2, size=.75, width=.15)
    p3
    ###########################
    
asis_output("##Spike.ab X4th.dose\n")
p2 = ggplot(df1, aes(x=time, y=titer, group=samp, color=dose4)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.5) + theme_classic()
p2


asis_output("##Spike.ab\n")
p1 = ggplot(df1, aes(x=time, y=titer, group=samp)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.25, color="blue") + theme_classic()
p1



asis_output("##difference between 6month and 4W spike.ab.titer split by solid/liquid\n")

#difference between 6month and 4W
tdiff = t3-t2
kruskal.test(tdiff ~ x.6mfu[,"Malignancy.category..Solid.or.Liquid"])
pander(by(tdiff, x.6mfu[,"Malignancy.category..Solid.or.Liquid"] , mysummary))





#asis_output("##6month and 4W spike.ab.titer split by highly immune suppressed\n")
#kruskal.test(tdiff ~ x.6mfu[,"is.highly.immune.suppressed"])
#pander(by(tdiff, x.6mfu[,"is.highly.immune.suppressed"] , mysummary))
#plotTitersByVar(t2, t3, x.6mfu[, "is.highly.immune.suppressed"])
#
#asis_output("##6month and 4W spike.ab.titer split by immunotherapy\n")
#kruskal.test(tdiff ~ as.factor(x.6mfu[, "Immunotherapy"]))
#pander(by(tdiff, as.factor(x.6mfu[, "Immunotherapy"]) , mysummary))
#plotTitersByVar(t2, t3, as.factor(x.6mfu[, "Immunotherapy"]))




## -----------------------------------------------------------------------------

spike1 = x[, "X4.week.Spike.Ab.Result..AU.mL"]
sinai1 = x[,"Sinai.neutralization.assay.4.week.WT.ID50"] = as.numeric(x[,"Sinai.neutralization.assay.4.week.WT.ID50"])
sinai2 = x[,"Sinai.neutralization.assay.4.week.WT.omicron.ID50"] = as.numeric(x[,"Sinai.neutralization.assay.4.week.WT.omicron.ID50" ])

asis_output("##X4.week.Spike.Ab.Result..AU.mL vs Sinai.neutralization.assay.4.week.WT.ID50\n")

cor(spike1, sinai1, method="spearman", use="pairwise.complete.obs")
cor.test(spike1, sinai1, method="spearman", use="pairwise.complete.obs")
plot(spike1, sinai1)

asis_output("##X4.week.Spike.Ab.Result..AU.mL vs Sinai.neutralization.assay.4.week.WT.omicron.ID50\n")

cor(spike1, sinai2, method="spearman", use="pairwise.complete.obs")
cor.test(spike1, sinai2, method="spearman", use="pairwise.complete.obs")
plot(spike1, sinai2)



## -----------------------------------------------------------------------------

#run Mcnemar's test on the paired test results

tab1 = table(x[,c("Baseline.Study.Spike.Ab.Result..Positive.or.Negative", "X4.week.Spike.Antibody.Result..Positive.or.Negative")])
print(tab1)
#pander(tab1)

ctab1 = ctable(x[,c("Baseline.Study.Spike.Ab.Result..Positive.or.Negative", "X4.week.Spike.Antibody.Result..Positive.or.Negative")], useNA="no")
pander(ctab1[[1]])

pander(mcnemar.test(tab1))
   
  


## ---- results="asis"----------------------------------------------------------

for (v in vars.AB){
    cat(paste0("###", v, "\n"))
    runStratifiedPairedTables(x,"Baseline.Study.Spike.Ab.Result..Positive.or.Negative", "X4.week.Spike.Antibody.Result..Positive.or.Negative", v)
    cat("\n\n")

}


## -----------------------------------------------------------------------------

#Baseline.Study.Spike.Ab.Result..AU.mL
#X4.week.Spike.Ab.Result..AU.mL
t1 = x[,"Baseline.Study.Spike.Ab.Result..AU.mL"]
t2 = x[,"X4.week.Spike.Ab.Result..AU.mL"]

wilcox.test(t2-t1)
t.test(t2-t1)


#what are the patients that are missing values?
pander(x[is.na(t2),c("Name", "MRN", "X4.week.Spike.Ab.Result..AU.mL")])



df1 = data.frame(rbind(
    cbind(time=0, samp = 1:length(t1), titer=t1),
    cbind(time=4, samp = 1:length(t1), titer=t2)
))
df1$time = factor(df1$time)
ggplot(df1, aes(x=time, y=titer, color = time, group=samp)) + geom_jitter(width=.05) + geom_line(alpha=.15, color="black") + theme_classic()








## ---- results="asis", warning=F, cache=F--------------------------------------


for (v in vars.AB){

    cat(paste0("###", v, "\n"))


    t1 = x[,"Baseline.Study.Spike.Ab.Result..AU.mL"]
    t2 = x[,"X4.week.Spike.Ab.Result..AU.mL"]
    group = x[,v]

    group[group == "n/a"] = NA

    plotTitersByVar(t1, t2, group)

    #if there are NAs, then also make a plot without them
    if (sum(is.na(group))>0){
        ix = which(!is.na(group))
        plotTitersByVar(t1[ix], t2[ix], group[ix])
    }
    #plotTitersByVar(log10(t1), log10(t2), group, ylab="log10(titers)")
    cat("\n\n")

    cat("####Baseline.Study.Spike.Ab.Result..AU.mL\n") 
    cat(pander(by(t1, group, MeanCI, na.rm=T)))
    cat("####X4.week.Spike.Ab.Result..AU.mL\n") 
    cat(pander(by(t2, group, MeanCI, na.rm=T)))



    runTiterTables(t1, t2, group)
    cat("\n\n")

    runStratifiedPairedVals(x,"Baseline.Study.Spike.Ab.Result..AU.mL", "X4.week.Spike.Ab.Result..AU.mL", v)
    cat("\n\n")

}


## -----------------------------------------------------------------------------

df2 = x[,c("Baseline.Study.Spike.Ab.Result..AU.mL", "X4.week.Spike.Ab.Result..AU.mL", "Baseline.T.cell.Assay.Result..mIU.mL", "X4.week.T.cell.Assay.Result..mIU.mL")]



asis_output("##spike ab baseline vs tcell assay baseline\n")

cor.test(x[,"Baseline.Study.Spike.Ab.Result..AU.mL"], x[,"Baseline.T.cell.Assay.Result..mIU.mL"], method="spearman", use="pairwise.complete.obs")
#plot(x$spike.4, x$X4.WEEK.B.1.1.529.ID50)
#ggplot(x, aes(x=spike.4, y=X4.WEEK.B.1.1.529.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(df2, x = "Baseline.Study.Spike.Ab.Result..AU.mL", y = "Baseline.T.cell.Assay.Result..mIU.mL",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman", label.y.npc="bottom", label.x.npc = "middle")


asis_output("##spike ab 4 week vs tcell assay 4 week\n")

cor.test(x[,"X4.week.Spike.Ab.Result..AU.mL"], x[,"X4.week.T.cell.Assay.Result..mIU.mL"], method="spearman", use="pairwise.complete.obs")
#plot(x$spike.4, x$X4.WEEK.B.1.1.529.ID50)
#ggplot(x, aes(x=spike.4, y=X4.WEEK.B.1.1.529.ID50)) + geom_point() + geom_smooth(method=lm, se=T)
ggscatter(df2, x = "X4.week.Spike.Ab.Result..AU.mL", y = "X4.week.T.cell.Assay.Result..mIU.mL",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman", label.y.npc="bottom", label.x.npc = "middle")



asis_output("##combined spike ab vs tcell assay\n")
df3 = rbind(
    data.frame(ab.spike = x[,"Baseline.Study.Spike.Ab.Result..AU.mL"], tcell = x[,"Baseline.T.cell.Assay.Result..mIU.mL"], time="baseline"),
    data.frame(ab.spike = x[,"X4.week.Spike.Ab.Result..AU.mL"], tcell = x[,"X4.week.T.cell.Assay.Result..mIU.mL"], time="four week")
)

ggplot(df3, aes(x=ab.spike, y=tcell, color=time)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    stat_cor(method="spearman", label.y.npc = "bottom",  label.x.npc = "middle", show.legend = FALSE) +
    #coord_cartesian(xlim =c(0, 1000), ylim= c(0,200)) +
    #facet_zoom(
    #    xlim =c(50, 1000),
    #    ylim= c(0,2000), 
    #    horizontal=F,
    #    shrink=T,
    #    show.area=T,
    #    zoom.size=.5
    #)+
    #theme_classic()
    theme_bw() +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme()







## -----------------------------------------------------------------------------
#restrict ourselves to the 35 patient that have baseline spike negative
x.baseSpikeNeg = x[x$Baseline.Study.Spike.Ab.Result..Positive.or.Negative == "negative",]

dim(x.baseSpikeNeg)



## ---- results="asis"----------------------------------------------------------

runDiffAssociationTest2(x.baseSpikeNeg, "IgA", "X4.week.Spike.Antibody.Result..Positive.or.Negative")
runDiffAssociationTest2(x.baseSpikeNeg, "IgM", "X4.week.Spike.Antibody.Result..Positive.or.Negative")
runDiffAssociationTest2(x.baseSpikeNeg, "IgG", "X4.week.Spike.Antibody.Result..Positive.or.Negative")



## -----------------------------------------------------------------------------
    #combine the two datasets

    #load the mix & match dataset
    x.mnm = read.csv("../covidBooster4th/Mix_and_Match_master_20220510.csv")
    #take out 
    x.mnm = x.mnm[is.na(x.mnm[,1]),]

    #exclude any row without an MRN
    x.mnm = x.mnm[!is.na(x.mnm$MRN),]

    dat1 = x.mnm[,c("Outcome", "Baseline.IgG", "Baseline.IGA", "Baseline.IgM")]
    colnames(dat1) = c("outcome", "igg", "iga", "igm")
    dat1$outcome[dat1$outcome == "Non-responder"] = "negative"
    dat1$outcome[dat1$outcome == "Responder"] = "positive"
    dat2 = x.baseSpikeNeg[,c("X4.week.Spike.Antibody.Result..Positive.or.Negative", "IgG", "IgA", "IgM")]
    colnames(dat2) = c("outcome", "igg", "iga", "igm")
    dat3 = rbind(dat1, dat2)   



## ---- results="asis"----------------------------------------------------------

runDiffAssociationTest2(dat3, "iga", "outcome")
runDiffAssociationTest2(dat3, "igm", "outcome")
runDiffAssociationTest2(dat3, "igg", "outcome")



## -----------------------------------------------------------------------------
#fix the new variables
#x[,"Baseline.Study.Spike.Ab.Result..Positive.or.Negative"]
#x[,"X4.week.Spike.Antibody.Result..Positive.or.Negative"]

fixDectec <- function(a){
    b = gsub(a, pat="dectected", rep = "detected")
    b[b == 0] = NA
    b[b == "#n/a" ] = NA
    b
}

x[,"Baseline.Nabs"] = fixDectec(x[,"Baseline.Nabs"])
x[,"X4wk.Nabs"] = fixDectec(x[,"X4wk.Nabs"])

#x[,"Baseline.Study.Spike.Ab.Result..AU.mL"]
#x[,"X4.week.Spike.Ab.Result..AU.mL"]

x$signal.inhibition.baseline = as.numeric(x[,"Baseline...signal.inhibition"])
x$signal.inhibition.4w = as.numeric(x[,"X4.week...signal.inhibition"])





## ---- results="asis"----------------------------------------------------------

runPropTest2(x, "Baseline.Nabs", "Baseline.Study.Spike.Ab.Result..Positive.or.Negative")
runPropTest2(x, "X4wk.Nabs", "X4.week.Spike.Antibody.Result..Positive.or.Negative")

runCorTest(x, "Baseline.Study.Spike.Ab.Result..AU.mL", "signal.inhibition.baseline")



p3 = ggscatter(x, x = "Baseline.Study.Spike.Ab.Result..AU.mL", y = "signal.inhibition.baseline",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman", label.y.npc="bottom", label.x.npc = "middle")
print(p3)

runCorTest(x, "X4.week.Spike.Ab.Result..AU.mL", "signal.inhibition.4w")

p4 = ggscatter(x, x = "X4.week.Spike.Ab.Result..AU.mL", y = "signal.inhibition.4w",
    add = "reg.line",  # Add regressin line
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    conf.int = TRUE # Add confidence interval
    ) + stat_cor(method = "spearman", label.y.npc="bottom", label.x.npc = "middle")
print(p4)





## -----------------------------------------------------------------------------


asis_output("##combined spike ab vs signal inhibition\n")
df3 = rbind(
    data.frame(ab.spike = x[,"Baseline.Study.Spike.Ab.Result..AU.mL"], sig.inhib = x[,"signal.inhibition.baseline"], time="baseline"),
    data.frame(ab.spike = x[,"X4.week.Spike.Ab.Result..AU.mL"], sig.inhib = x[,"signal.inhibition.4w"], time="four week")
)

ggplot(df3, aes(x=ab.spike, y=sig.inhib, color=time)) +
    geom_point() + 
    geom_smooth(method = "lm") +
    stat_cor(method="spearman", label.y.npc = "bottom",  label.x.npc = "middle", show.legend = FALSE) +
    #coord_cartesian(xlim =c(0, 1000), ylim= c(0,200)) +
    #facet_zoom(
    #    xlim =c(50, 1000),
    #    ylim= c(0,2000), 
    #    horizontal=F,
    #    shrink=T,
    #    show.area=T,
    #    zoom.size=.5
    #)+
    #theme_classic()
    theme_bw() +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme()






## -----------------------------------------------------------------------------

knit_exit()


## -----------------------------------------------------------------------------

#restrict dataset to patients that are negative at baseline
x.c = x[x$Baseline.Study.Spike.Ab.Result..Positive.or.Negative == "negative",]



## ---- results="asis"----------------------------------------------------------

for (v in vars.CD){
    runPropTest(x.c, "X4.week.Spike.Antibody.Result..Positive.or.Negative", v)
}



## -----------------------------------------------------------------------------

#looking at all levels
tab1 = table(baseline=x.c$Baseline.T.cell.Assay.Result..pos..bor..neg..or.inv, week4=x.c$X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv)
print(tab1)
pander(tab1)
pander(fisher.test(tab1))


#restricting to just the pos/neg levels
tcell.base = x.c$Baseline.T.cell.Assay.Result..pos..bor..neg..or.inv
tcell.4w = x.c$X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv
tcell.base[tcell.base %nin% c("pos", "neg")] = NA
tcell.4w[tcell.4w %nin% c("pos", "neg")] = NA

tab2 = table(tcell.base, tcell.4w)
print(tab2)
pander(tab2)
pander(fisher.test(tab2))



## -----------------------------------------------------------------------------

t.assay.base = x.c$Baseline.T.cell.Assay.Result..mIU.mL
t.assay.4w = x.c$X4.week.T.cell.Assay.Result..mIU.mL


cor.test(t.assay.base, t.assay.4w)


f = 5
plot(jitter(t.assay.base, factor=f), jitter(t.assay.4w, factor=f))
cbind(mrn=x.c$MRN,  t.assay.base,  t.assay.4w)




## -----------------------------------------------------------------------------

spike.res.4w = x.c$X4.week.Spike.Antibody.Result..Positive.or.Negative
tcell.4w = x.c$X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv

tab1 = table(spike.res.4w, tcell.4w)
print(tab1)
pander(tab1)
pander(fisher.test(tab1))

#restrict to just the pos/neg values
tcell.4w[tcell.4w %nin% c("pos", "neg")] = NA
tab2 = table(spike.res.4w, tcell.4w)
print(tab2)
pander(tab2)
pander(fisher.test(tab2))





## -----------------------------------------------------------------------------

#restrict dataset to patients that are below 50 at baseline
x.d = x[x$Baseline.Study.Spike.Ab.Result..AU.mL <= 50,]



## ---- results="asis"----------------------------------------------------------

for (v in vars.CD){
    runDiffAssociationTest(x.d, var1="X4.week.Spike.Ab.Result..AU.mL", v)
}



## -----------------------------------------------------------------------------

#looking at all levels
tab1 = table(baseline=x.d$Baseline.T.cell.Assay.Result..pos..bor..neg..or.inv, week4=x.d$X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv)
print(tab1)
pander(tab1)
pander(fisher.test(tab1))


#restricting to just the pos/neg levels
tcell.base = x.d$Baseline.T.cell.Assay.Result..pos..bor..neg..or.inv
tcell.4w = x.d$X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv
tcell.base[tcell.base %nin% c("pos", "neg")] = NA
tcell.4w[tcell.4w %nin% c("pos", "neg")] = NA

tab2 = table(tcell.base, tcell.4w)
print(tab2)
pander(tab2)
pander(fisher.test(tab2))



## -----------------------------------------------------------------------------

t.assay.base = x.d$Baseline.T.cell.Assay.Result..mIU.mL
t.assay.4w = x.d$X4.week.T.cell.Assay.Result..mIU.mL

cor.test(t.assay.base, t.assay.4w, method="kendall")

#plot(t.assay.base, t.assay.4w)

f = 5
plot(jitter(t.assay.base, factor=f), jitter(t.assay.4w, factor=f))
cbind(mrn=x.d$MRN,  t.assay.base,  t.assay.4w)



## -----------------------------------------------------------------------------

spike.res.4w = x.d$X4.week.Spike.Antibody.Result..Positive.or.Negative
tcell.4w = x.d$X4.week.T.cell.Assay.Result..pos..bor..neg..or.inv

tab1 = table(spike.res.4w, tcell.4w)
print(tab1)
pander(tab1)
pander(fisher.test(tab1))

#restrict to just the pos/neg values
tcell.4w[tcell.4w %nin% c("pos", "neg")] = NA
tab2 = table(spike.res.4w, tcell.4w)
print(tab2)
pander(tab2)
pander(fisher.test(tab2))





## -----------------------------------------------------------------------------
y = read.csv("wic_master_20211102.csv")

#clean up the header names
h1 = read.csv("wic_master_20211102.csv", header=F, nrow=1)
#h1 = make.names(gsub(h1[1,], pat="\\(.*\\)", rep=""))
h1 = make.names(gsub(h1, pat="\\.\\.", rep="."))
colnames(y) = h1

#trim whitespace and set to lowercase from every character
for (v in colnames(y)){
    if (is.character(y[,v])){
        y[,v] = tolower(trimws(y[,v]))
        #get rid of wierd white spaces
        y[,v] = trimws(y[, v], whitespace = "[\\h\\v]")
    }
}

#exclude patients from first column
table(y[,"Patients.to.exclude"] == 1)
#y = y[y[,"Patients.to.exclude"] == "",]
y = y[y[,"Patients.to.exclude"] != "1",]

#exclude any row without an MRN
y = y[!is.na(y$MRN),]

#fix the spike.ab numbers
fix.num.vars = c(
    "Spike.Ab.AU.mL",
    "Spike.Ab.AU.mL..4.6.mo.f.u."
)
y[,fix.num.vars]

for (v in fix.num.vars){
    y[,v] = gsub(y[,v], pat="[<,>]", rep="")
    y[,v] = gsub(y[,v], pat="au/ml", rep="")
    y[,v] = as.numeric(y[,v])
}


#fix the test result vars
fix.res.vars = c(
    "Spike.Antibody.Result",
    "Spike.Antibody.Result..4.6.mo.f.u."
)
for (v in fix.res.vars){
    y[y[,v] %nin% c("positive", "negative"),v] = NA
}





## ---- results="asis"----------------------------------------------------------
#run a quick summary for each variable
for (v in colnames(y)){
    cat("##", v, "\n")
    cat("\n")
    if (length(unique(y[,v])) < 40){
        cat(pander(table(y[,v], useNA="ifany")))
        cat("\n")
    }else{
        cat(pander(summary(y[,v])))
        cat("\n")
    }
}

vars.EF = c(
    "Type.of.vaccine",
    "Solid.liquid",
    "Immunosuppressive"
)

y[is.na(y[,"Spike.Antibody.Result..4.6.mo.f.u."]),]




## -----------------------------------------------------------------------------


tab1 = table(y[,c("Spike.Antibody.Result", "Spike.Antibody.Result..4.6.mo.f.u.")])
print(tab1)
#pander(tab1)

ctab1 = ctable(y[,c("Spike.Antibody.Result", "Spike.Antibody.Result..4.6.mo.f.u.")])
pander(ctab1[[1]])

pander(mcnemar.test(tab1))
   
  


## ---- results="asis"----------------------------------------------------------

for (v in vars.EF){
    cat(paste0("###", v, "\n"))
    runStratifiedPairedTables(y,"Spike.Antibody.Result", "Spike.Antibody.Result..4.6.mo.f.u.", v)
    cat("\n\n")
}



## -----------------------------------------------------------------------------

t1 = y[,"Spike.Ab.AU.mL"]
t2 = y[,"Spike.Ab.AU.mL..4.6.mo.f.u."]


#is there a difference betwee the two titer counts?
diff = t2 - t1
wilcox.test(diff)
t.test(t1, t2, paired=T)
#baseline
summary(t1)
#4-6mo followup
summary(t2)

df1 = data.frame(rbind(
    cbind(time=0, samp = 1:length(t1), titer=t1),
    cbind(time=4, samp = 1:length(t1), titer=t2)
))
df1$time = factor(df1$time)
ggplot(df1, aes(x=time, y=titer, color = time, group=samp)) + geom_jitter(width=.05) + geom_line(alpha=.15, color="black") + theme_classic()





## ---- results="asis", warning=F, cache=F--------------------------------------


for (v in vars.EF){

    cat(paste0("###", v, "\n"))


    t1 = y[,"Spike.Ab.AU.mL"]
    t2 = y[,"Spike.Ab.AU.mL..4.6.mo.f.u."]
    group = as.factor(y[,v])
    plotTitersByVar(t1, t2, group, timepoint="4-6mo")
    cat("\n\n")

    #plotTitersByVar(log10(t1), log10(t2), group, timepoint="4-6mo", ylab="log10(titers)")
    #cat("\n\n")


    cat("####Spike.Ab.AU.mL\n") 
    cat(pander(by(t1, group, MeanCI)))
    cat("####Spike.Ab.AU.mL..4.6.mo.f.u.\n") 
    cat(pander(by(t2, group, MeanCI)))


    runStratifiedPairedVals(y,"Spike.Ab.AU.mL", "Spike.Ab.AU.mL..4.6.mo.f.u.", v)
    cat("\n\n")

}


## -----------------------------------------------------------------------------

t1 = y[,"Spike.Ab.AU.mL"]
t2 = y[,"Spike.Ab.AU.mL..4.6.mo.f.u."]

y$rel.base.titer = 1
y$rel.4mo.titer = t2/t1 

wilcox.test(y$rel.4mo.titer)
t.test(t2-t1)



r1 = y$rel.base.titer
r2 = y$rel.4mo.titer 
df1 = data.frame(rbind(
    cbind(time=0, samp = 1:length(t1), titer=r1),
    cbind(time="4-6m", samp = 1:length(t1), titer=r2)
))
df1$time = factor(df1$time)
df1$titer = as.numeric(df1$titer)
ggplot(df1, aes(x=time, y=titer, color = time, group=samp)) + geom_jitter(width=.05) + geom_line(alpha=.15, color="black") + ylim(0, 7)





## -----------------------------------------------------------------------------




## ---- results="asis", warning=F, cache=F--------------------------------------


for (v in vars.EF){

    cat(paste0("###", v, "\n"))


    t1 = y[,"rel.base.titer"]
    t2 = y[,"rel.4mo.titer"]
    group = as.factor(y[,v])
    plotTitersByVar(t1, t2, group, ylim=c(0,7), timepoint="4-6mo")
    cat("\n\n")

    runTiterTables(y[,"Spike.Ab.AU.mL"], y[,"Spike.Ab.AU.mL..4.6.mo.f.u."], group)
    cat("\n\n")


    #runStratifiedPairedVals(y,"Spike.Ab.AU.mL", "Spike.Ab.AU.mL..4.6.mo.f.u.", v)
    #cat("\n\n")

}


## -----------------------------------------------------------------------------

#start with the patients that have a valid post vacc titer
x.3 = x[!is.na(x[,"Post.vacination.Spike.Ab.Result..AU.mL"]),]

#exclude the patients from the 3 titer timepoint analysis
x.3 = x.3[which(x.3$Patients.to.include.in.3.timepoint.analysis==1),]
dim(x.3)
table(x.3$Malig)

titer.vars = c(
    "Post.vacination.Spike.Ab.Result..AU.mL",
    "Baseline.Study.Spike.Ab.Result..AU.mL",
    "X4.week.Spike.Ab.Result..AU.mL"
)

plot3TitersByVar <- function(t1, t2, t3, group, ylim=NA, tps = c(0,16,20), tp.labs = c("t0", "t1", "t2"), ylab="titers"){
    df1 = data.frame(rbind(
        data.frame(time=tps[1], samp = 1:length(t1), titer=t1, group=group),
        data.frame(time=tps[2], samp = 1:length(t1), titer=t2, group=group),
        data.frame(time=tps[3], samp = 1:length(t1), titer=t3, group=group)
    ))
    #df1$time = factor(df1$time, levels=tps) 

    p1 = ggplot(df1, aes(x=time, y=titer, group=samp, color=group)) + geom_jitter(width=.05, size=.75) + geom_line(alpha=.25) + ylab(ylab) + theme_classic()
    if (!is.na(ylim)){
        p1 = p1 + ylim(ylim)
    }

    #calculate CIs
    #https://www.r-graph-gallery.com/4-barplot-with-error-bar.html
    ##add mean trends for each level
    df2 = lapply(unique(group), function(g){
        alpha=0.05
        t=qt((1-alpha)/2 + .5, sum(group==g)-1)   # tend to 1.96 if sample size is big enough
        rbind(
            data.frame(time=tps[1], samp = g, titer=mean(t1[group==g], na.rm=T), titer.ci = sd(t1[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g),
            data.frame(time=tps[2], samp = g, titer=mean(t2[group==g], na.rm=T), titer.ci = sd(t2[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g),
            data.frame(time=tps[3], samp = g, titer=mean(t3[group==g], na.rm=T), titer.ci = sd(t3[group==g], na.rm=T)/sqrt(sum(group==g))*t, group=g)
        )
    })
    df2 = do.call(rbind, df2)
    #df2$time = factor(df2$time, levels=tps) 

    p2 = p1 + geom_line(alpha=1, data=df2, size=1.5)+
        scale_x_continuous(breaks=tps, labels=tp.labs) 
    p3 = p2 + geom_errorbar(aes(ymin = titer - titer.ci, ymax = titer + titer.ci), data=df2, size=.75, width=1.15) +
        scale_x_continuous(breaks=tps, labels=tp.labs) 
    print(p2)
    print(p3)

    df3 = data.frame(rbind(
        data.frame(time=tps[1], titer=t1, group=group),
        data.frame(time=tps[2], titer=t2, group=group),
        data.frame(time=tps[3], titer=t3, group=group)
    ))
    df3$time = factor(df3$time)

    #b1 = ggplot(df3, aes(x=time, y=titer, fill=group)) +  geom_boxplot(outlier.shape=NA) +  geom_point(position=position_jitterdodge(jitter.width=.1), color="black") + ylab(ylab) + theme_classic()
    #    #stat_summary(position=position_dodge(), fun.y=mean, geom="point", size=2) +
    #    #stat_summary(position=position_dodge(), fun.data = mean_se, geom = "errorbar")
    #print(b1)

    b2 = ggplot(df3, aes(x=time, y=titer, fill=group)) +
        stat_summary(position=position_dodge(), fun.data ="mean_cl_boot", geom="crossbar")  + 
        geom_point(position=position_jitterdodge(jitter.width=.1), color="black") + 
        ylab(ylab) + 
        theme_classic()
    print(b2)


}


t1 = x.3[,titer.vars[1]]
t2 = x.3[,titer.vars[2]]
t3 = x.3[,titer.vars[3]]
group = x.3$Malignancy.category..Solid.or.Liquid

tps = c("Post Full Vx", "Pre Booster", "Post Booster")
tps = factor(tps, levels=tps)
plot3TitersByVar(t1, t2, t3, group, tp.labs = c("Post Full Vx", "Pre Booster", "Post Booster"))



## ---- results="asis"----------------------------------------------------------


cat(paste0("####", titer.vars[1]))
cat(pander(by(t1, group, MeanCI)))
cat(paste0("####", titer.vars[2]))
cat(pander(by(t2, group, MeanCI)))
cat(paste0("####", titer.vars[3]))
cat(pander(by(t3, group, MeanCI)))




## -----------------------------------------------------------------------------

knit_exit()

