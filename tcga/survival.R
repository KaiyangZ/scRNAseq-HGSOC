#!/usr/bin/env Rscript

# This script is for scoring stress-related transcriptional profile in TCGA samples, and check the association between the stress score and progression-free survival (PFS), desease-free interval (DFI)

## the deconvolutation of TCGA data was performed as described in Häkkinen et al., 2019 (https://doi.org/10.1101/854505)

# read bulk and deconvoluted TCGA-OV expression 

# set path to deconvoluted dir
decom.dir <- 'rnaseq-decom-tcga/ov/2020_04_01/'

# read bulk
read_expr <- function(file, ...) {
    tab <- read.table(file, ..., row.names = NULL)
    rownames(tab) <- tab[, 1]   # NB. put gene symbols as names, drop entrez IDs
    tab <- as.matrix(tab[, -(1:2)])
    return (tab)
}

B <- read_expr(file.path("rnaseq-decom-tcga/ov/2018_11_23/", 'tcga_ov-expr-vaharautio.tsv.gz'), sep = '\t', header = T)

# read decomposed 
G <- as.matrix(read.table(file.path(decom.dir, 'out/G.tsv.gz'),
    sep = '\t', header = T, row.names = 1L))
W <- as.matrix(read.table(file.path(decom.dir, 'out/W.tsv.gz'),
    sep = '\t', header = T, row.names = 1L))
Z <- read_expr(file.path(decom.dir, 'out/Z.tsv.gz'),
    sep = '\t', header = T)

#' read TCGA clinical data and filter samples based on stage and grade
clinical.pan <- read.table("TCGA_clinical_data_panCancer.csv", sep = "\t", header = T, row.names = 1, stringsAsFactors = F, comment.char = "")

# get ov data
clinical.ov.flt <- subset(clinical.pan, type == "OV" & clinical_stage %in% c("Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV") & histological_grade %in% c("G2", "G3", "G4"))

clinical.ov.flt$bcr_patient_barcode <- gsub("-", "\\.", clinical.ov.flt$bcr_patient_barcode)

# filter out recurrent tumors from the RNAseq data
B <- B[, substr(colnames(B), 1, 12) %in% clinical.ov.flt$bcr_patient_barcode & substr(colnames(B), 14,15) == "01"]
clinical.ov.flt <- subset(clinical.ov.flt, bcr_patient_barcode %in% substr(colnames(B), 1, 12) )

G = G[, colnames(B), drop = F]
Z = Z[, substr(colnames(Z), 1, 28) %in% colnames(B)]
W = W[gsub("\\.", "-", colnames(B)), ]

unit.ranks <- function(X) {
    for (j in seq_len(ncol(X))) {
        m <- !is.na(X[, j])
        X[m, j] <- ( rank( X[m, j] ) - 1. ) / ( sum(m) - 1. )
    }
    return (X)
}

# normalize & rank transform expression values
Zn <- t( t(Z) / c(t( c(G) * W )) )

# only take EOC profiles 
Zn.EOC = Zn[, grepl("\\.EOC", colnames(Zn))]
Zr.EOC <- t(unit.ranks( t(Zn.EOC[rowSums(Zn.EOC) > 0, ]) ))

# get gene sets
clusterPair.lik.genes.maxCore.c <- readRDS("int/clusterPair.lik.genes.maxCore.c.RDS")

# compute scores
source('ssgsea/ssgsea.R')
get.scores <- function(X, S) {
        ress <- NULL
        for (j in seq_len(length(S))) {
                set.seed(123L)
                signature <- names(S)[j]
                res <- ssgsea( X, S[[j]], alpha = .75, compute.p = 1000L )
                ress <- c( ress, `names<-`( list(res), signature ) )
        }
        return (do.call(cbind, ress))
}
ssgsea.Zr.EOC <- get.scores(Zr.EOC, clusterPair.lik.genes.maxCore.c)
saveRDS(ssgsea.Zr.EOC, file = "int/ssgsea.Zr.EOC.RDS")

ssgsea.EOC.es <- ssgsea.Zr.EOC[, grepl("gc[0-9]+\\.es", colnames(ssgsea.Zr.EOC), ignore.case = T)]
ssgsea.EOC.p <- ssgsea.Zr.EOC[, grepl("gc[0-9]+\\.p", colnames(ssgsea.Zr.EOC), ignore.case = T)]

saveRDS(ssgsea.EOC.es, file = "int/ssgsea.EOC.es.RDS")
saveRDS(ssgsea.EOC.p, file = "int/ssgsea.EOC.p.RDS")

# Kaplan-Meier estimate on PFS/DFI and stress score (gc6)
library(survival)
library(survminer)
library(RTCGA.clinical)

clinical.ov.flt$PFI.month <- as.numeric(clinical.ov.flt$PFI.time)/30
clinical.ov.flt$DFI.month <- as.numeric(clinical.ov.flt$DFI.time)/30

clinical.ov.flt$stress.p <- ssgsea.EOC.p$gc6.p[match(clinical.ov.flt$bcr_patient_barcode, substr(rownames(ssgsea.EOC.p), 1, 12) )]
clinical.ov.flt$stess.es <- ssgsea.EOC.es$gc6.es[match(clinical.ov.flt$bcr_patient_barcode, substr(rownames(ssgsea.EOC.p), 1, 12) )]

clinical.ov.flt$stress.cat <- ifelse(clinical.ov.flt$stress.p < 0.5, "low", ifelse(clinical.ov.flt$stress.p > 0.95, "high", NA))

clinical.ov.flt$stress.cat = factor(clinical.ov.flt$stress.cat, levels = c("high", "low"))

# fit for PFS ~ stress and DFI ~ stress
stress.fit.dfi <- survfit(Surv(DFI.month, as.logical(as.numeric(DFI))) ~ stress.cat, data = clinical.ov.flt )
stress.fit.pfi <- survfit(Surv(PFI.month, as.logical(as.numeric(PFI))) ~ stress.cat, data = clinical.ov.flt )

stress.km.dfi = ggsurvplot(stress.fit.dfi, data = clinical.ov.flt,
                            legend.title = NULL,
                            legend.labs = c("stress-high", "stress-low"),
                            pval = TRUE,
                            xlim = c(0, 60),
                            conf.int = FALSE,
                            # Add risk table
                            risk.table = TRUE,
                            tables.height = 0.2,
                            break.x.by = 10,
                            xlab = "Months",
                            ylab = "Disease-free survival probability",
                            palette = c("red3", "royalblue3"),
                            pval.size = 5,
                            fontsiz = 5,
                            ggtheme = theme_bw() + theme(axis.line = element_line(colour = "black"),text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
)

pdf("KM-stress_DFI.pdf", onefile = F, width = 7, height = 7)
stress.km.dfi
dev.off()

stress.km.pfi = ggsurvplot(stress.fit.pfi, data = clinical.ov.flt,
                            # Change legends: title & labels
                            legend.title = NULL,
                            legend.labs = c("stress-high", "stress-low"),
                            pval = TRUE,
                            xlim = c(0, 60),
                            conf.int = FALSE,
                            # Add risk table
                            risk.table = TRUE,
                            tables.height = 0.2,
                            break.x.by = 10,
                            xlab = "Months",
                            ylab = "Progression-free survival probability",
                            palette = c("red3", "royalblue3"),
                            pval.size = 5,
                            fontsiz = 5,
                            ggtheme = theme_bw() + theme(axis.line = element_line(colour = "black"),text = element_text(size=16),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
                        )

pdf("KM-stress_PFS.pdf", onefile = F, width = 7, height = 7)
stress.km.pfi
dev.off()

# cox-regression, using age and tumor purity as covariates

library(gridExtra)

clinical.ov.flt$dfiObj <- with(clinical.ov.flt, Surv(DFI.month, as.logical(as.numeric(DFI))))
clinical.ov.flt$pfiObj <- with(clinical.ov.flt, Surv(PFI.month, as.logical(as.numeric(PFI))))
clinical.ov.flt$age_at_initial_pathologic_diagnosis <- as.numeric(clinical.ov.flt$age_at_initial_pathologic_diagnosis)

clinical.ov.flt$stressScore <- ssgsea.EOC.es[match(clinical.ov.flt$bcr_patient_barcode, substr(rownames(ssgsea.EOC.es), 1, 12)), "gc6.es"]
clinical.ov.flt$stressScore.z = scale(clinical.ov.flt$stressScore)

# add signature3 status
sig3 <- read.table("TCGA_signature_proportions_OvCa_180316AV.txt", sep = "\t", stringsAsFactors = F, header = T)
clinical.ov.flt$sig3 <- sig3$sig3_pct[match(clinical.ov.flt$bcr_patient_barcode, gsub("-", "\\.", sig3$sample))]
clinical.ov.flt$sig3.cat <- ifelse(clinical.ov.flt$sig3 > 0, "positive", "negative")

saveRDS(clinical.ov.flt, file = "int/clinical.ov.flt.RDS")

# add tumor purity
clinical.ov.flt$purity <- W[match(clinical.ov.flt$bcr_patient_barcode, gsub("-", "\\.", substr(rownames(W), 1, 12))), "EOC"]

cox.stress.dfi <- coxph(as.formula(paste("dfiObj", paste(c("sig3.cat", "age_at_initial_pathologic_diagnosis", "purity", "stressScore.z"), collapse = "+"), sep = " ~ ")), data = clinical.ov.flt)
cox.stress.pfi <- coxph(as.formula(paste("pfiObj", paste(c("sig3.cat", "age_at_initial_pathologic_diagnosis", "purity", "stressScore.z" ), collapse = "+"), sep = " ~ ")), data = clinical.ov.flt)

# make forest plot
pdf("cox_multVars_forest-PFS.pdf", width = 4, height = 4, onefile = F)
ggforest(cox.stress.pfi, data = clinical.ov.flt)
dev.off()

pdf("cox_multVars_forest-DFI.pdf", width = 4, height = 4, onefile = F)
ggforest(cox.stress.dfi, data = clinical.ov.flt)
dev.off()

# fisher exact test to check if any stress group is enriched in either sig3 pos or sig3 neg patients
stress.sig3.fisher.pval <- fisher.test(matrix(c(nrow(subset(clinical.ov.flt, stress.cat == "high" & sig3.cat == "positive")), nrow(subset(clinical.ov.flt, stress.cat == "high" & sig3.cat == "negative") ), nrow(subset(clinical.ov.flt, stress.cat == "low" & sig3.cat == "positive")), nrow(subset(clinical.ov.flt, stress.cat == "low" & sig3.cat == "negative") )), 2, 2))

# stratify signature3 
clinical.ov.flt$sig3.cat = factor(clinical.ov.flt$sig3.cat, levels = c("negative", "positive"))

stress.sig3.fit.pfi <- survfit(Surv(PFI.month, as.logical(as.numeric(PFI))) ~ sig3.cat + stress.cat, data = clinical.ov.flt )
stress.sig3.fit.dfi <- survfit(Surv(DFI.month, as.logical(as.numeric(DFI))) ~ sig3.cat + stress.cat, data = clinical.ov.flt )

stress.sig3.km.PFI = ggsurvplot(stress.sig3.fit.pfi, data = clinical.ov.flt,
                            # Change legends: title & labels
                            legend.title = NULL,
                            legend.labs = c("stress-high & signature3 negative", "stress-low & signature3 negative", "stress-high & sigature3 positive", "stress-low & signature3 positive"),
                            pval = TRUE,
                            xlim = c(0, 60),
                            conf.int = FALSE,
                            risk.table = TRUE,
                            tables.height = 0.2,
                            break.x.by = 10,
                            xlab = "Months",
                            ylab = "Progression-free survival probability",
                            palette = c("red4", "royalblue4", "#FF0000B3", "#4876FFB3"),
                            pval.size = 5,
                            fontsiz = 5,
                            ggtheme = theme_bw() + theme(axis.line = element_line(colour = "black"),text = element_text(size=12),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
                        )

pdf("KM-stress_sig3-PFS.pdf", onefile = F, width = 12, height = 8)
stress.sig3.km.PFI
dev.off()

stress.sig3.km.DFI = ggsurvplot(stress.sig3.fit.dfi, data = clinical.ov.flt,
                            # Change legends: title & labels
                            legend.title = NULL,
                            legend.labs = c("stress-high & signature3 negative", "stress-low & signature3 negative", "stress-high & sigature3 positive", "stress-low & signature3 positive"),
                            pval = TRUE,
                            xlim = c(0, 60),
                            conf.int = FALSE,
                            risk.table = TRUE,
                            tables.height = 0.2,
                            break.x.by = 10,
                            xlab = "Months",
                            ylab = "Disease-free survival probability",
                            palette = c("red4", "royalblue4", "#FF0000B3", "#4876FFB3"),
                            pval.size = 5,
                            fontsiz = 5,
                            ggtheme = theme_bw() + theme(axis.line = element_line(colour = "black"),text = element_text(size=12),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),  axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15))
                        )

pdf("KM-stress_sig3-DFI.pdf", onefile = F, width = 9, height = 8)
stress.sig3.km.DFI
dev.off()
