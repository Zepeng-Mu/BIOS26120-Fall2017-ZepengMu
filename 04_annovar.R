library(ggplot2)
library(dplyr)
anno220 = read.table(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/220annotation/SRR765980_220.annovar.hg38_multianno.txt",
                     header = T,
                     stringsAsFactors = F,
                     sep = "\t")

annoHq = read.csv(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/30-220Annotation/SRR765980_30_220.annovar.hg38_multianno.txt",
                   header = T,
                   stringsAsFactors = F,
                  sep = "\t")


# count and plot percentage for mutation type
func = table(anno220$Func.refGene)
func1 = table(annoHq$Func.refGene)
combine = c(func1, func)
df = data.frame(count = combine, set = c(rep("HQ30 (30-220)", length(func1)), rep("HQ220 (>220)", length(func))))
df$type = c(row.names(func1), row.names(func))
ggplot(data = df, aes(fill = type, y = count, x = set)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "percentage")
ggplot(data = df, aes(fill = type, y = count, x = set)) +
  geom_bar(stat = "identity") +
  labs(y = "count")

# count and plot percentage for exonic mutation type
exo = table(anno220$ExonicFunc.refGene)
exo1 = table(annoHq$ExonicFunc.refGene)
combine1 = c(exo1, exo)
df1 = data.frame(count = combine1, set = c(rep("HQ30 (30-220)", length(exo1)), rep("HQ220 (>220)", length(exo))))
df1$type = c(row.names(exo1), row.names(exo))
df1 = df1[-c(1, 10, 11, 20), ]
df1$count = log10(df1$count)

ggplot(data = df1, aes(fill = type, y = count, x = set)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Log10(count)")

# retrieve OMIM annotated variants
sortOMIM = anno220[grep(anno220$CLNDSDB, pattern = "OMIM"), ]
sortOMIM = sortOMIM[which(sortOMIM$Polyphen2_HVAR_rankscore != "" & sortOMIM$Polyphen2_HVAR_rankscore != "."), ]
sortOMIM = sortOMIM[order(sortOMIM$Polyphen2_HVAR_rankscore, decreasing = T), ]
write.csv(sortOMIM, file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/HQ20.csv",
          quote = F,
          row.names = F)
sortOMIMUpd = sortOMIM[sortOMIM$Polyphen2_HVAR_rankscore >= 0.6 & sortOMIM$SIFT_converted_rankscore >= 0.7, ]

#retrieve GTEx reported eQTLs and get p-values
sortGTEx = anno220[anno220$GTEx_V6_gene != "" & anno220$GTEx_V6_gene != ".", ]
sortGTEx = sortGTEx[sortGTEx$Func.refGene != "exonic", ]
out = matrix(ncol = 3)
for(i in 1:nrow(sortGTEx)) {
  rs = sortGTEx$avsnp150[i]
  ENSID = strsplit(sortGTEx$GTEx_V6_gene[i], split = "|", fixed = T)[[1]]
  tissue = strsplit(sortGTEx$GTEx_V6_tissue[i], split = "|", fixed = T)[[1]]
  combine = cbind(rep(rs, length(ENSID)), ENSID, tissue)
  out = rbind(out, combine)
}

out = out[-1, ]
write.csv(out, file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/eQTLs.csv",
          row.names = F,
          quote = F)

tested = read.csv(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/GTEx_Portal_eQTLs.csv",
                  header = T,
                  stringsAsFactors = F)
tested = tested[order(tested$P.Value), ]

signifTested = filter(tested, P.Value <= 0.05 / nrow(tested))
length(unique(signifTested$SNP))
length(unique(signifTested$Gencode.Id))
length(unique(signifTested$Tissue))

database = read.csv(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/SNP_Catalog.csv",
                      header = T,
                      stringsAsFactors = F)

dbSelected = database[, c("STUDY", "DISEASE.TRAIT", "MAPPED_GENE", "STRONGEST.SNP.RISK.ALLELE", "SNPS", "MAPPED_TRAIT")]
remove(database)

database[which(dbSelected$SNPS %in% tested$SNP), ]
dt = tested[which(tested$SNP %in% dbSelected$SNPS), ]
