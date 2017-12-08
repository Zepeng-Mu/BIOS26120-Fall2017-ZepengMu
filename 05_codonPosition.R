library(seqinr)
library(ggplot2)

anno = read.table(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/annotation/SRR765980.annovar.hg38_multianno.txt",
                  header = T,
                  stringsAsFactors = F,
                  sep = "\t")

View(anno)
nonsynonymous = anno[anno$ExonicFunc.refGene == "nonsynonymous SNV", ]
nonsynonymousPos = integer()
for(line in nonsynonymous$AAChange.refGene) {
  nucleoMatch = regmatches(line, regexpr(pattern = "c.[ATGC][0-9]+[ATGC]", line))
  nucleoMatch = s2c(nucleoMatch)
  nucleoPos = as.integer(c2s(nucleoMatch[4:(length(nucleoMatch) - 1)]))
  
  aaMatch = regmatches(line, regexpr(pattern = "p.[QWERTYIPASDFGHKLCVNM][0-9]+[QWERTYIPASDFGHKLCVNM]", line))
  aaMatch = s2c(aaMatch)
  aaPos = as.integer(c2s(aaMatch[4:(length(aaMatch) - 1)]))
  
  codonPos = nucleoPos - (aaPos - 1) * 3
  
  nonsynonymousPos = append(nonsynonymousPos, codonPos)
}

synonymous = anno[anno$ExonicFunc.refGene == "synonymous SNV", ]
synonymousPos = integer()
for(line in synonymous$AAChange.refGene) {
  nucleoMatch = regmatches(line, regexpr(pattern = "c.[ATGC][0-9]+[ATGC]", line))
  nucleoMatch = s2c(nucleoMatch)
  nucleoPos = as.integer(c2s(nucleoMatch[4:(length(nucleoMatch) - 1)]))
  
  aaMatch = regmatches(line, regexpr(pattern = "p.[QWERTYIPASDFGHKLCVNMX][0-9]+[QWERTYIPASDFGHKLCVNMX]", line))
  aaMatch = s2c(aaMatch)
  aaPos = as.integer(c2s(aaMatch[4:(length(aaMatch) - 1)]))
  
  codonPos = nucleoPos - (aaPos - 1) * 3
  
  synonymousPos = append(synonymousPos, codonPos)
}

table(synonymousPos)
table(nonsynonymousPos)

summ = data.frame(
  Position = rep(c("1", "2", "3"), 2),
  Type = rep(c("Nonsynonymous", "Synonymous"), each = 3),
  Count = c(800, 756, 119, 69, 6, 1864)
)

ggplot(data = summ, aes(fill = Position, y = Count, x = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "percentage")
