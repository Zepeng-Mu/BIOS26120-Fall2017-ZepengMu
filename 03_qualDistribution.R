dt1 = read.table(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/SRR765980.flt.vcf",
                header = FALSE,
                stringsAsFactors = FALSE,
                skip = 3395)

myHist = hist(dt1$V6, breaks = 40, plot = FALSE)
myCol = ifelse(myHist$breaks < 30, "grey", ifelse(myHist$breaks >= 220, "hotpink2", "deepskyblue2"))
plot(myHist, col = myCol, border = FALSE, main = "", xlim = c(0, 250), ylim = c(0, 32000), xlab = "Quality", ylab = "Count")

myHist = hist(dt1$V6, breaks = 40, plot = FALSE)
myHist$counts = log10(myHist$counts)
myCol = ifelse(myHist$breaks < 30, "grey", ifelse(myHist$breaks >= 220, "hotpink2", "deepskyblue2"))
plot(myHist, col = myCol, border = FALSE, main = "", xlim = c(0, 250), xlab = "Quality", ylab = "Log10(Count)")
