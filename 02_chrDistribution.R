chr = read.table(file = "~/Documents/study/biology/courses/bioinformatics/final-project/data/countChr.txt",
                 header = FALSE,
                 stringsAsFactors = FALSE)

autosome = paste0(rep("Chr", 22), seq(1:22))
chrName = append(autosome, c("ChrX", "ChrY"))
row.names(chr) = chrName

mycol = rainbow(24)
pie(chr$V1, labels = row.names(chr), radius = 1, col = mycol, cex = 0.7)
