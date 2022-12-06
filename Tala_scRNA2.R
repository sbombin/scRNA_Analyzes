obj01[["percent.mt"]] <- PercentageFeatureSet(obj01, pattern = "^MT-")
obj03[["percent.mt"]] <- PercentageFeatureSet(obj03, pattern = "^MT-")


score.count <- mean(obj01$nCount_RNA) + 3*sd(obj01$nCount_RNA)
score.feature <- mean(obj01$nFeature_RNA) + 3*sd(obj01$nFeature_RNA)



filtr1  <- subset(obj01, subset = nFeature_RNA > 200 & nFeature_RNA < score.feature & percent.mt < 15 &
                    nCount_RNA >500 & nCount_RNA < score.count)

head(filtr1)

score.count3 <- mean(obj03$nCount_RNA) + 3*sd(obj03$nCount_RNA)
score.feature3 <- mean(obj03$nFeature_RNA) + 3*sd(obj03$nFeature_RNA)

filtr3  <- subset(obj03, subset = nFeature_RNA > 200 & nFeature_RNA < score.feature3 & percent.mt < 10 &
                    nCount_RNA >500 & nCount_RNA < score.count3)

object  <- subset(object, subset = orig.ident == "H23" | orig.ident == "CK10341", invert = TRUE)

object  <- subset(object, subset = orig.ident == "H23", invert = TRUE)