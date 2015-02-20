#### create MD5sum files

library(tools)

fileList = system(paste0("git ls-files pkg"), intern=TRUE)

md5sumList = md5sum(fileList)
str(md5sumList)

md5sumTable = data.frame(MD5=md5sumList, files=names(md5sumList), stringsAsFactors=FALSE)
rownames(md5sumTable)= NULL

md5sumTable$Files = sapply(md5sumTable$Files, function(word) sub(pattern="pkg/", replacement="*", x=word, perl=TRUE, fixed=TRUE))
