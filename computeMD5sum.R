#### create MD5sum files

library(tools)

## list of files in package
fileList = system(paste0("git ls-files pkg"), intern=TRUE)

## md5sum creation
md5sumList = md5sum(fileList)
str(md5sumList)

## output in data.frame
md5sumTable = data.frame(MD5=md5sumList, files=names(md5sumList), stringsAsFactors=FALSE)
rownames(md5sumTable)= NULL

## file name at the good format
md5sumTable$files = sapply(md5sumTable$files, function(word) sub(pattern="pkg/", replacement="*", x=word))

## remove the MD5sum for the MD5 file
index.MD5 = which(md5sumTable$files == "*MD5")
md5sumTable = md5sumTable[-index.MD5,]

## save the table in the MD5 file
write.table(x=md5sumTable, file="./pkg/MD5", quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
