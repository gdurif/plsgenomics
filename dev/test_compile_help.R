### generate LaTex file from .Rd files to check for errors

library(tools)

## list of files in package
fileList = system(paste0("git ls-files pkg/man"), intern=TRUE)

# generate LaTeX file
for(file in fileList) {
     Rd2latex(file, out = "/dev/null", defines = .Platform$OS.type,
         stages = "render", outputEncoding = "ASCII",
         fragment = FALSE, writeEncoding = TRUE)
}


 
# Avis : pkg/man/mrpls.Rd: Some input could not be re-encoded to ASCII
# Avis : pkg/man/mrpls.cv.Rd: Some input could not be re-encoded to ASCII
# Avis : pkg/man/spls.adapt.Rd: Some input could not be re-encoded to ASCII
# Avis : pkg/man/spls.adapt.tune.Rd: Some input could not be re-encoded to ASCII



Rd2latex("pkg/man/mrpls.Rd", out = "test.tex", defines = .Platform$OS.type,
         stages = "render", outputEncoding = "ASCII",
         fragment = FALSE, writeEncoding = TRUE)

Rd2latex("pkg/man/mrpls.cv.Rd", out = "test.tex", defines = .Platform$OS.type,
         stages = "render", outputEncoding = "ASCII",
         fragment = FALSE, writeEncoding = TRUE)

Rd2latex("pkg/man/spls.adapt.Rd", out = "test.tex", defines = .Platform$OS.type,
         stages = "render", outputEncoding = "ASCII",
         fragment = FALSE, writeEncoding = TRUE)

Rd2latex("pkg/man/spls.adapt.tune.Rd", out = "test.tex", defines = .Platform$OS.type,
         stages = "render", outputEncoding = "ASCII",
         fragment = FALSE, writeEncoding = TRUE)
