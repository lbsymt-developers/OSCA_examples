library(BiocFileCache)
path <- tempfile()
bfc <- BiocFileCache(path, ask = FALSE)
bfc <- BiocFileCache(ask=FALSE)
url <- file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
                 "GSE85nnn/GSE85241/suppl",
                 "GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz")

grubman.fname <- bfcrpath(bfc, url)
local.name <- URLdecode(basename(url))
unlink(local.name)
