ALL_DST <- c("PGW", "LN_", "LL_")
ALL_CFT <- c("ABCD",
             "ABST",
             "ABXT",
             "ABSS", "ABTT",
             "XBXT",
             "AASS", "AATT", "BBSS", "BBTT",
             "ABXX",
             "XBXX",
             "AAXX", "BBXX",
             "AXXX")

ALL_MDL <- do.call(what = paste, args = c(expand.grid(ALL_DST, ALL_CFT), sep = ""))
