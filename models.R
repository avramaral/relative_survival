# See this: https://stackoverflow.com/questions/66877063/how-to-define-global-variables-across-multiple-functions-in-r-not-in-a-package

ALL_DST <- c("PGW", "LN_", "LL_", "GG_", "GAM")
ALL_CFT <- c("ABCD",
             "ABXD",
             "ABCC", "ABDD",
             "XBXD",
             "AACC", "AADD", "BBCC", "BBDD",
             "ABYZ",
             "ABXZ", 
             "ABYY", "ABZZ",
             "XBXZ",
             "AAYY", "AAZZ", "BBYY", "BBZZ",
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
