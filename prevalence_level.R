
sp_descr <- read.delim("species_descr.csv", sep = ",", header = T)

sp_descr$Prevalence.Level <- ""
sp_descr$Prevalence.Numeric <- as.numeric(gsub("%","",as.character(sp_descr$Prevalence)))

q <- quantile(sp_descr$Prevalence.Numeric, probs=c(0.33, 0.66))

print(q)

sp_descr[sp_descr$Prevalence.Numeric<q[1],]$Prevalence.Level <- "l"
sp_descr[sp_descr$Prevalence.Numeric >= q[1] & sp_descr$Prevalence.Numeric < q[2],]$Prevalence.Level <- "m"
sp_descr[sp_descr$Prevalence.Numeric >= q[2],]$Prevalence.Level <- "h"


print(sp_descr)

write.csv(sp_descr,"species_descr_ext.csv", row.names = FALSE)