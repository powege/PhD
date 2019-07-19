rm(list=ls())
graphics.off()

### Script that identifies all POS with < 10X coverage in  > 10% of strains 


dt <- data.table(CHR = c(1,1,1,2,2,2,3,3,3),
                 POS = c(1,1,6,1,1,3,1,1,6))
x <- table(dt)
x <- as.data.table(x)
y <- subset(x, x$N >= 2)
y <- y[,1:2]

