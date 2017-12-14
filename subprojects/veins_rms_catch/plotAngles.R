library(ggplot2)
library(viridis)

d <- read.csv("angles.csv")
p <- ggplot(d, aes(x=x, y=phi*180/pi)) + geom_line()
print(p)
