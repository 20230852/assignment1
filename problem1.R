library(ggplot2)

# Problem 1
## (a) 
dir.create("data",  showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

set.seed(1)
n = 200
x = seq(0, 1, length.out = n)
y = sin(2*pi*x) + rnorm(n, sd = 0.15)

dat = data.frame(x = x, y = y)
write.csv(dat, "data/sim_data.csv", row.names = FALSE)
head(dat)


## (b) 
dat = read.csv("data/sim_data.csv")

p = ggplot(dat, aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, colour = "red") +
  labs(title = "Scatter plot and fitted regression curve (LOESS)",
       x = "x", y = "y") +
  theme_minimal()
p


## (c)
ggsave("plots/regression_fit.png", plot = p, width = 7, height = 5, dpi = 300)
file.exists("plots/regression_fit.png")


## (d)
# https://github.com/20230852/assignment1
