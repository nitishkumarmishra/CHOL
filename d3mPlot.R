## ppt function for beta distribution plot
ppt <- function (cases, control) 
{
    graphics::par(mar = c(2, 2, 2, 1))
    graphics::par(mfrow = c(1, 2))
    cases <- as.numeric(cases)
    control <- as.numeric(control)
    beanplot::beanplot(cases, control, side = "both", col = list(c("red","black"), c("blue", "black")), ylim = c(-0.1, 1), ll = 0.05, names = c("cases", "control"), horizontal = T, log = "", main = "Density of beta value")
    graphics::grid()
    stats::qqplot(cases, control, xlim = c(0, 1), ylim = c(0, 1), main = "Q-Q plot")
    graphics::abline(a = 0, b = 1, col = "red")
}
#### Beta distribution plot
ppt(cases, control)
