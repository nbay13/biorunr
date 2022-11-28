library(tidyr)
library(dplyr)

# initialize a random data.frame
df <- data.frame(
	var1 = c(rnorm(mean = 2, sd = 1, n = 10), var2 = rnorm(mean = 0, sd = 1, n = 20)),
	group = sample(c(1, 2), 30, replace = T)
)

# function for calculating median and IQR values
get.boxplot.data <- function(df, group_var, low = 0.25, high = 0.75){
	out <- data.frame(
		df %>% 
		gather(variable, value, -group_var) %>% 
		group_by(get(group_var), variable) %>% 
		summarise(mid = median(value), lq = quantile(value, low), uq = quantile(value, high))
	)
	colnames(out)[1] <- group_var
	return(out)
}

# get boxplot stats for each group and variable
bp_df <- get.boxplot.data(df, group_var = "group")

# same as melt() from reshape package
plot_df <- df %>% gather(variable, value, - group)

# this value seems to match boxplot width on my PC
x_width = 0.375
# example plot
# Note: Group variable must be converted to numeric for width calculations
a <- ggplot(plot_df, aes(x = factor(group), y = value, fill = factor(group))) + 
geom_blank() +
geom_segment(data = bp_df, aes(x = group - x_width, xend = group + x_width, y = mid, yend = mid),
 color = c("blue", "red")) +
geom_segment(data = bp_df, aes(x = group - x_width, xend = group + x_width, y = uq, yend = uq),
 color = "green") +
geom_segment(data = bp_df, aes(x = group - x_width, xend = group + x_width, y = lq, yend = lq),
 color = c("purple", "orange")) +
geom_point() +
theme_classic() +
theme(legend.position = "none") +
labs(x = "", y = "value")

# add a boxplot layer underneath to see if IQRs look right
b <- ggplot(plot_df, aes(x = factor(group), y = value, fill = factor(group))) + 
geom_boxplot(color = NA) +
geom_segment(data = bp_df, aes(x = group - x_width, xend = group + x_width, y = mid, yend = mid),
 color = c("blue", "red")) +
geom_segment(data = bp_df, aes(x = group - x_width, xend = group + x_width, y = uq, yend = uq),
 color = "green") +
geom_segment(data = bp_df, aes(x = group - x_width, xend = group + x_width, y = lq, yend = lq),
 color = c("purple", "orange")) +
geom_point() +
theme_classic() +
theme(legend.position = "none") +
labs(x = "", y = "value")

library(ggpubr)
ggarrange(a, b, nrow = 2)