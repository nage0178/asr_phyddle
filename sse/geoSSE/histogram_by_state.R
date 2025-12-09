method_diff <- read.csv("~/asr_phyddle/sse/geoSSE/bayes/difference_phy-bayes.csv")
#colnames(method_diff)

library(tidyr)

long <- method_diff %>%
  pivot_longer(
    cols = `X1`:`X8`,
    names_to = "state",
    values_to = "diff"
  )
long

head(long)
df <- as.data.frame(long)


state_names <- c(
  X1="A->A,A",
  X2="B->B,B",
  X3="AB->A,B",
  X4="AB->B,A",
  X5="AB->AB,A",
  X6="AB->A,AB",
  X7="AB->AB,B",
  X8="AB->B,AB"
)

state_labeller <- function(variable,value){
  return(state_names[value])
}

df <- df %>% group_by(state) %>%  mutate(mean = mean(diff))


plot1 <- ggplot(df, aes(x = diff)) + geom_histogram(aes(y = after_stat(count / sum(count)))) + 
  facet_wrap(~state, , labeller=labeller(state = state_names)) +
  theme_bw() + labs(x = "phyddle probability -\nBayesian inference probability", y = "proportion of inferences") + 
  geom_vline(aes(xintercept = mean, group = state), colour = 'red')


pdf("~/asr_writing/manuscript/figs/geosse_hist.pdf", width = 8, height = 5)
plot1
dev.off()

png("~/asr_writing/manuscript/figs/geosse_hist.png", width = 8, height = 5, units= "in", res = 500)
plot1
dev.off()
