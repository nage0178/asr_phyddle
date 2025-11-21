# truth <- read.csv("~/asr_phyddle/geoSSE/2_states/parent_l_r/larger/random_daughters/larger/estimate/sim.test_true.labels_cat.csv")
# est <- read.csv("~/asr_phyddle/geoSSE/2_states/parent_l_r/larger/random_daughters/estimate/sim.test_est.labels_cat.csv")

est <- read.csv("~/asr_phyddle/sse/geoSSE/estimate/sim.test_est.labels_cat.csv")
truth <- read.csv("~/asr_phyddle/sse/geoSSE/estimate/sim.test_true.labels_cat.csv")
est <- est[1:2500, ]
truth <- truth[1:2500,]

truth <- truth[, -1]
est <- est[, -1]
nstates <- 8
n_node <- 49
find_node <- function(row){
  which(max(row) == row)
}

point_est <- matrix(NA, nrow = dim(truth)[1], ncol= n_node)
for (i in 1: n_node) {
  col <- ((i-1) * 8 + 1):(i *8)
  point_est[, i] <- apply(est[, col], 1, find_node)
}

point_est_v <- as.factor(c(point_est))
truth_v <- as.factor(c(as.matrix(truth + 1)))

library(caret)
cm <- confusionMatrix(point_est_v, truth_v)
mean(point_est_v == truth_v)
# extract the confusion matrix values as data.frame
cm_d <- as.data.frame(cm$table)
cm_d <- cbind(cm_d, NA)
colnames(cm_d)[4] <- "proportion"

for (state in 1:nstates) {
  rows <- which(cm_d$Reference == state)
  cm_d$proportion[rows] <- cm_d$Freq[rows] / sum(cm_d$Freq[rows])
}


# # confusion matrix statistics as data.frame
# cm_st <-data.frame(cm$overall)
# # round the values
# cm_st$cm.overall <- round(cm_st$cm.overall,2)
# 
# # here we also have the rounded percentage values
# cm_p <- as.data.frame(prop.table(cm$table))
# cm_d$Perc <- round(cm_p$Freq*100,2)


library(ggplot2)     # to plot
library(gridExtra)   # to put more
library(grid)        # plot together

which(cm_d$Reference== 1)
cm_d$true_state[which(cm_d$Reference== 1)] <- "A->A,A"
cm_d$true_state[which(cm_d$Reference== 2)] <- "B->B,B"
cm_d$true_state[which(cm_d$Reference== 3)] <- "AB->A,B"
cm_d$true_state[which(cm_d$Reference== 4)] <- "AB->B,A"
cm_d$true_state[which(cm_d$Reference== 5)] <- "AB->AB,A"
cm_d$true_state[which(cm_d$Reference== 6)] <- "AB->A,AB"
cm_d$true_state[which(cm_d$Reference== 7)] <- "AB->AB,B"
cm_d$true_state[which(cm_d$Reference== 8)] <- "AB->B,AB"
cm_d$true_state <- factor(cm_d$true_state, levels = c( "A->A,A", "B->B,B", "AB->A,B", "AB->B,A", "AB->AB,A", "AB->A,AB", "AB->AB,B", "AB->B,AB"))

cm_d$inf_state <- NA
cm_d$inf_state[which(cm_d$Prediction== 1)] <- "A->A,A"
cm_d$inf_state[which(cm_d$Prediction== 2)] <- "B->B,B"
cm_d$inf_state[which(cm_d$Prediction== 3)] <- "AB->A,B"
cm_d$inf_state[which(cm_d$Prediction== 4)] <- "AB->B,A"
cm_d$inf_state[which(cm_d$Prediction== 5)] <- "AB->AB,A"
cm_d$inf_state[which(cm_d$Prediction== 6)] <- "AB->A,AB"
cm_d$inf_state[which(cm_d$Prediction== 7)] <- "AB->AB,B"
cm_d$inf_state[which(cm_d$Prediction== 8)] <- "AB->B,AB"
cm_d$inf_state <- factor(cm_d$inf_state, levels = c( "A->A,A", "B->B,B", "AB->A,B", "AB->B,A", "AB->AB,A", "AB->A,AB", "AB->AB,B", "AB->B,AB"))


pdf("~/asr_writing/manuscript/figs/confusion_geosse.pdf", width =4, height = 5)
# plotting the matrix
ggplot(data = cm_d, aes(y = inf_state , x = true_state, fill = proportion))+
  labs(x = "true ancestral states", y = "inferred ancestral states")+
  geom_tile() +
  geom_text(aes(label = round(proportion, digits = 2)), col = "red") +
  
  #geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'red', size = 8) +
  theme_light() +   theme(legend.position = "top") + theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=1), axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1))

dev.off()


bayes_conf <- read.csv("~/asr_phyddle/sse/geoSSE/bayes/bayes_confusion.csv")


b_d <- as.data.frame(cbind(bayes_conf, NA))
colnames(b_d)[4] <- "proportion"

for (state in 1:8) {
  rows <- which(b_d$Reference == state)
  b_d$proportion[rows] <- b_d$Freq[rows] / sum(b_d$Freq[rows])
}

b_d$true_state[which(b_d$Reference== 1)] <- "A->A,A"
b_d$true_state[which(b_d$Reference== 2)] <- "B->B,B"
b_d$true_state[which(b_d$Reference== 3)] <- "AB->A,B"
b_d$true_state[which(b_d$Reference== 4)] <- "AB->B,A"
b_d$true_state[which(b_d$Reference== 5)] <- "AB->AB,A"
b_d$true_state[which(b_d$Reference== 6)] <- "AB->A,AB"
b_d$true_state[which(b_d$Reference== 7)] <- "AB->AB,B"
b_d$true_state[which(b_d$Reference== 8)] <- "AB->B,AB"
b_d$true_state <- factor(b_d$true_state, levels = c( "A->A,A", "B->B,B", "AB->A,B", "AB->B,A", "AB->AB,A", "AB->A,AB", "AB->AB,B", "AB->B,AB"))

b_d$inf_state <- NA
b_d$inf_state[which(b_d$Prediction== 1)] <- "A->A,A"
b_d$inf_state[which(b_d$Prediction== 2)] <- "B->B,B"
b_d$inf_state[which(b_d$Prediction== 3)] <- "AB->A,B"
b_d$inf_state[which(b_d$Prediction== 4)] <- "AB->B,A"
b_d$inf_state[which(b_d$Prediction== 5)] <- "AB->AB,A"
b_d$inf_state[which(b_d$Prediction== 6)] <- "AB->A,AB"
b_d$inf_state[which(b_d$Prediction== 7)] <- "AB->AB,B"
b_d$inf_state[which(b_d$Prediction== 8)] <- "AB->B,AB"
b_d$inf_state <- factor(b_d$inf_state, levels = c( "A->A,A", "B->B,B", "AB->A,B", "AB->B,A", "AB->AB,A", "AB->A,AB", "AB->AB,B", "AB->B,AB"))
pdf("~/asr_writing/manuscript/figs/confusion_geosse_bayes.pdf", width =4, height = 5)
# plotting the matrix
ggplot(data = b_d, aes(y = inf_state , x = true_state, fill = proportion))+
  labs(x = "true ancestral states", y = "inferred ancestral states")+
  geom_tile() +
  geom_text(aes(label = round(proportion, digits = 2)), col = "red") +
  #geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'red', size = 8) +
  theme_light() +   theme(legend.position = "top") + theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=1), axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1))

dev.off()
both_method <- rbind(cm_d, b_d)
both_method <- cbind(both_method, c(rep("phyddle", nrow(cm_d)), rep("Bayesian", nrow(b_d))))
colnames(both_method)[7] <- "method"
pdf("~/asr_writing/manuscript/figs/confusion_geosse_both.pdf", width =10, height = 5)
ggplot(data = both_method, aes(y = inf_state , x = true_state, fill = proportion))+
  facet_wrap(~method)+
  labs(x = "true ancestral states", y = "inferred ancestral states")+
  geom_tile() +
  geom_text(aes(label = round(proportion, digits = 2)), col = "red") +
  #geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'red', size = 8) +
  theme_light()  + theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=1), axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

png("~/asr_writing/manuscript/figs/confusion_geosse_both.png", width =10, height = 5, units = "in", res = 500)
ggplot(data = both_method, aes(y = inf_state , x = true_state, fill = proportion))+
  facet_wrap(~method)+
  labs(x = "true ancestral states", y = "inferred ancestral states")+
  geom_tile() +
  geom_text(aes(label = round(proportion, digits = 2)), col = "red") +
  #geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'red', size = 8) +
  theme_light()  + theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=1), axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()
