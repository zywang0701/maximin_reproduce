### objective: hist 2*3 (I1), (I8), (I9)'s oracle hist + 1st weight

## (I1)
## do summary-task7.R
library(ggplot2)
library(ggpubr)
df = data.frame(matrix(mat.store.List.all[[1]][,1], ncol=1))
colnames(df) = c("est")
p1=ggplot(df, aes(x=est))+
  geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(xintercept = trueList.List.all[[1]]$true.val, color="red", size=1)+
  geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
  labs(x="est", title="Non-regularity + Instability")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(0.3, 1.0)
p1

df = data.frame(matrix(weight.sim.List.all[[1]][,4], ncol=1))
colnames(df) = c("weight")
p2=ggplot(df, aes(x=weight))+
  geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+
  # geom_vline(xintercept = trueList.List.all[[1]]$weight.true[,1], color="red", size=1)+
  # geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
  labs(x="4th group's weight")+
  theme(plot.title = element_text(hjust = 0.5))
p2

## (I8)
## do summary-irregular.R in task16-measure
df = data.frame(matrix(mat.store.List.all[[1]][,1], ncol=1))
colnames(df) = c("est")
p3=ggplot(df, aes(x=est))+
  geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(xintercept = trueList.List.all[[1]]$true.val, color="red", size=1)+
  geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
  labs(x="est", title="Non-regularity")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(c(-0.25, 0.45))
p3

df = data.frame(matrix(weight.sim.List.all[[1]][,1], ncol=1))
colnames(df) = c("weight")
p4=ggplot(df, aes(x=weight))+
  geom_histogram(aes(y=..density..), colour="black",fill="white", bins=60)+
  geom_density(alpha=.2, fill="#FF6666")+
  # geom_vline(xintercept = trueList.List.all[[1]]$weight.true[,1], color="red", size=1)+
  # geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
  labs(x="1st group's weight")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(c(-0.005, 1.0))
p4

# ## (I9)
# ## do summary-high.R in task13-subsampling, setting==4
# df = data.frame(matrix(mat.store.List.all[[1]][,1], ncol=1))
# colnames(df) = c("est")
# p3=ggplot(df, aes(x=est))+
#   geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
#   geom_density(alpha=.2, fill="#FF6666")+
#   geom_vline(xintercept = trueList.List.all[[1]]$true.val, color="red", size=1)+
#   geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
#   labs(x="est", title="Irregularity")+
#   theme(plot.title = element_text(hjust = 0.5))
# p3
# 
# df = data.frame(matrix(weight.sim.List.all[[1]][,1], ncol=1))
# colnames(df) = c("weight")
# p4=ggplot(df, aes(x=weight))+
#   geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
#   geom_density(alpha=.2, fill="#FF6666")+
#   # geom_vline(xintercept = trueList.List.all[[1]]$weight.true[,1], color="red", size=1)+
#   # geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
#   labs(x="1st group's weight")+
#   theme(plot.title = element_text(hjust = 0.5))
# p4

## (I10)
## do summary-high.R in task13-subsampling, setting==5
df = data.frame(matrix(mat.store.List.all[[1]][,1], ncol=1))
colnames(df) = c("est")
p5=ggplot(df, aes(x=est))+
  geom_histogram(aes(y=..density..), colour="black",fill="white", bins=50)+
  geom_density(alpha=.2, fill="#FF6666")+
  geom_vline(xintercept = trueList.List.all[[1]]$true.val, color="red", size=1)+
  geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
  labs(x="est", title="Normaltiy")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(c(-0.35, 0.35))
p5

df = data.frame(matrix(weight.sim.List.all[[1]][,1], ncol=1))
colnames(df) = c("weight")
p6=ggplot(df, aes(x=weight))+
  geom_histogram(aes(y=..density..), colour="black",fill="white", bins=60)+
  geom_density(alpha=.2, fill="#FF6666")+
  # geom_vline(xintercept = trueList.List.all[[1]]$weight.true[,1], color="red", size=1)+
  # geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
  labs(x="1st group's weight")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(c(0, 1))
p6

# df_new = data.frame(matrix(0, nrow=500, ncol=2))
# df_new[,1] = seq(1,500)
# df_new[,2] = df[,1]
# colnames(df_new) = c("ind", "weight")
# ggplot(df_new, aes(x=weight))+
#   geom_histogram(aes(y=..count..), colour="black",fill="white", bins=50)+
#   geom_density(alpha=.2, fill="#FF6666")+
#   # geom_vline(xintercept = trueList.List.all[[1]]$weight.true[,1], color="red", size=1)+
#   # geom_vline(xintercept = mean(df[,1]), color="blue", linetype="dashed", size=1)+
#   labs(x="1st group's weight")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# set.seed(1234)
# df1 <- data.frame(
#   sex=factor(rep(c("F", "M"), each=200)),
#   weight=round(c(rnorm(200, mean=55, sd=5),
#                  rnorm(200, mean=65, sd=5)))
# )
# head(df1)
# # Histogram with density plot
# ggplot(df1, aes(x=weight)) + 
#   geom_histogram(aes(y=..density..), colour="black", fill="white")+
#   geom_density(alpha=.2, fill="#FF6666") 


##########################################
ggarrange(p1,p3,p5,
          p2,p4,p6,
          ncol=3, nrow=2, heights = c(1.1, 1.0))
ggarrange(p1, p3, p5, ncol=3, nrow=1)
