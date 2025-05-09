setwd("~/Documents/GitHub/JLDispense/benchmarking")


library(ggplot2)


data = read.csv("benchmark.csv")

data= data[data$time<800,]

insts=c("P-200","P-100 Multichannel","PlateMaster","Nimbus","Cobra")
pistons=c(1,1,1,1,4)
channels=c(1,8,96,1,4)
piston_dict= as.list(setNames(pistons,insts))
channel_dict= as.list(setNames(channels,insts))
piston_vals=rep(0,nrow(data))
channel_vals=rep(0,nrow(data))
for (i in 1:nrow(data)){
  instruments = data$instruments[i]
  instruments = as.vector(strsplit(instruments,","))
  p = 0 
  c=0
  for (inst in instruments[[1]]){
    p = p + piston_dict[inst][[1]]
    c = c + channel_dict[inst][[1]]
  }
  piston_vals[i]= p
  channel_vals[i]=c
}

data$pistons= piston_vals
data$channels=channel_vals
model <- lm(time ~ I(wells^2)*pistons*channels,data=data)

summary(model)


plt = ggplot(data=data)+ 
  geom_point(aes(x=wells^3*pistons*channels,y=time))+
  theme_classic()
plt

data$pistons = piston_vals
