library(mht)
df = read.csv("tmpout.csv")
df = df[,-1] # remove index column in dataframe
df = df[,colnames(df)!="system"] # remove system tag if exists

X = as.matrix(df[,-1])
y = as.matrix(df[,1])

X = scale(X)

mod = bolasso(X, y, m=100, mu=round(10**seq(-0.4,-2,-0.1),digit=3), probaseuil=0.9)
mod$ind
mod$frequency
write.csv(mod$frequency, "bolasso_frequency.csv")
