library(mht)
df = read.csv("tmpout.csv")
df = df[,-1]
df = df[,colnames(df)!="system"]
X = as.matrix(df[,-1])
y = as.matrix(df[,1])

X = scale(X)

mod = bolasso(X, y, m=100, mu=round(10**seq(0,-2,-0.2),digit=3), probaseuil=0.9)
mod$ind
mod$frequency
write.csv(mod$frequency, "bolasso_frequency.csv")
