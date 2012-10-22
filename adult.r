library("arules")

adata <- read.table("adult.data", header=F, sep=",")

save(adata,file="adata.rda"); 
mydata <- load("adata.rda");
newdata <- get(mydata);



#finaldata <- as(newdata, "itemMatrix");



## find all closed frequent itemsets
closed <- apriori(newdata, parameter = list(target = "closed", support = 0.4))
## rule induction
rules <- ruleInduction(closed, newdata, confidence = 0.8, control = NULL)
summary(rules)
