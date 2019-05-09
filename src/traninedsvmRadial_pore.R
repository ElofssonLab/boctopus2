library(e1071)
args <- commandArgs(trailingOnly = TRUE)

testData <- read.table(args[1])
rowL     <- length(testData[1,])
testX    <- testData[,c(1:rowL)]

print("Start")
for (i in c(1:length(testData[,1])))
    print(testData[i,1])
print("End")

load("svmtrain42_porefacing.Rdata")
print("Start")
#tempResult <- predict(svmtrain36_porefacing, testX, decision.values = TRUE, probability = TRUE)
tempResult <- predict(svmtrain42_porefacing, testX, decision.values = TRUE, probability = TRUE)

for (i in c(1:length(tempResult)))
    print(tempResult[i])
print("End")

print("Start")
for (i in c(1:length(tempResult)))
    print(attr(tempResult, "decision.values")[i])
print("End")

print("Prob")
print(attr(tempResult, "probabilities"))
print("End")
