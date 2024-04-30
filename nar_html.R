library(Matrix)
library(foreach)
library(doParallel)
library(MASS)
library(xtable)
library(geosphere)
library(dplyr)
library(GNAR)
library(ggplot2)

sales_df = readr::read_csv("./sales.csv")

sales_df$date = as.Date(sales_df$date)
head(sales_df,5)


# Determine the number of unique products and time periods
N <- length(unique(sales_df$product_id))
T <- length(unique(sales_df$date))
print(T)
print(N)
#Could consider centering the data
Phi <- W1 <- matrix(0, nrow = N, ncol = N)  # Initialize Phi and W1 matrices
# Calculate W1

# Calculate Phi
for(i in 1:N){
  for (j in 1:N){
    # Condition: Check if the price between row i and j is less than 5
    if (abs(sales_df$price[i] - sales_df$price[j]) < 5) {
      Phi[i,j] = 1
    }
  }
}

# Set Phi as W1 for simplicity
Phi = W1

# Normalize Phi
Phi <- Phi / rowSums(Phi)


# Normalize W1
W1 <- W1 / rowSums(W1)

# Use the last 20 days as test data
X_train = matrix(0,N, T - 20)
Y_train = list()
for (i in 1:(T-20)){
  X_train[,i] = df_sales[with(df_sales, date==unique(df_sales$date)[i]),]$sales
  Y_train[[i]] = df_sales[with(df_sales, date==unique(df_sales$date)[i]),][,6:9]
}
T_train = dim(X_train)[2]-1
X_test = matrix(0,N, 21)
Y_test = list()
for (i in 1:21){
  X_test[,i] = df_sales[with(df_sales, date==unique(df_sales$date)[T-21+i]),]$sales
  Y_test[[i]] = df_sales[with(df_sales, date==unique(df_sales$date)[T-21+i]),][,6:9]
}
T_test = dim(X_test)[2]-1

gcv = function(X_train, Y_train, N, T_train, Sum1, Sum2, lambda_a = c(0.0001, 0.0005, 0.001, 0.005,0.01,0.05,0.1,0.5,1,5,10,50), lambda_b = c(0.0001, 0.0005, 0.001, 0.005,0.01,0.05,0.1,0.5,1,5,10,50)){
  l1 = length(lambda_a)
  l2 = length(lambda_b)
  gcv_matrix = matrix(rep(0, l1*l2), l1,l2)
  for (j in 1:l1){
    for (k in 1:l2){
      lambda1 = lambda_a[j]
      lambda2 = lambda_b[k]
      SX <- chol2inv(chol(Sum1+diag(c(rep(T_train*lambda1,N),rep(T_train*lambda2,N),rep(0,4)))))
      beta <- as.matrix(SX%*%Sum2)
      err = tr = 0
      for (i in 1:T_train){
        X = X_train[,i]
        Y = Y_train[[i]]
        Z = as.matrix(cbind(diag(X), diag(as.vector(W1%*%X)), Y))
        Xt = X_train[,i+1]
        err = err + norm(Xt- Z%*% beta,type = "F")^2
        tr = tr + N - sum(diag(Z %*% SX %*% t(Z)))
      }
      gcv_matrix[j,k] = 1/T_train*err/(1/T_train*tr)^2
    }
  }
  l_a = lambda_a[which(gcv_matrix == min(gcv_matrix), arr.ind = TRUE)[1]]
  l_b = lambda_b[which(gcv_matrix == min(gcv_matrix), arr.ind = TRUE)[2]]
  beta <- as.matrix(chol2inv(chol(Sum1+diag(c(rep(T_train*l_a,N),rep(T_train*l_b,N),rep(0,4)))))%*%Sum2)
  return(list(l_a = l_a, l_b = l_b, beta = beta))
}

Sum_ZTZ = 0
Sum_ZTX = 0

for (i in 1:T_train){
  X = X_train[,i]
  Y = Y_train[[i]]
  Z = as.matrix(cbind(diag(X), diag(as.vector(W1%*%X)), Y))
  Xt = X_train[,i+1]
  ZTZ <- t(Z)%*%Z
  ZTX <- t(Z)%*%Xt
  Sum_ZTZ = ZTZ + Sum_ZTZ
  Sum_ZTX =  ZTX + Sum_ZTX
  
}
ols = gcv(X_train, Y_train, N, T_train, Sum_ZTZ, Sum_ZTX)


mse_ols = 0

for (i in 1:T_test){
  X = X_test[,i]
  Y = Y_test[[i]]
  Z = as.matrix(cbind(diag(X), diag(as.vector(W1%*%X)), Y))
  Xt = X_test[,i+1]
  mse_ols = mse_ols + norm(Xt - Z%*%ols$beta,type = "F")^2
}
mse_ols = mse_ols/N/T_test

print(paste("OLS:",round(mse_ols,4)))


