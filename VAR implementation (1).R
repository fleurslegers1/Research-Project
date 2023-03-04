library(tidyverse)
library(lubridate)
library(vars)
library(nsprcomp)
library(modelr)
library(ggpubr)
library(zoo)
library(pense)
library(gridExtra)
library(dplyr)
library(parallel)
library(nsprcomp)
library(furrr)
library(matlib)
library(glmnet)
library(caret)
library(MASS)
library(glmnetUtils)


plan(multisession, workers = 8)

suppressMessages({
  
load_data <- function(data,
                      log_,
                      moving_av,
                      start_train,
                      end_train,
                      start_pred,
                      end_pred,
                      horizon_){
  # This function loads the training and test set
  # Options are: linearly interpolated, constant interpolated and modelled
  # Optional to take log of values and Moving Averages
  
  training_dates <- seq(start_train, end_train, by = 1)
  test_dates <- seq(start_pred, end_pred - horizon_, by = 1)
  
  if(data == "int_linear"){
    df <- read_csv("./R/Data sets/df_linear_N12.csv", col_names = TRUE)
    df <- subset(df, select = -1)
    colnames(df)[3] <- "Value"
  }
  
  else if(data == "int_constant"){
    df <- read_csv("./R/Data sets/df_constant_N12.csv", col_names = TRUE)
    df <- subset(df, select = -1)
    colnames(df)[3] <- "Value"
  }
  
  else if(data == "modelled"){
    df <- read_csv("./R/Data sets/df_modelled.csv", col_names = TRUE)
    df <- subset(df, select = -1) %>% dplyr::select(RWZI, Datum, Value)
    df <- df %>% mutate(Value = as.numeric(Value))
  }
  
  else{
    print("wrong input for data")
  }
  
  # Take 10-log of values
  if(log_){
    df$Value <- log10(df$Value)
  }
  
  # Remove STSs that are problematic (closed down, extremely high values)
  df <- df %>% filter(!(RWZI %in% c("WOERDEN", "KATWOUDE", "ECK-EN-WIEL", "VALBURG", "LIENDEN", "VEENENDAAL")))
  
  # Calculate Moving Average
  if(moving_av > 1){
    for(RWZI in unique(df$RWZI)){
      new_value <- rollmean(df[df$RWZI == RWZI, "Value"], moving_av, na.pad = TRUE, align = "right")
      df[df$RWZI == RWZI, "Value"] <- new_value
    }
  }
  
  # Remove NA's introduced by moving average
  df <- df[complete.cases(df),]
  
  # Select correct dates
  df_training <- df %>% filter(Datum %in% training_dates)
  df_test <- df %>% filter(Datum %in% test_dates)
  
  return(list(df = df, df_training = df_training, df_test = df_test))
}  

fill_adjacency <- function(rwzi_s,
                           n_rwzi,
                           type,
                           self,
                           n_col,
                           n_PC,
                           df_training){
  # based on the type of model in the input
  
  if(type == "PCA"){
    # create empty adjacency matrix
    Adj <- matrix(data = 0, nrow = n_rwzi, ncol = n_rwzi)
    
    # transform df to wide matrix
    df_pca <- subset(pivot_wider(df_training, names_from = RWZI, values_from = last(colnames(df_training))), select = -1)
    df_pca <- df_pca[complete.cases(df_pca),]
    
    # perform PCA
    # one principal component is calculated,
    # the n_col stps with highest coefficients are selected 
    PCs <- nsprcomp(df_pca, retx = TRUE, nrestart = 5, em_moving_avxiter = 300)
    s_dev <- PCs$sdev
    total_s_dev <- sum(s_dev)
    s_dev <- s_dev / total_s_dev
    PCs <- PCs$rotation
    
    s_dev_n_PC <- 0
    for(i in 1:n_PC){
      s_dev_n_PC <- s_dev_n_PC + s_dev[i]
    }
    print(paste0("Total proportion of explained variance by ", n_PC ," first PCs:", s_dev_n_PC))
    
    # get principal variables from PCs and their explained proportion of variance
    importance <- c()
    for(i in 1:n_rwzi){
      for(j in 1:n_PC){
        importance[i] <- s_dev[j] * abs(PCs[i,j])
      }
    }
    PVs <- as_tibble(importance) %>% mutate(index = seq(1,n_rwzi,by=1)) %>% arrange(desc(value)) %>% dplyr::select(index) %>% pull(1)
    PVs <- PVs[1:n_col]
    
    principal_rwzis <- colnames(df_pca)[PVs]
    print("Principal STPs:")
    print(principal_rwzis)
    
    # set columns of principal stps to 1
    Adj[, PVs] <- 1
    
    # if self == TRUE, set diagonal to 1
    if(self){
      diag(Adj) <- 1
    }
  }
  
  else if(type == "KNN"){
    # load df of distances between stps
    Dist <- read.delim("./R/Data sources/RWZI Distances.csv", sep = ";", row.names = 1, header = TRUE)
    rownames(Dist) <- gsub(" ", "-", toupper(rownames(Dist)))
    colnames(Dist) <- rownames(Dist)
    Dist <- Dist %>% dplyr::select(sort(names(.)))
    Dist <- as.matrix(Dist[, rwzi_s])
    
    # create empty adjacency matrix
    Adj <- matrix(data = 0, nrow = n_rwzi, ncol = n_rwzi)
    
    # choose n_col closest and set corresponding value in Adj to 1
    for(RWZI in 1:n_rwzi){
      indices <- match(names(Dist[rwzi_s[RWZI], order(Dist[rwzi_s[RWZI],])[2:(n_col + 1)]]), rwzi_s) # which column in Dist
      Adj[RWZI, indices] <- 1
    }
    
    if(self){
      diag(Adj) <- 1
    }
  }
  
  else if(type == "AR" | type == "LAST"){
    # create adjacency matrix with only ones on the diagonal 
    Adj <- matrix(data = 0, nrow = n_rwzi, ncol = n_rwzi)
    diag(Adj) <- 1
  }
  
  else if(type == "CORR"){
    # read in df with correlations between lagged time series of STPs,
    # with lags up to 14 days
    Dist <- read.delim("./R/Data sources/RWZI Correlations.csv", sep = ";", row.names = 1, header = TRUE)
    
    # only select the columns with p lag
    Dist <- Dist[,1:(p*n_rwzi)]
    Dist <- as.matrix(Dist)
    
    # create empty adjacency matrix
    Adj <- matrix(data = 0, nrow = n_rwzi, ncol = n_rwzi)
    
    # for each row, sort correlations and select n_col for adjacency matrix
    for(RWZI in 1:n_rwzi){
      
      # do a lot of stuff to get from the column names the column of Adj to put a 1
      names_of_columns <- names(Dist[rwzi_s[RWZI], order(Dist[rwzi_s[RWZI],])[1:n_col]])
      for(i in 1:14){
        names_of_columns <- gsub(paste0(".L", i), "", names_of_columns)
      }
      names_of_columns <- gsub("[.]", "-", names_of_columns)
      indices <- match(names_of_columns, rwzi_s)
      Adj[RWZI, indices] <- 1
    }
    
    if(self){
      diag(Adj) <- 1
    }
  }

  else{     
    # create adjacency matrix with all ones
    Adj <- matrix(data = 1, nrow = n_rwzi, ncol = n_rwzi)
  }

  return(Adj)
}

my_regression <- function(df){
  # this function finds x given A and b such that Ax = b
  # using the Moore-Penrose inverse (and SVD)
  
  # First column of df contains the target vector b
  b <- as.matrix(df[,1], ncol = 1)
  df$constant <- 1
  
  # Rest of df contains the matrix A
  A <- as.matrix(df[,-1])
  
  # Calculates reciprocal condition number
  #if(condition_number){
  #  condition_number <<- rcond(A)
  #}
  
  # x = pseudo-inverse(A) * b
  SVD <- svd(A, nu = nrow(A), nv = ncol(A))
  pseudo_Sigma <- matrix(data = 0, nrow = ncol(A), ncol = nrow(A))
  diag(pseudo_Sigma)[1:length(SVD["d"]$d)] <- 1/SVD["d"]$d      
  
  pseudo_A <- SVD$v %*% pseudo_Sigma %*% t(SVD$u)
  
  x <- pseudo_A %*% b
  
  return(x[,1])
}

regression <- function(type,
                       df_training,
                       p,
                       n_rwzi,
                       A_){
  # this function fills the coefficient matrix W
  # by performing my_regression for each STP (row of W)
  
  if(type == "LAST"){
    W <- matrix(data = 0, nrow = n_rwzi, ncol = n_rwzi + 1)
    diag(W) <- 1
  }
  
  else{
  
  # transform data to wide matrix
  df_temp <- subset(pivot_wider(df_training, names_from = RWZI, values_from = Value), select = -1)
  
  # add shifted values
  ncol <- ncol(df_temp)
  colnames <- colnames(df_temp)
  for(lag in 1:p){
    for(col in 1:ncol){
      temp <- c(rep(NA, lag), head(pull(df_temp, col), -lag))
      df_temp <- df_temp %>% add_column(name = temp, .name_repair = "unique")
    }
  }
  for(lag in 1:p){
    colnames <- c(colnames, paste0(colnames[1:ncol], ".L", lag))
  }
  colnames(df_temp) <- colnames
  
  # create coefficient matrix W
  W <- matrix(data = 0, nrow = n_rwzi, ncol = p * n_rwzi + 1)
  W <- as.data.frame(W)
  colnames(W) <- c(tail(colnames(df_temp), -n_rwzi), "const")

    #### LINEAR REGRESSION ####
    for(row in 1:n_rwzi){
      
      # find columns to use (from A)
      indices <- row
      for(i in 1:p){
        indices <- c(indices,  which(A_[row,] > 0) + i * n_rwzi)
      }
      
      # make df with only these columns (indices) and delete NA-rows (first p)                             
      df_temp_2 <- df_temp[complete.cases(df_temp), indices]                      
      colnames(df_temp_2) <- gsub("-", "", colnames(df_temp_2))
      
      # do least squares for each row of W
      coefficients <- my_regression(df_temp_2)
      
      for(coeff in 1:(length(coefficients) - 1)){
        W[row, indices[coeff + 1] - n_rwzi] <- coefficients[coeff]
      }
      
      W[row, ncol(W)] <- tail(coefficients, 1)
    }
  
  
  # Check if entire columns are NA
  for(col in ncol(W)){
    if(sum(is.na(W[,col])) == nrow(W)){
      print("column with only NAs found")
    }  
  }
  
  # Replace NA with 0
  for(row in 1:nrow(W)){
    for(col in 1:ncol(W)){
      if(is.na(W[row, col])){
        W[row, col] <- 0
      }
    }
  }
  }
  
  return(W)
}

RMSE <- function(x, y){
  # calculated the RMSE for two vectors of the same length
  
  if(length(x)!=length(y)){
    print("Error in RMSE: actual and predicted are of different size")
    return(NA)
  }
  rmse <- sum((x - y)^2, na.rm = TRUE)
  rmse <- 1.0/length(x) * rmse
  rmse <- rmse**(1/2)
  return(rmse)
}

predict <- function(type,
                    data,
                    log_,
                    p,
                    horizon_,
                    moving_av,
                    self,
                    n_col,
                    start_train,
                    end_train,
                    start_pred,
                    end_pred,
                    plotting_,
                    n_rwzi,
                    rwzi_s,
                    W,
                    df_training,
                    df_test,
                    n_PC){
  # this function makes predictions using the coefficient matrix W
  # both on the training and the test set,
  # it calculates the RMSE and makes plots
  # it also calculates the average of the STPs, and makes plots for these results as well
  
  #print(paste("p = ", p))
  #print(paste("n_col = ", n_col))
  
  #### PREDICT ON TRAINING DATA ####
  
  # add type to data
  df_temp <- df_training
  df_temp$type <- "actual"
  
    for(datum in as.character(tail(unique(df_training$Datum), -p))){
      
      # create vector with correct input measurements
      x <- c()
      
      for(j in 1:p){
        x <- rbind(x, df_training[df_training$Datum == as.Date(datum, "%Y-%m-%d") - j, "Value"])
      }
      
      # add constant
      x <- c(pull(x, 1),1)
      
      # ERROR CHECK
      if(length(x) != ncol(W)){
        print(datum)
      }
      
      # get predictions
      y <- as.matrix(W) %*% as.numeric(x)
      
      # add predictions to df
      columns <- cbind(rwzi_s, rep(datum, n_rwzi), y[,1], rep("predicted", n_rwzi))
      colnames(columns) <- colnames(df_temp)
      df_temp <- rbind(df_temp, columns)
  
    }
  
  # print root mean squared error
  temp <- full_join(df_temp,df_temp, by = c("RWZI", "Datum"))
  RMSE_training <<- RMSE(as.numeric(temp$Value.x), as.numeric(temp$Value.y))
  #print(paste("RMSE on training data: ", RMSE_training))
  
  #### PREDICT ON TEST DATA ####
  
  # add type to data
  df_temp_2 <- rbind(df_training, df_test)
  df_temp_2$type <- "actual"
  df_temp_2$which_day <- NA
  
  # for every day in start_pred to end_pred - horizon_, 
  # make horizon_ predictions
  df_for_RMSE <- c()
  
  which_prediction <- 1
  for(datum in as.character(seq(start_pred, end_pred - horizon_, by = 1))){
    
    # get actual values that will be used in the predictions 
    df_temp_3 <- df_temp_2[df_temp_2$Datum < as.Date(datum, "%Y-%m-%d") & df_temp_2$Datum >= as.Date(datum, "%Y-%m-%d") - p,]
    
    # save the values of the previous date
    y_prev = pull(df_temp_3[df_temp_3$Datum == as.Date(datum, "%Y-%m-%d") - 1, "Value"], 1)
    
    which_hor <- 1
    
    for(datum_in_horizon in as.character(seq(as.Date(datum, "%Y-%m-%d"), as.Date(datum, "%Y-%m-%d") + horizon_ - 1, by = 1))){
      
      # create vector with correct input measurements
      x <- c()
      
      for(j in 1:p){
        # get p previous values from df_temp_3
        x <- rbind(x, df_temp_3[df_temp_3$Datum == as.Date(datum_in_horizon, "%Y-%m-%d") - j, "Value"]) # CHECK
      }
      
      # add constant
      x <- c(pull(x, 1),1)
      
      # ERROR CHECK
      if(length(x) != ncol(W)){
        print(datum)
      }
      
      # get predictions
      y <- as.matrix(W) %*% as.numeric(x)
      
      # add restrictions:
      # difference between two consecutive predictions must be at most 0.759274950540281
      # 10^12 <= y <= 
      for(i in 1:length(y[,1])){
        if(!is.numeric(y[i,1])){
          View(y[i,1])
          print(datum_in_horizon)
          print(which_hor)
        }
        if(!is.numeric(y_prev[i])){
          View(y_prev[i])
          print(datum_in_horizon)
          print(which_hor)
        }
        if(y[i,1] - y_prev[i] > 0.76){
          y[i,1] <- y_prev[i] + 0.76
        }
        else if(y[i,1] - y_prev[i] < -0.76){
          y[i,1] <- y_prev[i] - 0.76
        }
        if(y[i,1] < 12){
          y[i,1] <- 12
        }
        if(y[i,1] > 15.25){
          y[i,1] <- 15.25
        }
      }
      
      y_prev <- y
      
      # add predictions to df
      columns <- cbind(rwzi_s, rep(datum_in_horizon, n_rwzi), y[,1], rep("predicted", n_rwzi), rep(which_hor, n_rwzi))
      colnames(columns) <- colnames(df_temp_2)
      df_temp_3 <- rbind(df_temp_3, columns)
      
      if(length(pull(df_temp_2[df_temp_2$Datum == as.Date(datum_in_horizon, "%Y-%m-%d") & df_temp_2$type == "actual", "Value"],1)) == n_rwzi){
        columns <- cbind(rwzi_s, rep(datum_in_horizon, n_rwzi), y[,1], pull(df_temp_2[df_temp_2$Datum == as.Date(datum_in_horizon, "%Y-%m-%d") & df_temp_2$type == "actual", "Value"],1), rep(which_hor, n_rwzi), rep(which_prediction, n_rwzi))
        df_for_RMSE <- rbind(df_for_RMSE, columns)
      }
      
      which_hor <- which_hor + 1
    }
    
    which_prediction <- which_prediction + 1
  }
  
  df_for_RMSE <- as.data.frame(df_for_RMSE)
  colnames(df_for_RMSE) <- c("RWZI", "Datum", "Predicted", "Actual", "Day_of_prediction", "Index_of_prediction")
  
  # save df_for_RMSE
  write.csv(as.data.frame(df_for_RMSE), file = paste0("./R/Outputs VAR/", type, "/n_PC = ", n_PC, "p = ", p, ", RMSE per prediction.cvs"))
  
  # if plotting_ = TRUE, also print/save the RMSEs for different horizon-days
  if(plotting_){
  RMSE_per_day <- c()
    for(i in 1:horizon_){
      actual <- as.numeric(df_for_RMSE[df_for_RMSE$Day_of_prediction == i, "Actual"])
      predicted <- as.numeric(df_for_RMSE[df_for_RMSE$Day_of_prediction == i, "Predicted"])
      RMSE_horizon <- RMSE(actual, predicted)
      RMSE_per_day <- c(RMSE_per_day, RMSE_horizon)
      print(paste("RMSE on testing data for", i, "days ahead: ", RMSE_horizon))
    }
  write.csv(as.data.frame(RMSE_per_day), file = paste0("./R/Outputs VAR/", type, "/n_PC = ", n_PC, "p = ", p, ", ", moving_av, self, n_col, ", RMSE per day.cvs"))
  }
  
  RMSE_test_rwzi <- list()
  
  for(rwzi in unique(df_for_RMSE$RWZI)){
    x <- df_for_RMSE %>% as_tibble() %>% filter(RWZI == rwzi) %>% dplyr::select(Predicted)  %>% pull(1) %>% as.numeric()
    y <- df_for_RMSE %>% as_tibble() %>% filter(RWZI == rwzi) %>% dplyr::select(Actual)  %>% pull(1) %>% as.numeric()
    rmse <- RMSE(x,y)
    RMSE_test_rwzi <- rbind(RMSE_test_rwzi, cbind("RWZI" = rwzi, "RMSE_test" = rmse))
  }
  
  RMSE_test <- RMSE(as.numeric(df_for_RMSE$Predicted), as.numeric(df_for_RMSE$Actual))
  #print(paste("RMSE on testing data: ", RMSE_test))
  
  #### MAKE PLOTS ####
  if(plotting_){
  
  pdf(file = paste0("./R/Figures/", type, "/n_PC = ", n_PC, 
                   ", p = ", p, 
                   ", MA = ", moving_av,
                   ", n_col = ", n_col,
                   ", log = ", log_,
                   ", self = ", self,
                   ", data = ", data_type,
                   ".pdf"))
  
  for(RWZI in rwzi_s){
    
    if(log_){
      ylab <- "Log(N12 per 100000)"
    }
    else{
      ylab <- "N12 per 100000"
    }
    
    # Make plot for training data + prediction on training data
    # add real data as thick red dots
    if(moving_av == 0){
      plot_training <- ggplot(data = df_temp[df_temp$RWZI == RWZI,], 
                              mapping = aes(x = Datum, "%Y-%m-%d", y = as.numeric(Value), color = type))+
        labs(title = paste(type, "on training data for ", RWZI), 
             subtitle = paste("p = ", p, 
                              ", data_type = ", data_type, 
                              ", log = ", log_, 
                              ", MA = ", moving_av),
             x = "Datum",
             y = ylab)+
        geom_point()+
        geom_line()
      #+
      #  geom_point(data = df_og, mapping = aes(x = Datum, "%Y-%m-%d", y = as.numeric(Value), color = "red", shape = 6))
    }
    
    else{
    plot_training <- ggplot(data = df_temp[df_temp$RWZI == RWZI,], 
                            mapping = aes(x = Datum, "%Y-%m-%d", y = as.numeric(Value), color = type))+
                            labs(title = paste(type, "on training data for ", RWZI), 
                                 subtitle = paste("p = ", p, 
                                                     ", data_type = ", data_type, 
                                                     ", log = ", log_, 
                                                     ", MA = ", moving_av),
                                 x = "Datum",
                                 y = ylab)+
      geom_point()+
      geom_line()
    }
    # Make data frame in form suitable for ggplot
    df_for_RMSE_2 <- df_for_RMSE %>% pivot_longer(cols = c(Actual, Predicted),
                                                names_to = "type",
                                                values_to = "Value")
    
    df_for_RMSE_2[df_for_RMSE_2$type == "Predicted","type"] <- paste0(pull(df_for_RMSE_2[df_for_RMSE_2$type == "Predicted","type"],1), pull(df_for_RMSE_2[df_for_RMSE_2$type == "Predicted","Index_of_prediction"],1))
    
    
    # Select some lines to plot
    indices <- c(1, seq(2, as.numeric(max(df_for_RMSE_2$Index_of_prediction)), by = 6), as.numeric(max(df_for_RMSE_2$Index_of_prediction)))
    df_for_RMSE_2 <- df_for_RMSE_2[df_for_RMSE_2$Index_of_prediction %in% indices, c("RWZI", "Datum", "type", "Value")] 
    
    
    # Add some training data
    df_to_add <- df_training[df_training$Datum > start_pred - horizon_,]
    df_to_add$Datum <- as.character(df_to_add$Datum)
    df_to_add$type <- "Actual"
      
    df_for_RMSE_2 <- rbind(df_for_RMSE_2, df_to_add) 
    df_for_RMSE_2 <- df_for_RMSE_2 %>% mutate(Datum = as.Date(Datum, "%Y-%m-%d"))
    df_for_RMSE_2 <- df_for_RMSE_2 %>% mutate(Value = as.numeric(Value))
    
    # Add last training data point to prediction
    for(i in indices){
      # get starting date - 1 of this prediction
      start_prediction_date <- df_for_RMSE_2 %>% filter(type == paste0("Predicted",i)) %>% dplyr::select(Datum) %>% pull(1) %>% min() - 1
      char_date <- as.character(start_prediction_date)
      
      # add this training data to prediction i
      columns <- cbind(rwzi_s, rep(char_date, n_rwzi), rep(paste0("Predicted",i), n_rwzi), pull(df_for_RMSE_2[df_for_RMSE_2$type == "Actual" & df_for_RMSE_2$Datum == start_prediction_date, "Value"],1))
      colnames(columns) <- c("RWZI", "Datum", "type", "Value")
      df_for_RMSE_2 <- rbind(df_for_RMSE_2, columns)
    }
    
    # Make plots
    plot_testing <- ggplot(data = df_for_RMSE_2[df_for_RMSE_2$Datum > start_pred - horizon_ & df_for_RMSE_2$RWZI == RWZI,], 
                           mapping = aes(x = Datum, y = as.numeric(Value), color = type))+
                           labs(title = paste(type, "on test data for ", RWZI), 
                                subtitle = paste("p = ", p, 
                                                 ", data_type = ", data_type, 
                                                 ", log = ", log_, 
                                                 ", MA = ", moving_av),
                                x = "Date",
                                y = ylab)+
      geom_point()+
      geom_line()+
      theme(legend.position = "none")
    
    # Save plots
    print(ggarrange(plot_training, plot_testing, nrow = 2, ncol = 1))
    
  }
  
  #### Calculate Average of NL ####
  inwoners_per_rwzi_2 <- read_delim(".R/Data sources/inwoners per rwzi_2.csv",
                                    delim = ";",
                                    escape_double = FALSE,
                                    trim_ws = TRUE)
  
  N <- sum(inwoners_per_rwzi_2$Inwoners)
  inwoners_per_rwzi_2$Inwoners <- inwoners_per_rwzi_2$Inwoners / N
  
  df_for_RMSE_2 <- df_for_RMSE_2 %>% mutate(Value = as.numeric(Value)) %>% left_join(inwoners_per_rwzi_2, by = "RWZI") %>% mutate(Value = Value * Inwoners)
  df_for_RMSE_2 <- df_for_RMSE_2 %>% distinct() %>% aggregate(Value ~ Datum + type, sum)
  
  df_temp$Value <- as.numeric(df_temp$Value)
  df_temp <- df_temp %>% left_join(inwoners_per_rwzi_2, by = "RWZI") %>% mutate(Value = Value * Inwoners)
  df_temp <- aggregate(Value ~ Datum + type, data = df_temp, sum)
  #
  df_temp_2 <- df_for_RMSE_2[complete.cases(df_for_RMSE_2),]
  #df_temp_2 <- df_temp_2 %>% mutate(Predicted = as.numeric(Predicted), Actual = as.numeric(Actual), Datum = as.Date(Datum))
  #df_temp_2 <- aggregate(Value ~ Datum + type, data = df_temp_2, mean)
  #
  plot_training <- ggplot(data = df_temp,
                          mapping = aes(x = Datum, y = Value, color = type))+
    labs(x = "Date",
         y = bquote(Log["10"] (load) ~("10"^{"-5"} ~persons^{"-1"})))+
    theme(text = element_text(size = 13))+
    geom_point()+
    geom_line()
  
  plot_testing <- ggplot(data = df_temp_2[df_temp_2$Datum > start_pred - horizon_, ],
                         mapping = aes(x = Datum, y = Value, color = type))+
    labs(x = "Date",
         y = bquote(Log["10"] (load) ~("10"^{"-5"} ~persons^{"-1"})))+
    geom_point()+
    geom_line()+
    theme(legend.position = "none",
          text = element_text(size = 13))
  
  print(ggarrange(plot_training, plot_testing, nrow = 2, ncol = 1))
  
  dev.off()
  }
  }
  
  #### RETURN ####
  
  RMSE_test <<- RMSE_test
  print(RMSE_test)
  RMSE_test_rwzi <<- RMSE_test_rwzi
  return(list(p = p, n_col = n_col, W = W, RMSE_training = RMSE_training, RMSE_test = RMSE_test, RMSE_test_rwzi = RMSE_test_rwzi))


predictive_power_plot <- function(){
  
  # this function makes a plot of the RMSE per datum, given some W
  # it also detects outliers and reports these
  # add type to data
  df_pp <- rbind(df_training, df_test)
  df_pp_2 <- df_pp
  df_pp_2$type <- "actual"
  
  # predict in-training errors
  for(datum in as.character(tail(unique(df_pp$Datum), -p))){
    
    # create vector with correct input measurements
    x <- c()
    
    for(j in 1:p){
      x <- rbind(x, df_pp[df_pp$Datum == as.Date(datum, "%Y-%m-%d") - j, "Value"])
    }
    
    # add constant
    x <- c(pull(x, 1),1)
    
    # get predictions
    y <- as.matrix(W) %*% as.numeric(x)
    
    # add predictions to df
    columns <- cbind(rwzi_s, rep(datum, n_rwzi), y[,1], rep("predicted", n_rwzi))
    colnames(columns) <- colnames(df_pp_2)
    df_pp_2 <- rbind(df_pp_2, columns)
    
  }
  
  # print root mean squared error
  temp <- full_join(df_pp_2[df_pp_2$type == "actual",c(1,2,3)], df_pp_2[df_pp_2$type == "predicted",c(1,2,3)], by = c("RWZI", "Datum"))
  temp <- temp[complete.cases(temp),]
  temp <- temp %>% mutate(Value.x = as.numeric(Value.x), Value.y = as.numeric(Value.y))
  temp$diff <- abs(temp$Value.x - temp$Value.y)
  temp <- temp %>% dplyr::select(RWZI, Datum, diff)
  
  pdf(file = paste0("./R/Figures/", type, "/Analysis/predictive_power_plot, MA = ", MA, 
                    ", p = ", p,
                    "n_col = ", n_col, ".pdf"))
  
  df_pp <- df_pp %>% filter(Datum %in% temp$Datum)
  df_pp <- df_pp %>% mutate(Value = as.numeric(Value))
  
  for(RWZI in rwzi_s){
    
    # Plot trainingsdata
    plot_training <- ggplot(data = df_pp[df_pp$RWZI == RWZI, ], mapping = aes(x = Datum, y = Value))+
      labs(title = paste("Data for ", RWZI), 
           x = "Datum",
           y = "Log(N12 per 100000)")+
      geom_point()+
      geom_line()
    
    # Make plot for training data + prediction on training data
    plotje <- ggplot(data = temp[temp$RWZI == RWZI,], mapping = aes(x = Datum, y = as.numeric(diff)))+
      labs(title = paste("Prediction error per date for ", RWZI), 
           subtitle = paste("p = ", p, 
                            ", MA = ", moving_av,
                            ", n_col = ", n_col),
           x = "Datum",
           y = "prediction error")+
      geom_point()+
      geom_line()
    
    print(ggarrange(plot_training, plotje, nrow = 2, ncol = 1))
  
  }
  
  # Calculate average of NL
  temp <- aggregate(diff ~ Datum, data = temp, mean)
  df_pp <- aggregate(Value ~ Datum, data = df_pp, mean)
  
  plot_training <- ggplot(data = df_pp, mapping = aes(x = Datum, y = Value))+
    labs(title = "Data, Average of NL", 
         x = "Datum",
         y = "Log(N12 per 100000)")+
    geom_point()+
    geom_line()
  
  plotje <- ggplot(data = temp, mapping = aes(x = Datum, y = diff))+
    labs(title = paste("Prediction oerror per date, average of NL"),
         subtitle = paste("p = ", p, 
                          ", MA = ", moving_av,
                          ", n_col = ", n_col),
         x = "Datum",
         y = "prediction error")+
    geom_point()+
    geom_line()
  
  print(ggarrange(plot_training, plotje, nrow = 2, ncol = 1))
  
  dev.off()
  
  return(plotje)
}

# choose type from: "VAR", "PCA", "KNN", "sparse"
VAR <- function(type = "KNN",
                data = "int_linear",
                log_ = TRUE,
                p = 1,
                horizon_ = 3,
                moving_av = 1,
                self = TRUE,
                n_col = 10,
                start_train = as.Date("2022-06-01", "%Y-%m-%d"),
                end_train = as.Date("2022-08-21", "%Y-%m-%d"),
                start_pred = as.Date("2022-08-22", "%Y-%m-%d"),
                end_pred = as.Date("2022-10-01", "%Y-%m-%d"),
                plotting_ = TRUE,
                n_PC = 3){
  
  if(type == "LAST"){
    p <- 1
  }
  
  #### LOAD DATA ####
  # load data of types "int_linear", "int_constant" and "weekly" 
  # if log, perform log10-transformoving_avtion
  # if moving_av > 1, use moving averages over previous measurements
  results_load_data <- load_data(data, log_, moving_av, start_train, end_train, start_pred, end_pred, horizon_)
  df <- results_load_data$df
  df_training <- results_load_data$df_training
  df_test <- results_load_data$df_test
  
  rwzi_s <- unique(df$RWZI)
  n_rwzi <- length(rwzi_s)
  
  #### MAKE ADJACENCY ####
  # moving_avke adjacency matrix based on type of VAR
  # for types "AR", VAR", "PCA", "KNN"
  A_ <- fill_adjacency(rwzi_s, n_rwzi, type, self, n_col, n_PC, df_training)

  
  #### LM REGRESSION ####
  # perform regression
  W <- regression(type, df_training, p, n_rwzi, A_)

  
  #### PREDICT TRAINING DATA ####
  # predicts on training and test data
  # makes plots
  # calculates RMSEs
  # calculates averages for NL
  output <- predict(type, 
                    data, 
                    log_, 
                    p, 
                    horizon_,
                    moving_av, 
                    self, 
                    n_col, 
                    start_train, 
                    end_train,
                    start_pred, 
                    end_pred,
                    plotting_, 
                    n_rwzi,
                    rwzi_s, 
                    W, 
                    df_training, 
                    df_test,
                    n_PC)
  
  return(output)
}

#### Functions for parameter loop ####

function_within_loop <- function(data_type, 
                                 moving_av, 
                                 type, 
                                 log_,
                                 self,
                                 start_train,
                                 end_train,
                                 start_pred,
                                 end_pred, 
                                 horizon_,
                                 plotting_,
                                 order_p,
                                 n_PC){
  
  print(paste0("p = ", order_p))
  
  if(type %in% c("KNN", "CORR", "PCA")){
    cols <- as.list(c(seq(1, 100, by = 10), seq(120, 240, by = 20)))
    output <- cols %>% map(function(x) {VAR(type, data_type, log_, order_p, horizon_ = horizon_, moving_av, self, x, 
                                            start_train, end_train, start_pred, end_pred, plotting_ = plotting_, n_PC = n_PC)})
  }
  
  else{
    output <- VAR(type, data_type, log_, order_p, horizon_ = horizon_, moving_av, self, 310, 
                  start_train, end_train, start_pred, end_pred, plotting_ = plotting_, n_PC = n_PC)
  }
  
  
  # output wegschrijven
  saveRDS(output, file = paste0("./R/Outputs VAR/", type, 
                                "/n_PC = ", n_PC, 
                                "p = ", order_p, 
                                ", moving_av = ", moving_av,
                                ", log = ", log_, 
                                ", self = ", self, 
                                ", training from ", start_train,
                                " to ", end_train,
                                ", start prediction on ", start_pred,
                                ", end prediction on ", end_pred,
                                ", horizon = ", horizon_,
                                ".RdS"))
  
  return(output)
}

find_optimal_parameter <- function(data_type = "modelled",
                                   moving_av = 0,
                                   type = "KNN",
                                   log_ = FALSE,
                                   self = TRUE,
                                   start_train = as.Date("2022-01-01", "%Y-%m-%d"),
                                   end_train = as.Date("2022-10-01", "%Y-%m-%d"),
                                   start_pred = as.Date("2022-10-02", "%Y-%m-%d"),
                                   end_pred = as.Date("2022-08-21", "%Y-%m-%d"),
                                   horizon_= 7,
                                   plotting_ = FALSE,
                                   n_PC = 3){

  results <- as.list(1:7) %>% future_map(function(x) {function_within_loop(data_type = data_type,                      
                                                       moving_av = moving_av, 
                                                       type = type, 
                                                       log_ = log_,
                                                       self = self,
                                                       start_train = start_train,
                                                       end_train = end_train,
                                                       start_pred = start_pred,
                                                       end_pred = end_pred,
                                                       horizon_ = horizon_,
                                                       plotting_ = plotting_,
                                                       order_p = x,
                                                       n_PC = n_PC)})
  
  return(results)
}

save_RMSE_per_date <- function(data_type = "int_linear",
                               moving_av = 7,
                               type = "PCA",
                               log_ = TRUE,
                               self = FALSE,
                               horizon = 7,
                               plotting_ = FALSE,
                               n_PC = 41){
  dates <- list(1:2)
  
  results <- dates %>% future_map(function(i) {
    
    # get start_train, end_train, start_pred, end_pred
    start_pred <- as.Date("2022-06-24") + i
    end_pred <- as.Date("2022-06-24") + i + 7
    start_train <- as.Date("2022-06-24") + i - 49
    end_train <- as.Date("2022-06-24") + i - 1
    
    # do VAR
    a <- VAR(data = "int_linear",
            moving_av = 7,
            type = "PCA",
            log_ = TRUE,
            p = 7,
            n_col = 60,
            self = FALSE,
            start_train = start_train,
            end_train = end_train,
            start_pred = start_pred,
            end_pred = end_pred,
            horizon_ = 7,
            plotting_ = FALSE,
            n_PC = 41)
    
    test_error <- a$RMSE_test
    
    return(test_error)
    
  })
  
  saveRDS(output, file = paste0("./R/Outputs VAR/test results.RdS"))
}

#### MAIN ####

# a <- VAR(
#   data = "modelled",
#   moving_av = 0,
#   type = "PCA",
#   log_ = FALSE,
#   p = 2,
#   n_col = 31,
#   self = FALSE,
#   start_train = as.Date("2022-01-01", "%Y-%m-%d"),
#   end_train = as.Date("2022-06-23", "%Y-%m-%d"),
#   start_pred = as.Date("2022-06-24", "%Y-%m-%d"),
#   end_pred = as.Date("2022-10-20", "%Y-%m-%d"),
#   horizon_ = 7,
#   plotting_ = TRUE,
#   n_PC = 4
# )
# 
# b <- a$RMSE_test_rwzi %>% as_tibble() %>% mutate(RWZI = unlist(RWZI), RMSE_test = as.numeric(unlist(RMSE_test)))
# write.csv(as.data.frame(b), file = "./R/Figures/PCA/modelled, results per RWZI.csv")


})