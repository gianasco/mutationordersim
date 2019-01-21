
#summary Statistics 

library(stringr) # to test functions 


#read data from test file 

data <- read.delim(“CosmicCompleteExport.tsv”,header = TRUE)



#finding minimum function
my.min <- function(x) {
  if(any(is.na(x))){     # remove NA or NAN from data
    x<- x[-which(is.na(x)==TRUE)]
  }
  sortArray<-sort(x)
  sortArray[1]
}




#finding maximum function
my.max <- function(x) {
  if(any(is.na(x))){     # remove NA or NAN from data
    x<- x[-which(is.na(x)==TRUE)]
  }
  sortArray<-sort(x) # order list in increment order 
  sortArray[length(sortArray)]
}

#finding mean function
my.mean <- function(x) {
  if(any(is.na(x))){     # to remove NA or NAN from data
    x<- x[-which(is.na(x)==TRUE)]
  }
  return(sum(x) / length(x))
}


#finding medain function
my.median <- function(x) {
  if(any(is.na(x))){     # to remove NA or NAN from data
    x<- x[-which(is.na(x)==TRUE)]
  }
  sorted.x <- sort(x)
  
  if (length(x) %% 2 == 0)
  {
    indices <- c(length(x) / 2, length(x) / 2 + 1)
    return(my.mean(sorted.x[indices]))
  }
  else
  {
    index <- ceiling(length(x) / 2)
    return(sorted.x[index])
  }
}

#finding varianc function
my.var <- function(x) {
  if(any(is.na(x))){      # to remove NA or NAN from data
    x<- x[-which(is.na(x)==TRUE)]
  }
  m <- my.mean(x)
  return(sum((x - m) ^ 2) / (length(x) - 1))
}

#finding standrad deviation
my.sd <- function(x) {
  if(any(is.na(x))){    # to remove NA or NAN from data
    x<- x[-which(is.na(x)==TRUE)]
  }
  return(sqrt(my.var(x)))
}


#testing custom functions

#expect_equal(my.min(data$rna_act_3utr)-min(data$rna_act_3utr,na.rm = TRUE),0)
#expect_equal(my.max(data$rna_act_3utr)-max(data$rna_act_3utr,na.rm = TRUE),0)
#expect_equal(my.mean(data$rna_act_3utr)-mean(data$rna_act_3utr,na.rm = TRUE),0)
#expect_equal(my.sd(data$rna_act_3utr)-sd(data$rna_act_3utr,na.rm = TRUE),0)
#expect_equal(my.var(data$rna_act_3utr)-var(data$rna_act_3utr,na.rm = TRUE),0)
#expect_equal(my.median(data$rna_act_3utr)-median(data$rna_act_3utr,na.rm = TRUE),0)



#summary function with both custom functions and inbuild functions to compare
sumStat <- function(x) {
  if(is.numeric(x)){
    c(min = min(x,na.rm = TRUE), my.min = my.min(x),
      max = max(x,na.rm = TRUE), my.max = my.max(x),
      mean = mean(x,na.rm = TRUE), my.mean = my.mean(x),
      median = median(x,na.rm = TRUE), my.medain = my.median(x),
      std = sd(x,na.rm = TRUE), my.sd = my.sd(x),
      var=var(x,na.rm = TRUE), my.var = my.var(x))
  }
  else{
    c(min = NA, my.min =NA,
      max = NA, my.max= NA,
      mean = NA, my.mean=NA,
      median = NA, my.medain = NA,
      std = NA, my.sd = NA,
      var=NA, my.var =NA) 
    
  }
}

#stastical summary of data
summary<- sapply(data, sumStat)
summary
#writing to csv file
write.csv(rbind(data,summary), file="summaryDefault.csv")

