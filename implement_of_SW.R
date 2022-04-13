#!/usr/bin/env Rscript 

### Usage: Rscript --vanilla hw1.R <input file> <score file> <path_of_output>
### Example: Rscript --vanilla hw1.R input.txt blosum62.txt /usr/yuhang/course/cbb752/hw/output.txt
### Note: Smith-Waterman Algorithm 


args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) { 
  stop("At least two arguments must be supplied (inputFile, scoreFile).n", call.= FALSE) 
} else if (length(args)>=3) { 
  # default gap penalties    
  args[4] = -2    
  args[5] = -1    
} 

## Specifying author and email 
p <- c(person("Yuhang", "Chen", role = "aut", 
              email = "yuhang.chen@yale.edu")) 

## define a function to get letter of sequence by index
get_letter <-  function(word, i){
    a <- substring(word,i,i)
    return(a)
}

## define a function to enumerate letters of a given string
enumerate_letter <-  function(word){
  tmp <- rep("?",length(word))
  for (i in 1:nchar(word)) {
    tmp[i] <- substring(word,i,i)
  }
  return(tmp)
}

## define a function, when there are more than one maximum, only return the index of first one
max_index <- function(vector){
  tmp <- which(vector==max(vector))
  if (length(tmp)>1) {
    return(tmp[1])
  } else{
    return(tmp)
  }
}

## Implement Smith-Waterman Algorithm 
runSW <- function(inputFile, scoreFile, openGap, extGap, outputPath) {
  openGap <- as.numeric(openGap)
  extGap <- as.numeric(extGap)
  
  ## read in the files
  suppressWarnings(input <- read.table(inputFile,header = F))
  suppressWarnings(sim_mat <- read.table(scoreFile,header=T))
                      
  query <- input$V1[1] ## in default, query sequence is the first one,
  target <- input$V1[2]  ## target sequence is the second one.
  
  query_char <- enumerate_letter(query)
  target_char <- enumerate_letter(target)
  
  #------------ calculation -------------#
  
  ## initialization: Creating a score matrix of (n+1)*(m+1), 
  ## here n is the length of query and m is the length of target sequence
  m <- nchar(query)
  n <- nchar(target)
  score_mat <- matrix(data = 0,nrow = n+1, ncol = m+1)
  
  ## source_mat is the matrix to record the source of value of each cell
  ## 0: the scores of all 3 neighbors are negative, curren cell takes 0 
  ## 1: take the upper-diag value
  ## 2: take the left value
  ## 3: take the top value
  source_mat <- matrix(data = 0,nrow = n+1, ncol = m+1)
  
  ## filling the matrix
  ## first filling the first row and the first column with 0
  for (i in 2:(n+1)) {
    for (j in 2:(m+1)) {
      
      ## then consider the value of current cell, which depends on the 
      ## match/mismatch of the same position of the sequences
      target_letter <- get_letter(target, i-1)
      query_letter <- get_letter(query, j-1)
      ## derive the similarity value according the letters
      current_score <- sim_mat[query_letter, target_letter] 
      
      ## consider about the value of the neighbor cells of current cell
      left <- score_mat[i,j-1]
      top <- score_mat[i-1,j]
      upper_diag <- score_mat[i-1,j-1]
      
      ## enumerate all the possibilities of introducing gaps
      left_nbs <- rep(-99, j-1)
      for (k in 1:(j-1)) {  ## consider all the left neighbors for current cell in the same row
        left_nbs[k] <- score_mat[i,k] + openGap + extGap * (j-k-1)
      }
      left_indel <- max(left_nbs)
      
      top_nbs <- rep(-99, i-1)
      for (s in 1:(i-1)) {  ## consider all the top neighbors for current cell in the same column
        top_nbs[s] <- score_mat[s,j] + openGap + extGap * (i-s-1)
      }
      top_indel <- max(top_nbs)
      
      max_score <- max(upper_diag+current_score, left_indel, top_indel)
      score_mat[i,j] <- ifelse(max_score <= 0, 0, max_score)
      
      ## write the source of score for each cell
      source_token <- max_index(c(upper_diag+current_score, left_indel, top_indel))
      source_mat[i,j] <- ifelse(max_score <= 0, 0, source_token)
      
    }
  }
  
  colnames(score_mat) <- c(" ",query_char)
  rownames(score_mat) <- c(" ",target_char)
  colnames(source_mat) <- c(" ",query_char)
  rownames(source_mat) <- c(" ",target_char)
  
  ## tracing back
  ## start from the cell which has the highest score. 
  ## Then the next cell moving to could be the left, the top, and the upper diagonal cells. 
  ## If it is a match, then go diagonally; 
  ## if it is a mismatch, then go to the neighbor according to the source matrix
  highest <- max(score_mat)
  highest_idx <- which(score_mat == highest,arr.ind = T)
  row_idx <- highest_idx[1]
  col_idx <- highest_idx[2]
  
  q_letter <- colnames(score_mat)[col_idx] 
  t_letter <- rownames(score_mat)[row_idx]
  
  ## to store the matching results
  query_matched <- ""
  target_matched <- ""
  
  while (score_mat[row_idx,col_idx] >= 1) {  ## stop tracing back when current cell equals 0
    
    ## if it is a match (by match, here means the 2 letters are the same), go upper diagonally
    if (q_letter == t_letter) {
      new_q <- colnames(score_mat)[col_idx]
      new_t <- rownames(score_mat)[row_idx]
      row_idx <- row_idx - 1
      col_idx <- col_idx - 1
    } else if (q_letter != t_letter) { ## if it is a mismatch, check the source of value of the current cell
      
      if (source_mat[row_idx,col_idx] == 1) { # score is from upper-diag neighbor
        new_q <- colnames(score_mat)[col_idx]
        new_t <- rownames(score_mat)[row_idx]
        row_idx <- row_idx - 1
        col_idx <- col_idx - 1
      } else if (source_mat[row_idx,col_idx] == 2) { # score is from left neighbor
        new_q <- colnames(score_mat)[col_idx]
        new_t <- "-"
        row_idx <- row_idx
        col_idx <- col_idx - 1
      } else if (source_mat[row_idx,col_idx] == 3) { # score is from the top neighbor
        new_q <- "-"
        new_t <- rownames(score_mat)[row_idx]
        row_idx <- row_idx - 1
        col_idx <- col_idx
      }
    }
    
    query_matched <- c(new_q,query_matched)
    target_matched <- c(new_t,target_matched)
    q_letter <- colnames(score_mat)[col_idx]
    t_letter <- rownames(score_mat)[row_idx]
  }
  
  q_matched <- paste(query_matched, collapse="")
  t_matched <- paste(target_matched, collapse="")
  
  #------------ write output ------------#
  #path <- "/Users/yuhangchen/Downloads/HW1_cbb752b22_programming_supp_files/test.out.txt"
  path <- outputPath
  
  ## write query and target sequences
  dash1 <- "-----------"
  cat(dash1,file = path)
  cat("\n",file = path, append = T)
  seq <- "|Sequences|"
  cat(seq,file = path, append = T)
  cat("\n",file = path, append = T)
  cat(dash1,file = path,append = T)
  cat("\n",file = path, append = T)
  seq1 <- "sequence1"
  cat(seq1,file = path, append = T)
  cat("\n",file = path, append = T)
  cat(query,file = path, append = T)
  cat("\n",file = path, append = T)
  seq2 <- "sequence2"
  cat(seq2,file = path, append = T)
  cat("\n",file = path, append = T)
  cat(target,file = path, append = T)
  cat("\n",file = path, append = T)

  ## write score matrix  
  dash2 <- "--------------"
  cat(dash2,file = path,append = T)
  cat("\n",file = path, append = T)
  sm <- "|Score Matrix|"
  cat(sm,file = path, append = T)
  cat("\n",file = path, append = T)
  cat(dash2,file = path,append = T)
  cat("\n",file = path, append = T)
  # print colnames
  cat("\t",file = path, append = T)
  cat("\t",file = path, append = T)
  for (i in 1:nchar(query)) {
    cat(query_char[i],file = path, append = T)
    cat("\t",file = path, append = T)
  }
  cat("\n",file = path, append = T)
  
  # print the row with all 0
  cat("\t",file = path, append = T)
  for (i in 1:(nchar(query)+1)) {
    cat(0,file = path, append = T)
    cat("\t",file = path, append = T)
  }
  cat("\n",file = path, append = T)
  
  # print other rows
  for (i in 1:nchar(target)) {
    cat(target_char[i],file = path, append = T)
    cat("\t",file = path, append = T)
    for (j in 1:(nchar(query)+1)) {
      cat(score_mat[i+1,j],file = path, append = T)
      cat("\t",file = path, append = T)
    }
    cat("\n",file = path, append = T)
  }
  
  ## write alignment info
  dash3 <- "----------------------"
  cat(dash3,file = path,append = T)
  cat("\n",file = path, append = T)
  bla <- "|Best Local Alignment|"
  cat(bla,file = path,append = T)
  cat("\n",file = path, append = T)
  cat(dash3,file = path,append = T)
  cat("\n",file = path, append = T)
  cat("Alignment Score:",file = path,append = T)
  cat(highest,file = path,append = T)
  cat("\n",file = path, append = T)
  cat("Alignment Results:",file = path,append = T)
  cat("\n",file = path, append = T)
  
  ## consider about '|'
  qq <- enumerate_letter(q_matched)
  tt <- enumerate_letter(t_matched)
  dat <- as.data.frame(rbind(qq,tt))
  dat <- t(dat)
  record <- rep(-1,nrow(dat))
  for (r in 1:nrow(dat)) {
    ## for a position, if both of the letters are not "-", then using 1 representing a '|'
    record[r] <- ifelse(dat[r,"qq"]==dat[r,"tt"], 1,0)
  }
  lines <- rep(".",nrow(dat))
  for (t in 1:nrow(dat)) {
    lines[t] <- ifelse(record[t]==1,"|"," ")
  }
  lines <- paste(lines,collapse = "")
  
  ## consider about the spaces on the both sides of unalignment
  max_idx <- max(col_idx,row_idx)
  # to see which one has the spaces on the left; 
  # lspaces: 1 -> target seq has spaces on the left;
  # lspaces: 2 -> query seq has spaces on the left
  lspaces <- max_index(c(col_idx,row_idx))
  
  if (lspaces == 1) { # target seq has spaces on the left
    
    ## write the unmatched letters on the left of query
    left_char_query <- substring(query,1,col_idx-1)
    cat(left_char_query,file = path,append = T)
    
    ## write the matched letters of query
    cat(paste("(",q_matched,")",sep = ""),file = path,append = T)
    
    ## write the unmatched letters on the right of query
    right_char_query <- substring(query,highest_idx[2],nchar(query))
    cat(right_char_query,file = path,append = T)
    cat("\n",file = path, append = T)
    
    ## write the lines
    left_lines <- rep(" ",max_idx)
    left_lines <- paste(left_lines,collapse = "")
    cat(left_lines,file = path,append = T)
    cat(lines,file = path,append = T)
    right_lines <- rep(" ", nchar(query) - highest_idx[2]+2)
    right_lines <- paste(right_lines,collapse = "")
    cat(right_lines,file = path,append = T)
    cat("\n",file = path, append = T)
    
    ## write the spaces on the left of target
    left_space_target <- rep(" ",abs(row_idx-col_idx))
    left_space_target <- paste(left_space_target,collapse = "")
    cat(left_space_target,file = path,append = T)
    
    ## write the unmatched letters on the left of target
    left_char_target <- substring(target,1,row_idx-1)
    cat(left_char_target,file = path,append = T)
    
    ## write the matched letters of target
    cat(paste("(",t_matched,")",sep = ""),file = path,append = T)
    
    ## write the unmatched letters on the right of target
    right_char_target <- substring(target,highest_idx[1],nchar(target))
    cat(right_char_target,file = path,append = T)
    
    ## write the spaces on the right of target
    right_space_target <- rep(" ", nchar(query) - highest_idx[2]+1)
    right_space_target <- paste(right_space_target,collapse = "")
    cat(right_space_target,file = path,append = T)
    cat("\n",file = path, append = T)
    
  } else if (lspaces == 2) { # query seq has spaces on the left
    ## write the spaces on the left of query
    left_space_query <- rep(" ",abs(row_idx-col_idx))
    left_space_query <- paste(left_space_query,collapse = "")
    cat(left_space_query,file = path,append = T)
    
    ## write the unmatched letters on the left of query
    left_char_query <- substring(query,1,col_idx-1)
    cat(left_char_query,file = path,append = T)
    
    ## writh the matched letters of query
    cat(paste("(",q_matched,")",sep = ""),file = path,append = T)
    
    ## write the unmatched letters on the right of query
    right_char_query <- substring(query,highest_idx[2],nchar(query))
    cat(right_char_query,file = path,append = T)
    
    ## write the spaces on the right of query
    right_space_query <- rep(" ", nchar(target) - highest_idx[1]+1)
    right_space_query <- paste(right_space_query,collapse = "")
    cat(right_space_query,file = path,append = T)
    cat("\n",file = path, append = T)
    
    ## write the lines
    left_lines <- rep(" ",max_idx)
    left_lines <- paste(left_lines,collapse = "")
    cat(left_lines,file = path,append = T)
    cat(lines,file = path,append = T)
    right_lines <- rep(" ", nchar(query) - highest_idx[2]+2)
    right_lines <- paste(right_lines,collapse = "")
    cat(right_lines,file = path,append = T)
    cat("\n",file = path, append = T)
    
    ## write the unmatched letters on the left of target
    left_char_target <- substring(target,1,row_idx-1)
    cat(left_char_target,file = path,append = T)
    
    ## write the matched letters of target
    cat(paste("(",t_matched,")",sep = ""),file = path,append = T)
    
    ## write the unmatched letters on the right of target
    right_char_target <- substring(target,highest_idx[1],nchar(target))
    cat(right_char_target,file = path,append = T)
    cat("\n",file = path, append = T)
  } 

} 

## Run the main function and generate results 
#input <- "/Users/yuhangchen/Downloads/HW1_cbb752b22_programming_supp_files/input.txt"
#sim_mat <- "/Users/yuhangchen/Downloads/HW1_cbb752b22_programming_supp_files/blosum62.txt"
#path <- "/Users/yuhangchen/Downloads/HW1_cbb752b22_programming_supp_files/real.output.txt"

#runSW(inputFile = input,scoreFile = sim_mat,openGap = -2,extGap = -1,outputPath = path)

suppressWarnings(runSW(inputFile=args[1], scoreFile=args[2], openGap=args[4], 
                       extGap=args[5],outputPath = args[3]))

