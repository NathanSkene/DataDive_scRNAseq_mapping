convert_zeisel_tibble_to_exp_matrix <- function(tibbleIN,level=5){
  if(!level %in% 1:5){stop("level should be 1,2,3,4 or 5")}
  if(level==1){
    specificityA = data.frame(tibbleIN %>% select(Gene,Lvl1,Expr_sum_mean))
    specificityB = dcast(specificityA,Lvl1~Gene,fun.aggregate = mean)
  }else if(level==2){
    specificityA = data.frame(tibbleIN %>% select(Gene,Lvl2,Expr_sum_mean))
    specificityB = dcast(specificityA,Lvl2~Gene,fun.aggregate = mean)
  }else if(level==3){
    specificityA = data.frame(tibbleIN %>% select(Gene,Lvl3,Expr_sum_mean))
    specificityB = dcast(specificityA,Lvl3~Gene,fun.aggregate = mean)
  }else if(level==4){
    specificityA = data.frame(tibbleIN %>% select(Gene,Lvl4,Expr_sum_mean))
    specificityB = dcast(specificityA,Lvl4~Gene,fun.aggregate = mean)
  }else if(level==5){
    specificityA = data.frame(tibbleIN %>% select(Gene,Lvl5,Expr_sum_mean))
    specificityB = dcast(specificityA,Lvl5~Gene,fun.aggregate = mean)
  }
  colnames(specificityB)[1] = "ct"
  rownames(specificityB) = specificityB$ct
  specificityB = specificityB[,-1]
  specificityB = t(specificityB)
  return(specificityB)
}
