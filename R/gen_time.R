genTimeFn<-function(object){
  
  n=apply(exp(-m(object)),c(2,6),cumprod)
  
  apply((stock.wt(object)%*%mat(object)%*%n)%*%ages(stock.wt(object)),c(2,6),sum)%/%
    apply(stock.wt(object)%*%mat(object)%*%n,c(2,6),sum)}
