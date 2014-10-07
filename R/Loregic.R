edgelist2triplets<-function(grn){
  targets=unique(grn[,2]);triplets=data.frame(RF1=NULL,RF2=NULL,target=NULL)
  for (t in 1:length(targets)){
    x=grn[grn[,2]==targets[t],1];if(length(x)<2) next
    curr_ppt=data.frame(t(combn(x,2)),targets[t]);colnames(curr_ppt)=c('RF1','RF2','target')
    triplets=rbind(triplets,curr_ppt)}
return(triplets)}

loregic=function(grn,bdata){ 
  which.max.all=function(x,thld=-Inf){return(which(x==max(x)&x>thld))}
  if(ncol(grn)==2) motif=edgelist2triplets(grn) else motif=grn
  logic_gate_names=c('T=0','AND','RF1*~RF2','RF1','~RF1*RF2','RF2','XOR','OR','NOR','XNOR','~RF2','RF1+~RF2','~RF1','~RF1+RF2','NAND','T=1')
  int2bin=function(x, b=2){# convert integer vector x to base-b 
    xi=as.integer(x);if(any(is.na(xi) | ((x-xi)!=0))) print(list(ERROR="no integer input", x=x))
    ndigits=1+floor(logb(max(x), base=2));Base.b=array(NA, dim=c(length(x), ndigits))
    for(i in 1:ndigits){Base.b[, ndigits-i+1]=(x %% b);x=(x %/% b)}
    if(length(x)==1) Base.b[1, ] else Base.b}
  vecmatch=function(curr_data,tt2to1,flag_add=T){
    curr_element_counts=sapply(seq_len(nrow(curr_data)),function(k)which(apply(tt2to1,1,identical,as.numeric(curr_data[k,]))))
    n=floor(nrow(tt2to1)/2)
    if(flag_add){for(i in 1:n){if(!(2*i-1)%in%curr_element_counts && !(2*i)%in%curr_element_counts) curr_element_counts=c(curr_element_counts,2*i-1,2*i)}}
    curr_stats=rep(0,nrow(tt2to1));names(curr_stats)=1:nrow(tt2to1);curr_stats[as.numeric(names(table(curr_element_counts)))]=table(curr_element_counts) 
    return(curr_stats)}
  logfun_Zvalue=t(int2bin(x=0:15))
  logfun=vector("list",dim(logfun_Zvalue)[2]);names(logfun)=logic_gate_names
  for (i in 1:length(logfun)){logfun[[i]]=cbind(c(0,0,1,1),c(0,1,0,1),logfun_Zvalue[,i])}
  tt2to1=int2bin(x=0:7)
  logfun_level=sapply(1:length(logfun),function(i)sapply(1:4,function(j)which(apply(tt2to1,1,identical,logfun[[i]][j,]))))
  logfun_vects=matrix(0,8,16);for(i in 1:ncol(logfun_vects)){logfun_vects[logfun_level[,i],i]=1};colnames(logfun_vects)=logic_gate_names
  m=nrow(motif);motif_patt=matrix(0,m,length(logfun));colnames(motif_patt)=logic_gate_names
  for (i in 1:m){
    g1=intersect(motif[i,1],rownames(bdata));g2=intersect(motif[i,2],rownames(bdata));g3=intersect(motif[i,3],rownames(bdata))
    if(length(g1)==0||length(g2)==0||length(g3)==0) next else curr_data=na.omit(t(bdata[c(g1,g2,g3),]))
    if(nrow(curr_data)==0) next
    curr_stats=vecmatch(curr_data,tt2to1);curr_logfun_vects=1:ncol(logfun_vects)
    curr_succ_prob=rep(0,4);for(k in c(1,3,5,7)){
        curr_max=as.numeric(names(which.max.all(curr_stats[k:(k+1)])))
        mk=sum(curr_stats[k:(k+1)]);nk=curr_stats[curr_max[1]];curr_succ_prob[ceiling(k/2)]=(nk+1)/(mk+2)
        if(length(curr_max)==1)curr_logfun_vects=intersect(curr_logfun_vects,which(logfun_vects[curr_max,]==1))#curr_logfun_vects[,curr_logfun_vects[curr_max,]==1]
    }
    motif_patt[i,curr_logfun_vects]=prod(curr_succ_prob)
  }
  return(cbind(motif,motif_patt))
}