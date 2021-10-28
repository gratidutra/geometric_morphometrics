readmulti_nts<-function(filelist){   
  n<-length(filelist)
  names<-gsub(".nts","",filelist, ignore.case=T)
  landdata<-nind<-NULL
  for (i in 1:n){
    ntsfile<-scan(file=filelist[i],what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\"",quiet=TRUE)
    comment <- grep("\'", ntsfile)
    if (length(comment) != 0){
      ntsfile<-scan(file=file,what="char",quote="",sep="\n",strip.white=TRUE,comment.char="\'",quiet=TRUE)
    }
    header<-unlist(strsplit(ntsfile[1],"\\s+"))
    if(header[1]!=1){
      stop("NTS file not a rectangular matrix. First value in parameter line must be '1'.") }
    header<-casefold(header,upper=TRUE)
    dimval<-unlist(grep("DIM=",header))
    if(length(dimval)==0){
      stop("Header does not contain 'DIM=' designator.") }  
    labval<-unlist(grep("L",header))
    r.lab<-ifelse(is.element("2",labval)==TRUE,T,F)
    c.lab<-ifelse(is.element("3",labval)==TRUE,T,F)
    header<-sub("L","",header)
    header<-as.numeric(sub("DIM=","",header))
    missdata<-ifelse(header[4]!=0,T,F)
    if(missdata==TRUE){missval<-ifelse(dimval==6,header[5],header[6]) } 
    p<-header[2];k<-header[3]
    nind<-rbind(nind,p)
    if (min(nind)!=max(nind)) {
      stop("Number of landmarks not the same in all files.") } 
    tmp<-unlist(strsplit(ntsfile[-1],"\\s+"))
    rowlab<-NULL; 
    if(r.lab==TRUE){
      rowlab<-tmp[1:p]
      tmp<-tmp[-(1:length(rowlab))]   }
    if(c.lab==TRUE){ tmp<-tmp[-(1:k)] }
    if(missdata==TRUE){tmp[grep(missval,as.integer(tmp))] <- NA}
    options(warn=-1)
    data<-matrix(as.numeric(tmp),ncol=k,byrow=TRUE)
    landdata<-rbind(landdata,data)
  }
  coords<-arrayspecs(landdata,p,k)
  if(sum(which(is.na(landdata)==TRUE))>0){cat("NOTE.  Missing data identified.")}
  dimnames(coords)[[3]]<- as.list(names)
  return(coords=coords)
}
