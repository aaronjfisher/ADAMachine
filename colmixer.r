colmixer<-function(s){#sequence of values (numeric)
  mag<-max(s)-min(s)
	scaled<-(s-min(s))*4/(mag)
	red<-c(); green<-c(); blue<-c() #will hold weights for the rgb function
	for (j in 1:length(s)){
		if (scaled[j]<=1) {red[j]<-1; green[j]<-0;blue[j]<-scaled[j]}
		else if (scaled[j]<=2) {red[j] <- 1-(scaled[j]-1); green[j]<-0; blue[j]<-1} #red goes down
		else if (scaled[j]<=3) {red[j] <- 0; green[j]<- .85*(scaled[j]-2); blue[j]<-1} #green goes up
		else if (scaled[j]>3) {red[j] <- 0; green[j]<- .85; blue[j]<- 1-(scaled[j]-3)} #blue goes down
	}
	return(rgb(red,green,blue))
}
