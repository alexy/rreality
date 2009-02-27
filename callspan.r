#!/usr/bin/env r
#library(RODBC)
#chan <- odbcConnect("mit_reality")
#sqlTables(chan)
#sqlColumns(chan,"callspan")


# remove outliers
# z <- z[z$duration>=0 & z$duration <= 11000,]

# remove empty columns
# which(!is.na(z$number))
# which(!is.na(z$remote))
# z$number <- NULL
# z$remote <- NULL

# sort by person_oid, then starttime within each person:
# z <- z[order(z$person_oid,z$starttime),]

# test it's sorted within each person_oid,
# could use table instead of rle:
load("~/R/reality/callspanz.rda") # z
print(dim(z))

sum(sapply(rle(z$person_oid)$values,function(i){x <- z[z$person_oid==i,]$starttime; sum(x!=sort(x))}))

# compute inter-call intervals
interstart <- c(NA,diff(z$starttime))
z.p.rle <- rle(z$person_oid)
# head(x,-1) removes the last element:
z.p.first <- head(c(1,cumsum(z.p.rle$lengths)+1),-1)
interstart[z.p.first] <- NA
z$interstart <- interstart
sum(z$interstart!=0,na.rm=T) # almost all zeroes?! duplicates!

# remove duplicate rows indicated by interstart==0
z.meat <- c("starttime","endtime","duration","person_oid","callid","description","direction","status")
zi0 <- z[which(z$interstart==0),z.meat]
zi1 <- z[which(z$interstart==0)-1,z.meat]
zir <- apply(zi0!=zi1,1,sum)
z.duprownames <- names(which(zir==0))
z.keep <- setdiff(rownames(z),z.duprownames)
zorig <- z
z <- z[z.keep,]

day.dt.secs <- function(t) difftime(t,trunc(t,units="days"),units="secs")

posix.day  <- function(t) with(as.POSIXlt(t), yday)
posix.week <- function(t) floor(posix.day(t)/7)
posix.wday <- function(t) with(as.POSIXlt(t), wday)

z$daystart <- day.dt.secs(z$starttime)
z$week <- posix.week(z$starttime)

zw <- split(z,z$week)

nz <- function(v) v[v>0]
nanan <- function(x) is.na(x) | is.nan(x)

not.nanan <- function(x,fun) if (length(x[!nanan(x)])==0) 0 else { res <- do.call(fun,list(x[!nanan(x)],na.rm=T)); if (nanan(res)) 0 else res }

tbl.to.df <- function(tbl) {
	df <- as.data.frame(matrix(tbl,nrow=1))
	colnames(df) <- names(tbl)
	df
}

mean.da   <- function(x) not.nanan(x,mean)
median.da <- function(x) not.nanan(x,median)
sd.da     <- function(x) not.nanan(x,sd)
mad.da    <- function(x) not.nanan(as.numeric(x),mad)

# is.weekend <- function(x) x in 1:

cool.customer <- function(.df) { 
	call.type      <- tbl.to.df(table(.df$description))
	colnames(call.type) <- c("Packet","Sms","Voice")

	call.direction <- tbl.to.df(table(.df$direction))
	colnames(call.direction)[2] <- "Missed"
	
	call.status    <- tbl.to.df(table(.df$status))
	colnames(call.status)[c(1,6)] <- c("Ok","Nodelivery")
	
	some <- data.frame(
		n=dim(.df)[1],
		nz=sum(.df$duration>0,na.rm=T), # sum(nz(duration))
		tcl=sum(.df$duration,na.rm=T),
		
		#nz.bd=sum(.df[.df$wday]$duration>0,na.rm=T),
		
		duration.mean=mean.da(.df$duration),
		duration.median=median.da(.df$duration),
		duration.sd=sd.da(.df$duration),
		duration.mad=mad.da(.df$duration),
		
		# found experimentally NaN entries, throw off svm --
		# so we remove them here -- better mean? or median?
		
		nz.duration.mean=mean.da(nz(.df$duration)), 		nz.duration.median=median.da(nz(.df$duration)),
		nz.duration.sd=sd.da(nz(.df$duration)),
		nz.duration.mad=mad.da(nz(.df$duration)),
		
		start.mean=mean.da(day.dt.secs(.df$starttime)),
		start.median=median.da(day.dt.secs(.df$starttime)),
		start.sd=sd.da(day.dt.secs(.df$starttime)),
		start.mad=mad.da(day.dt.secs(.df$starttime)),

		# remove NaN -- how come they show up?
		interstart.mean=mean.da(.df$interstart),
		interstart.median=median.da(.df$interstart),
		interstart.sd=sd.da(.df$interstart),
		interstart.mad=mad.da(.df$interstart)
	)	
	cbind(some, call.type, call.direction, call.status)
}

library(plyr)
#zw.cc <- lapply(zw,function(.df) ddply(.df,.(person_oid),cool.customer))
#save(zw.cc,file="~/R/Reality/zw.cc.rda")
load("~/R/Reality/zw.cc.rda")

Reduce2 <- function(f,L) { res <- NULL; for(e in L) if (is.null(res)) res <- e else res <- f(res,e) }

zr <- Reduce(rbind,zw.cc)

# classification! svm likes factors
zr$person_oid <- as.factor(zr$person_oid)

# svm cannot rescale Failed and Pending on zr[1:100,]:
zrr <- subset(zr,select=-c(Failed,Pending))

library(e1071)
m1  <- svm(person_oid~.,data=zrr[1:100,])
m2  <- svm(person_oid~.,data=zrr[1:200,])
m3  <- svm(person_oid~.,data=zrr[1:300,])
m4  <- svm(person_oid~.,data=zrr[1:400,])
# m11 <- svm(person_oid~.,data=zrr[1:1100,])
# m12 <- svm(person_oid~.,data=zrr[1:1200,])
# m13 <- svm(person_oid~.,data=zrr[1:1300,])
# m14 <- svm(person_oid~.,data=zrr[1:1400,])
# 
# 
# library(party)
# ct1  <- ctree(person_oid~.,data=zrr[1:100,])
# ct5  <- ctree(person_oid~.,data=zrr[1:500,])
# ct11 <- ctree(person_oid~.,data=zrr[1:1100,])

range <- 1:100
z100 <- zrr$person_oid[range]
p100 <- predict(m2,zrr[range,])

inc <- function(x,by=1) eval.parent(x <- x+by)

# this demonstrates indexing by factors
# alternative is to use as.numeric right away
# and positions directly instead of as.character names
# NB how do we pass pips by reference?
inc.pips.facts <- function(pips,a,b) {
	comp.fact <- a[a==b]
	# is as.numeric applied directly to factors,
	# they collapse into contiguous integer range!
	comp <- as.numeric(as.character(comp.fact))
	comp.rle <- rle(sort(comp))
	#inc(pips[as.character(comp.rle$values)],comp.rle$lengths)
	pips[as.character(comp.rle$values)] <- 
		pips[as.character(comp.rle$values)] <- comp.rle$lengths
	pips
}

inc.pips.ints <- function(pips,a,b) {
	comp.fact <- a[a==b]
	# shrinks to contiguous range:	
	comp <- as.numeric(comp.fact)
	comp.rle <- rle(sort(comp))
	#inc(pips[comp.rle$values],comp.rle$lengths)
	pips[comp.rle$values] <- 
		pips[comp.rle$values]+comp.rle$lengths
	pips
}

# NB doesn't work! -- see inc.pips.{ints,facts}
# inc.pips <- function(pips,a,b) {
#   comp.fact <- a[a==b]
#   for (i in comp.fact) pips[i] <- pips[i]+1
# }


match <- function(model,df,step=100) {
	maxblock <- dim(df)[1] %/% step - 1
	levels <- levels(df$person_oid)
	pips <- rep(0,length(levels))
	names(pips) <- levels
	
	dfx <- subset(df,select=-person_oid)

	hits <- rep(0,maxblock)

	for (i in 1:maxblock) {
		base <- (i-1)*step
		range <- (base+1):(base+100) 
		a <- df$person_oid[range]
		b <- predict(model,dfx[range,])
		pips <- inc.pips.ints(pips,a,b)
		hits[i] <- sum(a==b)
	}
	list(hits=hits,pips=pips)
}


# one way to accumulate hits per person
zrr.p.levels <- levels(zrr$person_oid)

maxrow <- 1500
step <- 100

match.m <- match(m4,zrr,step)

real <- zrr$person_oid[1:(maxrow+step)]
real.rle <- rle(sort(as.numeric(real)))

precs <- sort(match.m$pips/real.rle$lengths,decreasing=T)
precs.non0.len <- length(precs[precs>0])
precs.non0 <- precs[1:(precs.non0.len+1)]
print(precs)
print(precs.non0.len)
postscript("mit-reality-svm-m2.eps",horizontal=F,onefile=F,height=3.5,width=3.2)
plot(precs.non0,ann=F)
title(main="SVM on Call Timings")
title(xlab="person ID", ylab="fraction correct")
dev.off()


stop("enough for now")


# now let's see how precision corresponds to the
# amount of data

zp <- split(z,z$person_oid)

zp.length <- sapply(zp,dim)[1,]
zp.length.order <- order(zp.length,decreasing=T)
zp.length.rank <- 1:length(zp.length)
names(zp.length.rank) <- names(zp.length)[zp.length.order]
print(zp.length.rank[names(prec)][1:41])


# adding weekday to zw
for (i in seq_along(zw)) zw[[i]]$wday <- posix.wday(zw[[i]]$starttime)

# weekend dataframe separate
zw.we <- lapply(zw, function(.df) ddply(.df[.df$wday %in% c(6,0),],.(person_oid),cool.customer))
zw.bd <- lapply(zw, function(.df) ddply(.df[.df$wday %in% 1:5,],.(person_oid),cool.customer))

zrw <- Reduce(rbind, zw.we)
zrw$person_oid <- as.factor(zrw$person_oid)
zrwr <- subset(zrw,select=-c(Failed,Pending,Nodelivery))

m1w <- svm(person_oid~.,data=zrwr[1:100,])
m2w <- svm(person_oid~.,data=zrwr[1:200,])

match.m1w <- match(m1w,zrwr,maxrow,step)
match.m2w <- match(m2w,zrwr,maxrow,step)

print("match.m1w:",match.m1b)
print("match.m2w:",match.m2b)

zrb <- Reduce(rbind, zw.bd)
zrb$person_oid <- as.factor(zrb$person_oid)
zrbr <- subset(zrb,select=-c(Failed,Pending))

m1b <- svm(person_oid~.,data=zrbr[1:100,])
m2b <- svm(person_oid~.,data=zrbr[1:200,])

match.m1b <- match(m1b,zrbr,maxrow,step)
match.m2b <- match(m2b,zrbr,maxrow,step)

print("match.m1b:",match.m1b)
print("match.m2b:",match.m2b)

# kmeans clustering on time averages
zrr.mat <- zrr
zrr.mat$start.mean <- as.numeric(zrr.mat$start.mean)
zrr.mat$start.median <- as.numeric(zrr.mat$start.median)
zrr.km10 <- kmeans(zrr.mat,10)
zrr.km20 <- kmeans(zrr.mat,20,iter=100) # didn't converge in 10
zrr.km10.pc <- data.frame(p=zrr$person_oid,c=zrr.km10$cluster)
zrr.km10.pc <- with(zrr.km10.pc, zrr.km10.pc[order(p,c),])
