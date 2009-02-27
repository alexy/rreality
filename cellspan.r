qhr <- function(t) floor(as.numeric(day.dt.secs(t)/900))

v <- sqlFetch(chan,"cellspan")
v1$qhr <- qhr(v1$starttime)
v1q <- split(v1,v1$qhr)
names(v1q)==0:99 # we exceed 96 because of EDT on 10/31->11/1
sapply(v1q,dim)[1,]

# equivalent methods to get the cell:
length(v1q[["0"]]$cell)
length(q1$"0"$cell)
length(q1$`0`$cell)

# find dominating towers
c1 <- sapply(q1, function(x) sort(table(x$cell),decr=T))

# simplifying names
names(v)[which(names(v)=="person_oid")]    <- "person"
names(v)[which(names(v)=="celltower_oid")] <- "cell"

# preparing ddply function:
names(sort(table(v1[v1$hr==1,]$cell),decr=T))[1]

# discretize hour and quarter-hour chunks
# NB: we get 25 hours and 100 chunks max 
# -- due to EDT 10/31-11/1
v$hr <- hr(v$starttime)
v$qhr <- qhr(v$starttime)

# get durations and their outliers
v$delta <- v$endtime - v$starttime
sort(v$delta,decr=T)[1:100]

vh <- ddply(v, .(person,hr), function (.dh) names(sort(table(.dh$cell),decr=T))[1])
names(vh)[3] <- "cell"

vh <- vh[order(vh$hr,vh$cell),]
vhc <- ddply(vh,.(hr,cell),function(.df) length(.df$cell))
names(vhc)[3] <- "count"

# for each hour, sort by deacreasing count
vhc <- vhc[order(vhc$hr,-vhc$count),]
# take top 5 counts from that:
vh5 <- ddply(vhc,.(hr),function(.df) .df[1:5,])

# most frequent cells overall, and the rarest
cell <- sort(v$cell)
cell.rle <- rle(cell)
crle <- data.frame(lengths=cell.rle$lengths,values=cell.rle$values)
crle <- crle[order(-crle$lengths),]
l <- crle$lengths
sum(l>1)
sum(l>2)
barplot(l[l>10000])

# lag() doesn't really change a timeseries, embed() does:
sum(embed(s1,2)[,1] != head(e1,-1))
which(head(e1,-1) != embed(s1,2)[,1])

# check ordering
v1so <- order(v1$starttime)
which(v1so != 1:length(v1so))
# integer(0)
v1eo <- order(v1$endtime)
which(v1eo != 1:length(v1eo))
# 111 elements

v1 <- v[v$person==1,]
v1cu <- sort(unique(v1$cell))
v1m <- matrix(0,nrow=length(v1cu),ncol=length(v1cu))
cellmap <- 1:length(v1cu)

namemap <- function(a,x) a[as.character(x)]
names(cellmap) <- v1cu

Reduce(function(from,to) {to <- namemap(cellmap,to); v1m[from,to] <<- v1m[from,to]+1; to}, tail(v1$cell,-1), init=namemap(cellmap,v1$cell[1]))

which(v1m==max(v1m))

# scalar -- single argument
# matcoords <- function(n,nrow) { row <- n %% nrow; col <- n %/% nrow; c(row,col) }

# vector -- multiple arguments, e.g. as might be returned by which
matcoords <- function(x,nrow) { sapply(x, function(n) { row <- n %% nrow; col <- n %/% nrow + 1; c(row,col) }) }

l1 <- v1m
# how can we pass * to Reduce shorter?
dim(l1) <- Reduce(function(x,y) x*y,dim(v1m))
names(l1) <- 1:length(l1)
l1 <- sort(l1[],decr=T)
sum(l1[l1>0]) # ==length(v1$cell) !
sum(l1>0) # smaller!

# from vh5, looks like office
# office <- Reduce(function(acc,hr) {append(acc,vh5[vh5$hr==hr,]$cell}, 11:17,init=c()) # strange results: 11 => 1 38 89 112 69
# office <- Reduce(function(acc,hr) cat(hr,"=>",vh5[vh5$hr==hr,]$cell,"\n"), 11:17, init=c())

office <- unique(vh5[vh5$hr %in% 11:17,]$cell)

coordsmat <- function(row,col,nrow) (col-1)*nrow+row

# NB is it a reasonable hash in R?
# NB how to make xlist a parameter,
# so that it persists in parent environment?
slist <- function(row,col,oid) {
	label <- paste(row,col)
	if (is.null(xlist[[label]])) {
		xlist[[label]] <<- c(oid)
	} else { 
		xlist[[label]] <<- append(xlist[[label]],oid)
	}
}

v1cell <- v1$cell
names(v1cell) <- rownames(v1)
rowlinks  <- rownames(v1)

v1m <- matrix(0,nrow=length(v1cu),ncol=length(v1cu))
xlist <- list() # for slist

Reduce(function(from.i,to) {
		from <- from.i[1]; i <- from.i[2]
		# backlink <- names(to) ## empty -- dropped?
		# backlink <- rownames(v1)[i] ## slooow!
		backlink <- rowlinks[i]
		to <- namemap(cellmap,to) 
		v1m[from,to] <<- v1m[from,to]+1
		slist(from,to,backlink)
		c(to,i+1)
	}, 
	tail(v1cell,-1), 
	init=c(namemap(cellmap,v1cell[1]),1)
	)
	
# NB manually verify, for
#> cellmap[496]
#15375 
#  496 
#> cellmap[490]
#15369 
#  490 
# that we have 1834 pairs (15375, 15369) in v1$cell

p1 <- xlist[[1]]
v1[as.character(as.numeric(p1)-1),]
