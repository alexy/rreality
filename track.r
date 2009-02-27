load("callspan-v.rda")
vt <- with(v,v[order(starttime,cell,person),])
vt.s <- as.integer(vt$start)

sum(vt.s!=sort(vt.s)) # check vt is really sorted by starttime
d <- subset(vt,select=c(starttime,cell,person))
names(d)[1] <- "start"

julian.posix <- function(pt) as.date(as.character(trunc(pt,units="days")),order="ymd")
julian.day   <- function(pt) as.numeric(julian.posix(pt))

height <- function(.df)dim(.df)[1]

right <- function(d,delta) {
  dmax <- height(d)
  left <- rep(0,dmax) 
  right <- rep(0,dmax)
  
  for (i in 1:dmax) {
    j <- 1
    from  <- d$start[i]
    while (((i+j) <= dmax) && ((d$start[i+j] - from) < delta)) j <- j + 1
    j <- j - 1
    right[i] <- j
  }
  
  for (i in dmax:1) {
    j <- 1
    from <- d$start[i]
    while (((i-j) >= 1) && ((from - d$start[i-j]) < delta)) j <- j + 1
    j <- j - 1
    left[i] <- j
  }
  list(left=left,right=right)
}


# peers1 doesn't record the original person in h
# -- he can be looked up in d[i,] for h[[i]] --
# -- and also doesn't record row names for q's
peers1 <- function(d,x) {
  dmax <- height(d)
  h <- list()
  for (i in 1:dmax) {
    p <- as.character(d$person[i])
    c <- d$cell[i]
    g <- list()
    for (j in (i-x$left[i]):(i+x$right[i])) {
      q <- as.character(d$person[j])
      if (d$cell[j]==c && q!=p) # 3=="3" => TRUE!
      if (is.null(g[[q]])) g[[q]] <- 1 else g[[q]] <- g[[q]]+1
    }
    if (length(g)>0) h[[i]] <- g
    if (i %% 1000 == 0) print(i)
  }
  h
}


peers2 <- function(d,x) {
  # as.numeric below caused 9e+05 as.character later,
  # leading to NA when indexing by such a rowname
  rnames <- as.integer(rownames(d))
  dmax <- height(d)
  h <- list()
  for (i in 1:dmax) {
    p <- d$person[i]
    c <- d$cell[i]
    g <- data.frame()
    for (j in (i-x$left[i]):(i+x$right[i])) {
      q <- d$person[j]
      if (d$cell[j]==c && q!=p) {
        j.row <- rnames[j]
        g <- rbind(g,c(q,j.row))
      }
    }
    if (length(g)>0) {
      i.row <- rnames[i]
      g <- rbind(c(p,i.row),g)
      names(g) <- c("person","row")
      # for OCaml, we won't reorder;
      # we don't really use order in x.pairs anyways,
      # restoring the p from the original d;
      # if we keep it first, it's faster
      h[[i]] <- with(g,g[order(person,row),])
    }
    if (i %% 1000 == 0) print(i)
  }
  h
}


x30 <- right(d,30)

# peers1
# ------
#ph <- peers(d,x30) # 2 hours! 
#load("track-ph.rda")
pl <- sapply(ph,length)
max(pl) # 8
which(pl==max(pl))
names(pl) <- 1:length(pl)
pl2 <- pl[pl==2]
pl1 <- pl[pl==1]
pl1[1:10]
pl1i <- as.integer(names(pl1))
pl1p <- d$person[pl1i] # original person in the pair
pl1q <- as.integer(sapply(ph[pl1i],names)) # second person on the pair

# we may drop the selfsame person from each person's list,
# but this is not what drop=T does here -- it eliminates empty factors:
dp <- split(d,d$person) 
dpr <- lapply(dp,rownames)
# we can always get a position from character contents via which;
# so we don't get much by naming with positions:
# dprn <- lapply(dpr,function(li) { len <- length(li); names(li) <- 1:len; li })
# -- we'll name with contents instead, making position the new content:
dprn <- lapply(dpr,function(li) { len <- length(li); il <- 1:len; names(il) <- li; il })

dr <- rownames(d)
drn <- as.integer(dr)
drnn <- 1:length(dr)
names(drnn) <- dr

# for peers2, let's see if we ever get a neighbor more than once,
# and need to concatenate instead of assigning...
# perhaps here 2M+ -Inf's and unlist slow things down to ~5 min:
# NB replace unlist by a custom list max?
phmax <- sapply(ph,function(v) max(unlist(v)))

# peers2
# ------
qh <- peers2(d,x30)
ql <- sapply(qh,function(.df) { if (is.null(.df)) 0 else height(.df) }
save(qh,file="track-qh.rda")
sum((pl[pl>0]+1)!=ql[ql>0])
# [1] 15532 -- the difference is with multiple neighbor records in delta t
names(ql) <- 1:length(ql)
ql2 <- ql[ql==2]
ql2i <- as.integer(names(ql2))
sum(pl1i!=ql2i) # =>0


pseq <- function(x,i) { 
  # when row.{p.q} where stored as.numeric instead of as.integer,
  # they were floats, and as.character could cause e.g. "9e+05"
  # -- which then fetched NA when used as index below; hence we
  # use format(as.character(...),scientific=F) for that case.
  # Once the d[] is re-made with as.integer for rows,
  # the format(as.numeric,...) can be replaced by as.character back 
  xc <- format(as.numeric(x[i,]),scientific=F,trim=T)
  #xc <- as.character(x[i,])

  row <- dprn[[xc[1]]][xc[2]]
  start <- dp[[xc[1]]]$start[row] # or [row,]$start
  c(row,start)
}


x.pairs <- function(x,i) {
  if (height(x)<=2) list(x)
  else {
    r1 <- drn[i]
    p1 <- d$person[i]
    # if we do not reorder in peers2 above,
    # the first row is always the x1
    # -- this will be the OCaml way
    x1 <- x[x$row==r1,]
    if (x1$person!=p1) stop("x1 person != p1 at i=>",i)
    x0 <- x[x$person!=p1,]
    x0h <- height(x0)
    res <- list()
    for (j in 1:x0h) {
      x2 <- x0[j,]
      if (x1$person<x2$person) xx <- rbind(x1,x2)
      else if (x1$person>x2$person) xx <- rbind(x2,x1)
      else stop("x1$person==x2$person at i=>",i)
      res[[j]] <- xx
    }
    res
  }
}

ql9  <- ql[ql>0]
ql9i <- as.integer(names(ql9))


do.x <- function(x) {
  seq.start.list <- lapply(1:2,function(i) pseq(x,i))
  seq.starts <- do.call(rbind,seq.start.list)
  rows <- x$row
  rhash <- paste(sort(rows),collapse=" ")
  if (is.null(rows_hash[[rhash]])) {
    # NB replacing rbind by cbind above didn't change interleaving order --
    # had to transpose rows to append sequence numbers first, times second!
    v <- c(x$person,rows,seq.starts)
    if (length(v)!=8) stop("v is too long at i=",i)
    last <<- last + 1
    # qq <- rbind(qq,v)
    qq[[last]] <<- v
    rows_hash[[rhash]] <<- 1
    if (last %% 1000 == 0) print(last)
  } 
  else { 
    rows_hash[[rhash]] <<- rows_hash[[rhash]] + 1
  }
}

diff.ts <- function(ts) {
  len <- length(ts)
  if (len>1) c(0,sapply(2:length(ts),function(i) ts[i]-ts[i-1]))
  else if (len==1) 0 
  else numeric(0)
}

do.track.df <- function(.df) {
  dt     <- .df$t.p-.df$t.q
  # those give use delta for the 1st cell in a trajectory,
  # not the last -- which we really want; we need 
  dt.p   <- diff.ts(.df$t.p)
  dt.q   <- diff.ts(.df$t.q)
  dseq.p <- diff.ts(.df$seq.p)
  dseq.q <- diff.ts(.df$seq.q)
  pos    <- 1:height(.df)
  cbind(pos=pos,subset(.df,select=c(seq.p,seq.q,t.p,t.q)),dt=dt,
        dt.p=dt.p,dt.q=dt.q,dseq.p=dseq.p,dseq.q=dseq.q)
}

join.tracks <- function(.df) {
  td <- with(.df,.df[dseq.p == 1 & dseq.q == 1,])
  if (height(td)==0) NULL
  else td 
}

rle.absolutes <- function(a.rle) {
  len <- length(a.rle$lengths)
  # c(1,(cumsum(a.rle$lengths)+1)[1:(len-1)])
  head(c(1,cumsum(a.rle$lengths)+1),-1)
}

df.rle <- function(a.rle) {
  absolut <- rle.absolutes(a.rle)
  data.frame(len=a.rle$lengths,val=a.rle$values,pos=absolut)
}

seq.rle <- function(.df) {
  p <- df.rle(rle(.df$dseq.p))
  q <- df.rle(rle(.df$dseq.q))
  list(p=p,q=q)
}

# NB: rewrite this to recognize segments of 1 
# interrupted by 2s or 3s as a whole!  But,
# .df[val<=2 & len>1,] -- gets 0s and 2s mostly
seq.filter <- function(.df) {
  with(.df,.df[val==1 & len>1,])
}


traject.rles <- function(tr) {
  p <- seq.filter(tr$p)
  q <- seq.filter(tr$q)
  
  if (height(p)==0) p <- NULL
  if (height(q)==0) q <- NULL
  
  if (is.null(p) && is.null(q)) NULL
  else list(p=p,q=q)
}

# here I was going to write a function to intersect integer ranges:
#merge.rles <- function(tj) {}
# -- similar to OCaml's lSet; but when I've googled it, I found 
# -- Perl's Set::IntSpan::Fast, 
# -- and rseek.org brought back IRanges!

library(IRanges)
common.seqs <- function(tj) {
  p <- tj$p
  q <- tj$q
  ir.p <- IRanges(p$pos,p$pos+p$len-1)
  ir.q <- IRanges(q$pos,q$pos+q$len-1)
  intersect(ir.p,ir.q)
}

weed.1.widths <- function(.df) {
  .df <- .df[.df$width>1,]
  if (height(.df)==0) NULL
  else .df
}

# an immediate dataframe-rbind way is slow; 
# we make a list instead, then do.call(rbind,it)
qq <- list()
last <- 0
rows_hash <- list()
count  <- 0
# qh[2890] is a multiple; ql9i[34]==2890
for (i in ql9i) { # or: pl1i; ql2i; ql2i[1:10]; ql9i[1:35] for testing 
  x <- qh[[i]]
  x.h <- height(x)
  if (x.h>2) {
    devnull <- lapply(x.pairs(x,i),do.x)
  } else {
    do.x(x)
  }
}
qd <- as.data.frame(do.call(rbind,qq))
names(qd) <- c("p","q","row.p","row.q","seq.p","seq.q","t.p","t.q")

pq <- sapply(1:height(qd),function(i)paste(qd$p[i],qd$q[i],collapse=" "))
qd <- cbind(qd,pq=pq)
qt <- split(qd,pq,drop=T)
qtl <- sapply(qt,height)

dq <- lapply(qt,do.track.df)
tr <- lapply(dq,seq.rle)

tj0     <- lapply(tr,traject.rles)
tj0null <- sapply(tj0,is.null)
tj      <- Filter(function(x)!is.null(x) && !is.null(x$p) && !is.null(x$q),tj0)

tc.ir <- lapply(tj,common.seqs)
# the names of the dataframe are just 1:height below,
# might as well use for something and compare to ff's:
tc1 <- Filter(function(e)height(e)>0,lapply(tc.ir,as.data.frame))
tc  <- Filter(not.null,lapply(tc1,weed.1.widths))
tt.maxw <- sapply(tt,function(.df)max(.df$width))
max(tt.maxw)
which(tt.maxw==max(tt.maxw))

# size as height shows the number of segments
tt.nsegs <- sapply(tt,height)
tt.decs <- sort(tt.nsegs,decreasing=T)

hist(tt.nsegs,breaks=100)

tt.noot <- tt
tt.noot[["35 88"]] <- NULL
tt.noot[["57 86"]] <- NULL
tt.noot[["16 37"]] <- NULL
tt.noot.nsegs <- sapply(tt.noot,height)
hist(tt.noot.nsegs,breaks=1:50)
tt.noot.decs <- sort(tt.noot.nsegs,decreasing=T)
head(tt.noot.decs,10)
tail(tt.noot.decs,10)
tt.nsegs.rle <- rle(sort(tt.nsegs))
tt.nsegs.rle

tt.avlen <- sapply(tt.noot,function(.df)mean(.df$width))

##### ----- verification of RLE method 
##### with direct dataframe-filtering


ranges.append <- function(.r) {
  .r$ranges[[length(.r$ranges)+1]] <- c(.r$first,.r$last)
  .r
}
rangify <- function(v) {
    .r <- list(ranges=list())
    if (length(v)==0) .r
    else if (length(v)==1) {
      # this can be handled by Reduce/ranges.append, but here is just:
      .r$ranges[[1]] <- c(v[1],v[1])
      .r
    } else {
      .r <- Reduce(
      function(.r,point) { 
        if (point==.r$last+1) { .r$last <- point; .r }
        else { 
          .r <- ranges.append(.r)
          .r$first <- point; .r$last <- point; 
          .r
        }
      },tail(v,-1),list(ranges=list(),first=v[1],last=v[1]))
      .r <- ranges.append(.r)
    }
    .df <- as.data.frame(do.call(rbind,.r$ranges))
    names(.df) <- c("start","end")
    width <- .df$end - .df$start + 1
    cbind(.df,width=width)
}

#debug.count <- 0
df.rangify <- function(.df) {
#  debug.count <<- debug.count+1
#  print(debug.count)
  r <- rangify(.df$pos)
  r <- r[r$width>1,]
  if (height(r)==0) NULL
  else r
}

df <- lapply(dq,join.tracks)
not.null <- function(x)!is.null(x)
df <- Filter(not.null,df)

# at this point we found that some positions were NA
# and then went to find out the 
# as.numeric=>"9e+05"=>as.character problem...

df.pos <- lapply(df,function(.df).df$pos)
df.null <- sapply(df.pos,function(v)which(is.na(v)))
df.null.len <- sapply(df.null,length)
sum(df.null.len) # this must be 0 if things are right!
    # [1] 3
    which(df.null.len!=0)
    # 16 37 68 89 
    #   115   719 
    # observe dq[["68 89"]] -- NA
    # qt -- NA -- qd["93807",] contains row.p: 9e+05
    which(rownames(d)=="900000")
    # [1] 775153
    qh[[775153]]
    #   person    row
    # 1     68 900000
    # 2     89 603476
    which(is.na(qd$seq.p))
    # 252063
    which(is.na(qd$seq.q))
    # 60560
    qd[60560,]
    #        p  q   row.p row.q seq.p seq.q        t.p t.q    pq
    # 60560 73 89 1643399 6e+05  7312    NA 1097617331  NA 73 89
    which(rownames(d)=="600000")
    pseq(qh[[576898]],2)
    #    600000            
    #     12289 1097617320 

    # BTW: qd rownames are just 1:height(qd), as it's hand-made;
    # thus qd["x",] === qd[x,]
    
# NB: in order to fix this properly, we need rows in qh to be integers,
# so they work in as.character in pseq; or format there;
# in order for the rows in qh to be as.integer, we've changed peers2
# to do rnames <- as.integer(rownames(d)) -- a next rerun would do it!
# for now, we've manually fixed qd as shown above, and reran from dq on
# now we have ff similar to tt (vertical mirror images)

ff <- lapply(df,df.rangify)

# NB: truehist
# next steps: 
# -- find timings of each trajectory
# -- find cells traveled in each trajectory
# -- find most frequent edges, pairs of edges, etc.;
# compare with the overall edge distribution
# -- periodicity (classes/commute together)
# -- most common cells for each individual 
# and how they figure in their trajectories
# -- average length of a track and their number
# as a metric of significant co-travel
# -- restore the terrain by restoring triangles 
#   A
#  / \
# B---C
# -- with single edges connecting nodes
# then mosaic them into a planar graph
# -- groups of more than two!

tt[[1]]
# 114 116
names(tt[1]) # "10 21"
qd.rows <- rownames(dq[["10 21"]][114:116,])
#[1] "38969" "38971" "38972"
qd[qd.rows,] # times!
d[as.character(unique(qd[qd.rows,]$row.p)),] # cells!

##### analysis of the average trajectory co-travel
tt.avlen <- sapply(tt.noot,function(.df)mean(.df$width))
tt.avlen.decr <- sort(tt.avlen,decr=T)
plot(tt.avlen.decr)
# we can see an outlier dominating
tt.avlen.decr[1:10]
tt.avlen.decr[1:10]
# 24 26 16 18 35 51 51 88  6 15  9 78 34 69 47 82 50 73 22 35 
# 18.40  4.00  4.00  4.00  4.00  4.00  3.50  3.50  3.50  3.36 
tt[["24 26"]]
plot(tail(tt.avlen.decr,-1))

# NB next: with ff, we have access to timings easily
# let's see what the average step size and total duration
# of a trajectory are, and measure the gaps between segments
# -- to glue likely candidates
# -- timings can help checking where seq jumps are small
# NB single co-occurrences are different from trajectories!

# the contents of vt$starttime is class POSIXct
# but making it back via as.POSIXct doesn't work regardless of tz
# so we use as.POSIXlt without a tz below:
time.of.int <- function(it) as.POSIXlt(it,origin="1970-01-01")

# this operates on a list element, ff[i]
ff.timings <- function(ldf) {
  pair <- names(ldf)
  meat <- ldf[[1]]
  mydf <- df[[pair]]
  # might record both p's and q's, but now doing p's only:
  meat <- cbind(meat,st.int=NA,et.int=NA,dt1.p=NA,dt1.q=NA,gap=NA)
  for (j in 1:height(meat)) {
    row.start <- which(mydf$pos==meat$start[j])
    # row.end <- which(mydf$pos==meat$end[j])
    row.end <- row.start + meat$width[j] - 1 # faster
    meat$st.int[j] <- as.integer(mydf$t.p[row.start])
    meat$et.int[j] <- as.integer(mydf$t.p[row.end])
    meat$dt1.p[j] <- mydf$dt.p[row.start]
    meat$dt1.q[j] <- mydf$dt.q[row.start]
    if (j>1) meat$gap[j] <- meat$st.int[j]-meat$et.int[j-1]
  }
  meat
}

ft <- lapply(1:length(ff),function(i)ff.timings(ff[i]))
names(ft) <- names(ff)
ft <- lapply(ft,function(.df)cbind(.df,st=time.of.int(.df$st.int),et=time.of.int(.df$et.int)))
# dt=difftime(.df$et,.df$st,units="mins") -- trying to do difftime wants origin again?!
ft <- lapply(ft,function(.df)cbind(.df,dt=.df$et.int-.df$st.int)) 

ft.nseg <- sapply(ft,height)
names(ft.nseg) <- names(ft)
ft.nseg.za <- sort(ft.nseg,decr=T)
ft.nseg.za[1:10]

ft.av.hops <- sapply(ft,function(.df)mean(.df$end-.df$start))
names(ft.av.hops) <- names(ft)
ft.av.hops.za <- sort(ft.av.hops,decr=T)
ft.av.hops.za[1:10]

ft.av.dt <- sapply(ft,function(.df)mean(.df$dt))
names(ft.av.dt) <- names(ft)
ft.av.dt.za <- sort(ft.av.dt,decr=T)
ft.av.dt.za[1:10]

ft.max.dt <- sapply(ft,function(.df)max(.df$dt))
names(ft.max.dt) <- names(ft)
ft.max.dt.za <- sort(ft.max.dt,decr=T)
ft.max.dt.za[1:10]

ft.dt <- data.frame(max.dt=ft.max.dt,av.dt=ft.av.dt,nseg=ft.nseg,av.hops=ft.av.hops)
ft.dt <- with(ft.dt,ft.dt[order(-max.dt,-nseg),])
ft.dt[1:40,]