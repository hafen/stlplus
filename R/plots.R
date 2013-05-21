plot.stl2 <- function(x, scales=list(y=list(relation="sliced")), type="l", as.table=TRUE, strip=FALSE, strip.left=TRUE, between=list(y=0.5), layout=NULL, ...) {

   if(x$pars$fc.number == 0) {
      d <- data.frame(time=rep(time.stl2(x), 4), values=c(getraw(x), seasonal(x), trend(x), remainder(x)), ind=factor(rep(c(1:4), each=x$n)))
      levels(d$ind) <- c("raw", "seasonal", "trend", "remainder")
      d$which <- "isnotNA"
      
      if(is.null(layout)) layout <- c(1, 4)

      # if there are any NA values, plot the part of the seasonal and trend 
      # components with a different color
      if(any(is.na(getraw(x)))) {
         d$values[d$ind=="seasonal" & is.na(getraw(x))] <- NA
         d$values[d$ind=="trend" & is.na(getraw(x))] <- NA
         d <- rbind(d, 
            data.frame(time=rep(time.stl2(x), 2), values=c(seasonal(x), trend(x)), ind=rep(c("seasonal", "trend"), each=x$n), which="isxNA")
         )
         d$values[d$ind=="seasonal" & d$which=="isxNA" & !is.na(getraw(x))] <- NA
         d$values[d$ind=="trend" & d$which=="isxNA" & !is.na(getraw(x))] <- NA      
      }

   } else {
      fc.number <- x$pars$fc.number
      nvar <- 3 + fc.number
      fc.name <- as.character(x$pars$fc$fc.name)
      if(fc.number == 1) {
         fcdat <- x$fc[,1]
      } else {
         fcdat <- stack(x$fc[,1:fc.number])$values
      }
      
      d <- data.frame(time=rep(time.stl2(x), nvar), values=c(getraw(x), seasonal(x), fcdat, remainder(x)), ind=factor(rep(c(1:nvar), each=x$n)))
      levels(d$ind) <- c("raw", "seasonal", fc.name, "remainder")
      d$which <- "isnotNA"
      
      if(is.null(layout)) layout <- c(1, nvar)
      
      # if there are any NA values, plot the part of the seasonal and trend 
      # components with a different color
      if(any(is.na(getraw(x)))) {
         d$values[d$ind=="seasonal" & is.na(getraw(x))] <- NA
         for(i in 1:fc.number) {
            d$values[d$ind==fc.name[i] & is.na(getraw(x))] <- NA
         }
         d <- rbind(d, 
            data.frame(time=rep(time.stl2(x), fc.number + 1), values=c(seasonal(x), fcdat), ind=rep(c("seasonal", fc.name), each=x$n), which="isxNA")
         )
         d$values[d$ind=="seasonal" & d$which=="isxNA" & !is.na(getraw(x))] <- NA
         for(i in 1:fc.number) {
            d$values[d$ind==fc.name[i] & d$which=="isxNA" & !is.na(getraw(x))] <- NA            
         }
      }      
   }
   

	p <- xyplot(values ~ time | ind,
		data=d,
      groups=which,
		type=type,
		layout=layout,
		scales=scales,
		as.table=as.table,
		strip=strip, strip.left=strip.left,
		between=between,
		...
	)
	p
}

.midmean <- function(x) {
	q <- quantile(x, probs=c(0.25, 0.75), na.rm=TRUE)
	isInner <- x<q[2] & x > q[1]
	mean(x[isInner], na.rm=TRUE)
}

plot.cycle <- function(x, layout=c(x$pars$n.p, 1), col="#0080ff", xlab="Time", ylab="Seasonal", panel= function(x, y, ...) {
   panel.segments(x, rep(.midmean(y), length(x)), x, y, col=col)
}, ...) {

   seas <- seasonal(x)
   t <- time.stl2(x)

   cycleSubIndices <- x$data$sub.labels

   if(x$pars$periodic) {
      
      p <- xyplot(seas ~ t | cycleSubIndices, 
         layout=layout, 
         type="l",
         xlab=xlab,
         ylab=ylab,
   		...
      )
   } else {
      
      p <- xyplot(seas ~ t | cycleSubIndices, 
         layout=layout, 
         type="l",
         panel=panel,
         xlab=xlab,
         ylab=ylab,
   		...
      )
   }
   p
}


plot.seasonal <- function(x, col=c("darkgray", "black"), lwd=2, xlab="Time", ylab="Centered Seasonal + Remainder", ...) {
   
   dat <- by(data.frame(v1=seasonal(x) + remainder(x), v2=seasonal(x), t=time.stl2(x), v3=x$data$sub.labels), list(x$data$sub.labels), function(dd) {
      mn <- mean(dd$v1, na.rm=TRUE)
      data.frame(a=dd$v1 - mn, b=dd$v2 - mn, t=dd$t, sub.labels=dd$v3)
   })

   dat <- do.call(rbind, dat)

   p <- xyplot(a + b ~ t | sub.labels, data=dat,
		type=c("p", "l"),
		col=col,
		lwd=lwd,
		distribute.type=TRUE,
      as.table=TRUE,
      xlab=xlab,
      ylab=ylab,
		...
   )
   p
}

# old plot.seasonal...
# function (x, col = c("black", "red"), lwd = 2, xlab = "Time", 
#     ylab = "Centered Seasonal + Remainder", ...) 
# {
#     cycleSubIndices <- rep(c(1:x$pars$n.p), ceiling(x$n/x$pars$n.p))[1:x$n]
#     vals <- stack(by(seasonal(x) + remainder(x), list(cycleSubIndices), 
#         function(x) {
#             x - mean(x, na.rm = TRUE)
#         }))
#     vals2 <- stack(by(remainder(x), list(cycleSubIndices), identity))
#     vals$ind <- as.numeric(as.character(vals$ind))
#     vals$ind <- factor(vals$ind, labels = x$sub.labels)
#     if (class(x$t) == "Date") {
#         t <- stack(by(as.integer(time.stl2(x)), list(cycleSubIndices), 
#             identity))
#         t$values <- as.Date(t$values, origin = "1970-01-01")
#     }
#     else {
#         t <- stack(by(time.stl2(x), list(cycleSubIndices), identity))
#     }
#     vals$t <- t$values
#     vals$seasonal <- vals$values - vals2$values
#     vals <- subset(vals, !is.na(values))
#     p <- xyplot(values + seasonal ~ t | ind, data = vals, type = c("p", 
#         "l"), col = col, lwd = lwd, distribute.type = TRUE, as.table = TRUE, 
#         xlab = xlab, ylab = ylab, ...)
#     p
# }

plot.rembycycle <- function(x, col="darkgray", locol="black", lolwd=2, xlab="Time", ylab="Remainder", ...) {

   vals2 <- data.frame(values=remainder(x), ind=x$data$sub.labels, t=time.stl2(x))

   # vals2 <- stack(by(remainder(x), list(cycleSubIndices), identity))
   #    vals2$ind <- as.numeric(as.character(vals2$ind))
   #    vals2$ind <- factor(vals2$ind, labels=x$sub.labels)

   # if(class(x$t) == "Date") {
   #    t <- stack(by(as.integer(time.stl2(x)), list(cycleSubIndices), identity))
   #  t$values <- as.Date(t$values, origin="1970-01-01") 
   # } else {
   #       t <- stack(by(time.stl2(x), list(cycleSubIndices), identity))      
   # }
	
   vals2 <- subset(vals2, !is.na(values))

   p <- xyplot(values ~ t | ind, 
      data=vals2, 
		type="p",
      panel=function(x, y, ...) {
         panel.abline(h=0, lty=2, col="darkgray")
         panel.xyplot(x, y, ...)
         panel.loess(x, y, col=locol, lwd=lolwd)
      },
      xlab=xlab,
      ylab=ylab,
      col=col,
      as.table=TRUE,
      ...
   )
   p
}

plot.trend <- function(x, xlab="Time", ylab="Trend", span=0.3, type=c("p", "l"), scales=list(y=list(relation="free")), lwd=c(1, 1), col=c("darkgray", "black", "darkgray"), layout=c(1, 2), between=list(y=0.5), strip=FALSE, strip.left=TRUE, as.table=TRUE, ...) {

   dat <- rbind(
      data.frame(x=time.stl2(x), y=remainder(x) + trend(x), type="p", pan="Trend"),
      data.frame(x=time.stl2(x), y=trend(x), type="l", pan="Trend"),
      data.frame(x=time.stl2(x), y=remainder(x), type="p", pan="Remainder"),
      data.frame(x=time.stl2(x), y=predict(loess(remainder(x) ~ c(1:length(remainder(x))), span=span, weights=x$data$weights), newdata=c(1:length(remainder(x)))), type="l", pan="Remainder")
   )

   p <- xyplot(y ~ x | pan,
      groups=type,
      data=dat,
		type=type,
		col=col,
      as.table=as.table,
      layout=layout,
      lwd=lwd,
      scales=scales,
      distribute.type=TRUE,
      between=between,
      strip=strip, strip.left=strip.left,
      xlab=xlab, ylab=ylab,
      ...
   )
   p
}

