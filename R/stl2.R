# x <- as.numeric(co2)
# t <- as.vector(time(co2))
# n.p <- 12; l.window <- 13; t.window <- 19; s.window <- 35; s.degree <- 1; t.degree=1
# sub.labels <- substr(month.name, 1, 3)
# outer <- 1; inner <- 2; details <- FALSE; robust <- FALSE
# l.degree <- t.degree
# s.jump <- 10; t.jump <- 1; l.jump <- 1

# x <- dailycount$Freq/1000
# n.p <- 7; l.window <- NULL; t.window <- 31; s.window <- 51; s.degree <- 1; t.degree=1
# sub.labels <- substr(month.name, 1, 3)
# outer <- 1; inner <- 1; details <- TRUE; robust <- FALSE
# l.degree <- t.degree

stl2 <- function(x, t=NULL, n.p, s.window, s.degree=1, t.window=NULL, t.degree=1, fc.window=NULL, fc.degree=NULL, fc.name=NULL, l.window=NULL, l.degree=t.degree, s.jump=ceiling(s.window/10), t.jump=ceiling(t.window/10), l.jump=ceiling(l.window/10), fc.jump=NULL, critfreq=0.05, s.blend=0, t.blend=0, l.blend=t.blend, fc.blend=NULL, inner=2, outer=1, sub.labels=NULL, sub.start=1, zero.weight=1e-6, details=FALSE, ...) UseMethod("stl2")

stl2.ts <- function(x, t=as.numeric(time(x)), n.p=frequency(x), s.window, s.degree=1, t.window=NULL, t.degree=1, fc.window=NULL, fc.degree=NULL, fc.name=NULL, l.window=NULL, l.degree=t.degree, s.jump=ceiling(s.window/10), t.jump=ceiling(t.window/10), l.jump=ceiling(l.window/10), fc.jump=NULL, critfreq=0.05, s.blend=0, t.blend=0, l.blend=t.blend, fc.blend=NULL, inner=2, outer=1, sub.labels=NULL, sub.start=1, zero.weight=1e-6, details=FALSE, ...) {

   if (is.matrix(x)) 
       stop("only univariate series are allowed")

   if(missing(n.p)) n.p <- frequency(x)
   if(!is.null(t)) {
      if(length(t) != length(x)) {
         stop("t must be same length as time series")
      }
   } else {
      t <- as.vector(time(x))
   }

   stl2.default(x, t=t, n.p=n.p, s.window=s.window, s.degree=s.degree, t.window=t.window, t.degree=t.degree, fc.window=fc.window, fc.degree=fc.degree, fc.name=fc.name, l.window=l.window, l.degree=l.degree, s.jump=s.jump, t.jump=t.jump, l.jump=l.jump, fc.jump=fc.jump, critfreq=0.05, s.blend=s.blend, t.blend=t.blend, l.blend=l.blend, fc.blend=NULL, inner=inner, outer=outer, sub.labels=sub.labels, sub.start=sub.start, details=details, ...)
}


stl2.zoo <- function(...) {
   stl2.ts(...)
}



stl2.default <- function(x, t=NULL, n.p, s.window, s.degree=1, t.window=NULL, t.degree=1, fc.window=NULL, fc.degree=NULL, fc.name=NULL, l.window=NULL, l.degree=t.degree, s.jump=ceiling(s.window/10), t.jump=ceiling(t.window/10), l.jump=ceiling(l.window/10), fc.jump=NULL, critfreq=0.05, s.blend=0, t.blend=0, l.blend=t.blend, fc.blend=NULL, inner=2, outer=1, sub.labels=NULL, sub.start=1, zero.weight=1e-6, details=FALSE, ...) {

   if(missing(n.p)) stop("must specify periodicity of seasonal (either explicitly or through a time series object)")
	n.p <- as.integer(n.p)
	
	if(n.p < 4)
	   stop(paste("Parameter n.p was set to ", n.p, ".  Must be at least 4.", sep=""))

   # family <- ifelse(robust, "symmetric", "gaussian")

	Y <- as.vector(x)
	n <- length(Y)
	nextodd <- function(x) {
		x <- round(x)
      x2 <- ifelse(x%%2==0, x+1, x)
      # if(any(x != x2))
      #    warning("A smoothing span was not odd, was rounded to nearest odd.  Check final object parameters to see which spans were used.")
		as.integer(x2)
	}
	
	
	wincheck <- function(x) {
      x <- nextodd(x)
	   if(any(x <= 0)) stop("Window lengths must be positive.")
      x
	}
	
	degcheck <- function(x) {
	   if(! all(x==0 | x==1 | x==2)) stop("Smoothing degree must be 0, 1, or 2")
	}
		
	get.t.window <- function(t.dg, s.dg, n.s, n.p, omega) {
      if(t.dg == 0) t.dg <- 1
      if(s.dg == 0) s.dg <- 1
      
      coefs_a <- data.frame(a = c(0.000103350651767650, 3.81086166990428e-06
      ), b = c(-0.000216653946625270, 0.000708495976681902))
      coefs_b <- data.frame(a = c(1.42686036792937, 2.24089552678906
      ), b = c(-3.1503819836694, -3.30435316073732), c = c(5.07481807116087, 
      5.08099438760489))
      coefs_c <- data.frame(a = c(1.66534145060448, 2.33114333880815
      ), b = c(-3.87719398039131, -1.8314816166323), c = c(6.46952900183769, 
      1.85431548427732))
      
      # estimate critical frequency for seasonal
      betac0 <- coefs_a$a[s.dg] + coefs_a$b[s.dg] * omega
      betac1 <- coefs_b$a[s.dg] + coefs_b$b[s.dg] * omega + coefs_b$c[s.dg] * omega^2
      betac2 <- coefs_c$a[s.dg] + coefs_c$b[s.dg] * omega + coefs_c$c[s.dg] * omega^2
      f_c <- (1 - (betac0 + betac1 / n.s + betac2 / n.s^2)) / n.p
      
      # choose 
      betat0 <- coefs_a$a[t.dg] + coefs_a$b[t.dg] * omega
      betat1 <- coefs_b$a[t.dg] + coefs_b$b[t.dg] * omega + coefs_b$c[t.dg] * omega^2
      betat2 <- coefs_c$a[t.dg] + coefs_c$b[t.dg] * omega + coefs_c$c[t.dg] * omega^2

      betat00 <- betat0 - f_c

      n.t <- nextodd((-betat1 - sqrt(betat1^2 - 4*betat00*betat2)) / (2 * betat00))
      
      n.t
   }
   	
	y_idx <- !is.na(Y)
	noNA <- all(y_idx)
	
	csSmooth <- function(d) {
      nn <- length(d$y)

		# check to make sure there is at least one valid value
		if(all(is.na(d$y))) {
			a <- rep(NA, 1:(nn+2))
		} else {
         cs.ev <- seq(1, length(d$y), by=s.jump)
         if(tail(cs.ev, 1) != nn) cs.ev <- c(cs.ev, nn)
         cs.ev <- c(0, cs.ev, nn + 1)
         # print(m)
         a <- .loess_stl2(y=d$y, span=s.window, degree=s.degree, weights=d$w, m=cs.ev, noNA=noNA, jump=s.jump, at=c(1:(nn+2)))		   
		}
      c(a, rep(NA, csLength - nn - 2))
   }

	if(is.null(l.window)) {
	   l.window <- nextodd(n.p)
	} else {
	   l.window <- wincheck(l.window)
	}

   if(is.null(sub.labels)) {
      sub.labels <- paste("subseries", 1:n.p)
   } else {
      if(length(sub.labels) != n.p) stop("sub.labels must be of length n.p")
      if(length(unique(sub.labels)) != n.p) stop("sub.labels must be unique")
   }
   tmp <- c(sub.start:n.p, rep(1:n.p, ceiling(n/n.p)))[1:n]
   sub.labels <- factor(tmp, labels=sub.labels)

	periodic <- FALSE
	if (is.character(s.window)) {
		if (is.na(pmatch(s.window, "periodic"))) 
			stop("unknown string value for s.window")
		else {
			periodic <- TRUE
			s.window <- 10 * n + 1
			s.degree <- 0
			s.jump <- ceiling(s.window/10)
		}
	} else {
	   s.window <- wincheck(s.window)
	}
	
	degcheck(s.degree); degcheck(t.degree); degcheck(l.degree)

	if (is.null(t.window)) {	   
      # t.window <- nextodd(ceiling(1.5 * n.p/(1 - 1.5/s.window)))
		t.window <- get.t.window(t.degree, s.degree, s.window, n.p, critfreq)
	} else {
	   t.window <- wincheck(t.window)
	}
	
	if(is.null(s.jump) || length(s.jump)==0) s.jump <- ceiling(s.window/10)
	if(is.null(t.jump) || length(t.jump)==0) t.jump <- ceiling(t.window/10)
	if(is.null(l.jump) || length(l.jump)==0) l.jump <- ceiling(l.window/10)

   # cat(s.degree, " ", t.degree, " ", l.degree, "\n")
   # cat(s.window, " ", t.window, " ", l.window, "\n")

   # Trend vector - initialize to 0 or NA, depending on what's in Y
   trend <- 0

   # start and end indices for after adding in extra n.p before and after
   st <- n.p + 1
   nd <- n + n.p

   # cycleSubIndices will keep track of what part of the seasonal each observation belongs to
   cycleSubIndices <- rep(c(1:n.p), ceiling(n/n.p))[1:n]

	if(any(by(Y, list(cycleSubIndices), function(x) all(is.na(x)))))
		stop("There is at least one subseries for which all values are missing.")

   C <- rep(NA, n + 2*n.p)

   dtls <- NULL

   w <- rep(1, n)

for(o_iter in 1:outer) {

   for(iter in 1:inner) {

      # step 1: detrending...
      Y.detrended <- Y - trend
      
      csLength <- ceiling(n/n.p) + 2
      
      # step 2: smoothing of cycle-subseries
      # mapReduce(map=cycleSubIndices, matrix(csSmooth(Y.detrended, w), nrow=1), data=data.frame(Y.detrended, w, cycleSubIndices))
      # if(periodic) {
      #    for(i in 1:n.p) {
      #       cycleSub <- Y.detrended[cycleSubIndices==i]
      #       cycleSub.length <- length(cycleSub)
      #       cs1 <- head(cycleSubIndices, n.p)
      #       cs2 <- tail(cycleSubIndices, n.p)
      #       C[c(cs1, cycleSubIndices, cs2)==i] <- rep(mean(cycleSub, na.rm=TRUE), cycleSub.length + 2) * c(1, w[cycleSubIndices==i],1)
      #    }
      # } else {
      #    C <- as.vector(
      #       do.call(rbind,
      #          by(data.frame(y=Y.detrended, w=w), list(cycleSubIndices), function(x) csSmooth(x))       
      #       )
      #    )[1:(n+2*n.p)]         
      # }
      
      
      for(i in 1:n.p) {
         cycleSub <- Y.detrended[cycleSubIndices==i]
         subWeights <- w[cycleSubIndices==i]
         cycleSub.length <- length(cycleSub)
         
         cs1 <- head(cycleSubIndices, n.p)
         cs2 <- tail(cycleSubIndices, n.p)
         
         notEnoughData <- length(cycleSub[!is.na(cycleSub)]) < s.window/2
         # if(notEnoughData || periodic) {
         if(periodic) {
            C[c(cs1, cycleSubIndices, cs2)==i] <- rep(weighted.mean(cycleSub, w=w[cycleSubIndices == i], na.rm=TRUE), cycleSub.length + 2)
         } else {
            cs.ev <- seq(1, cycleSub.length, by=s.jump)
            if(tail(cs.ev, 1) != cycleSub.length) cs.ev <- c(cs.ev, cycleSub.length)
            cs.ev <- c(0, cs.ev, cycleSub.length + 1)
            tmps <- .loess_stl2(y=cycleSub, span=s.window, degree=s.degree, m=cs.ev, weights=w[cycleSubIndices == i], blend=s.blend, jump=s.jump,  at=c(0:(cycleSub.length+1)))
           C[c(cs1, cycleSubIndices, cs2)==i] <- tmps
           # approx(x=cs.ev, y=tmps, xout=c(0:(cycleSub.length+1)))$y
         }
      }
      
      # Step 3: Low-pass filtering of collection of all the cycle-subseries
      # moving averages
      ma3 <- .ma(C, n.p)

      l.ev <- seq(1, n, by=l.jump)
      if(tail(l.ev, 1) != n) l.ev <- c(l.ev, n)
      L <- .loess_stl2(y=ma3, span=l.window, degree=l.degree, m=l.ev, weights=w, y_idx=y_idx, noNA=noNA, blend=l.blend, jump=l.jump, at=c(1:n))
      
		# L <- predict(loess(ma3 ~ c(1:n), degree=l.degree, span=l.window/n, family=family), newdata=c(1:n))
      
      # Step 4: Detrend smoothed cycle-subseries
      seasonal <- C[st:nd] - L

      # Step 5: Deseasonalize
      D <- Y - seasonal

      # Step 6: Trend Smoothing

      t.ev <- seq(1, n, by=t.jump)
      if(tail(t.ev, 1) != n) t.ev <- c(t.ev, n)
      trend <- .loess_stl2(y=D, span=t.window, degree=t.degree, m=t.ev, weights=w, y_idx=y_idx, noNA=noNA, blend=t.blend, jump=t.jump, at=c(1:n))
      # if(blend$t) {
      #    # TODO: validate blend parameters
      #    trend0 <- .loess_stl2(y=D, span=t.window, degree=t.degree, m=t.ev, weights=w, y_idx=y_idx, noNA=noNA)   
      # }
      
      # trend <- predict(loess(D ~ c(1:length(D)), degree=t.degree, span=t.window/length(D), family=family), newdata=c(1:length(D)))

      if(details) {
         dtls <- c(dtls, list(list(trend=trend, seasonal=seasonal, C=C, D=D, L=L, Y.detrended=Y.detrended, weights=w)))
      }
   }

   if(outer > 1) {
      mid1 <- floor(n/2+1)
      mid2 <- n - mid1+1
      R <- Y - seasonal - trend
      R.abs <- abs(R)
      h <- 3 * sum(sort(R.abs)[mid1:mid2])
      h9 = .999 * h
      h1 = .001 * h
      w <- (1 - (R.abs / h)^2)^2
      w[R.abs <= h1] <- 1
      w[R.abs >= h9] <- 0
      w[w == 0] <- zero.weight
      w[is.na(w)] <- 1
   }
}
   ## post-seasonal smoothing, if specified
   fc <- NULL
   fc.number <- 0
   fc.res <- NULL
   if(!is.null(fc.window)) {
      fc.number <- length(fc.window)
      fc.window <- wincheck(fc.window)
      if(is.null(fc.degree)) fc.degree <- 1
      if(length(fc.degree) < fc.number)
         fc.degree <- c(fc.degree, rep(fc.degree[length(fc.degree)], fc.number - length(fc.degree)))
      fc.cumulative <- rep(0, n) # keep sum of all previous fc smoothings

      degcheck(fc.degree)

      if(is.null(fc.name))
         fc.name <- paste("fc.", fc.window, sep="")

      if(is.null(fc.jump))
         fc.jump <- ceiling(fc.window/10)

      if(is.null(fc.blend))
         fc.blend <- rep(t.blend, fc.number)
         
      if(length(fc.blend) < fc.number)
         fc.blend <- c(fc.blend, rep(0, fc.number - length(fc.blend)))

      if(length(fc.jump) < fc.number)
         fc.jump <- c(fc.jump, rep(fc.jump[length(fc.jump)], fc.number - length(fc.jump)))
         
      fc.res <- data.frame(matrix(nrow=n, ncol=fc.number, data=0))
      for(ii in 1:fc.number) {      
         fc.ev <- seq(1, n, by=fc.jump[ii])
         if(tail(fc.ev, 1) != n) fc.ev <- c(fc.ev, n)
         tmp <- .loess_stl2(y=Y - seasonal - fc.cumulative, span=fc.window[ii], degree=fc.degree[ii], m=fc.ev, weights=w, y_idx=y_idx, noNA=noNA, blend=fc.blend[ii], jump=fc.jump[ii], at=c(1:n))
         
         fc.cumulative <- fc.cumulative + tmp
         fc.res[,ii] <- tmp
      }
      
      if(any(is.null(fc.name)) || length(fc.name) < fc.number)
         fc.name <- paste("trend", c(1:fc.number))
         
      names(fc.res) <- fc.name
      fc.res$remainder <- Y - seasonal - apply(fc.res, 1, sum)
      
      fc <- data.frame(fc.window=fc.window, fc.degree=fc.degree, fc.name=fc.name, fc.jump=fc.jump, fc.blend=fc.blend)
   }
   
   # compute remainder
   R <- Y - seasonal - trend
	
   a <- data.frame(raw=Y, seasonal=seasonal, trend=trend, remainder=R, weights=w, sub.labels=sub.labels)
   pars <- list(deg=data.frame(s.degree=s.degree, t.degree=t.degree, l.degree=l.degree), win=data.frame(s.window=s.window, t.window=t.window, l.window=l.window), blend=data.frame(s.blend=s.blend, t.blend=t.blend, l.blend=l.blend), jump=data.frame(s.jump=s.jump, t.jump=t.jump, l.jump=l.jump), fc.number=fc.number, fc=fc, n.p=n.p, inner=inner, outer=outer, periodic=periodic)
   b <- list(data=a, fc=fc.res, pars=pars, time=t, n=nrow(a), details=dtls)
   class(b) <- "stl2"
   b
}

