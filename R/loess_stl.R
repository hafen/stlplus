.loess_stl2 <- function(x=NULL, y, span, degree, weights=NULL, m=c(1:length(y)), y_idx=!is.na(y), noNA=all(y_idx), blend=0, jump=ceiling(span/10), at=c(1:length(y))) {

	nextodd <- function(x) {
		x <- round(x)
      x2 <- ifelse(x%%2==0, x+1, x)
		as.integer(x2)
	}

   n <- length(y[y_idx])
   if(is.null(x)) x <- c(1:length(y))

   if(is.null(weights)) weights <- rep(1, length(y))
   n_m <- length(m)
   
   if((span %% 2) == 0) {
      span <- span + 1
      warning(paste("Span must be odd! Changed span from ", span-1, " to ", span, sep=""))
   }
   
   s2 <- (span + 1)/2
   # set up indices in R - easier
   if(noNA) {
      if((diff(range(x))) < span) {
         l_idx <- rep(1, n_m)
         r_idx <- rep(n, n_m)
      } else{
         l_idx <- c(rep(1, length(m[m<s2])), m[m>=s2 & m <= n-s2]-s2+1, rep(n - span + 1, length(m[m>n-s2])))
         r_idx <- l_idx + span - 1         
      }
      aa <- abs(m - x[l_idx])
      bb <- abs(x[r_idx] - m)
      max_dist <- ifelse(aa > bb, aa, bb)
      # max_dist <- apply(cbind(abs(m - x[l_idx]), abs(x[r_idx] - m)), 1, max)
   } else {
      span3 <- min(span, n)
      x2 <- x[y_idx]
      
      # another approach
      a <- ann(ref=as.matrix(x2), target=as.matrix(m), tree.type="kd", k=span3, eps=0, verbose=FALSE)$knnIndexDist[,1:span3]

      l_idx <- apply(a, 1, min)
      r_idx <- apply(a, 1, max)
      
      # # if there are NAs, getting nearest neighbors is a little more difficult...
      # # speed this up later, but something that works will suffice for now
      # x2 <- x[y_idx]
      # # idea is to get the distances of all observations to each point in m and get the span closest      
      # a <- by(m, list(m), function(a) {
      #    order(abs(x2 - a))[1:span]
      # })
      # l_idx <- sapply(a, min)
      # r_idx <- sapply(a, max)
      
      max_dist <- apply(cbind(abs(m - x2[l_idx]), abs(x2[r_idx] - m)), 1, max)
   }
   if(span >= n)
      # max_dist <- max_dist * (span/n)
      max_dist <- max_dist + (span-n)/2
   
   out <- .C("loess_stl2", as.double(x[y_idx]), as.double(y[y_idx]), as.integer(n), as.integer(degree), as.integer(span), as.double(weights[y_idx]), as.integer(m), as.integer(n_m), as.integer(l_idx-1), as.integer(r_idx-1), as.double(max_dist), result=double(n_m), slopes=double(n_m), PACKAGE="stl2")
   
   res1 <- out$result
   # do interpolation
   if(jump > 1)
      res1 <- .stl2_interp(m, out$result, out$slope, at)
      # res1 <- approx(x=m, y=out$result, xout=at)$y
   
   if(blend > 0 && blend <= 1 && degree >= 1) {
      if(degree == 2)
         sp0 <- nextodd((span+1)/2)
      if(degree == 1)
         sp0 <- span

      n.b <- as.integer(span/2)

      blend <- 1 - blend # originally programmed backwards - easier to fix this way
      # indices for left and right blending points
      # take into account if n_m is too small
      mid <- median(m)
      bl_idx <- m <= n.b + jump  & m < mid
      br_idx <- m >= n - n.b - jump + 1 & m >= mid
      left <- m[bl_idx]
      right <- m[br_idx]
      bl_idx_interp <- at <= max(left)
      br_idx_interp <- at >= min(right)
      left_interp <- at[bl_idx_interp]
      right_interp <- at[br_idx_interp]
      # left_interp <- at[bl_idx_interp]
      # right_interp <- at[br_idx_interp]
      l_idx2 <- l_idx[bl_idx | br_idx]
      r_idx2 <- r_idx[bl_idx | br_idx]
      max_dist2 <- max_dist[bl_idx | br_idx]

      m2 <- c(left, right)
      n_m2 <- length(m2)
      
      # speed this up later by only getting the loess smooth at the tails.
      # right now, a lot of unnecessary calculation is done at the interior
      # where blending doesn't matter

      tmp <- .C("loess_stl2", as.double(x[y_idx]), as.double(y[y_idx]), as.integer(n), as.integer(0), as.integer(sp0), as.double(weights[y_idx]), as.integer(m2), as.integer(n_m2), as.integer(l_idx2-1), as.integer(r_idx2-1), as.double(max_dist2), result=double(n_m2), slopes=double(n_m2), PACKAGE="stl2")
      
      if(jump > 1) {
         res2_left <- .stl2_interp(left, head(tmp$result, length(left)), head(tmp$slope, length(left)), left_interp)
         res2_right <- .stl2_interp(right, tail(tmp$result, length(right)), tail(tmp$slope, length(right)), right_interp)
      } else {
         res2_left <- head(tmp$result, length(left))
         res2_right <- tail(tmp$result, length(right))
      }
      # res2 <- approx(x=m, y=tmp$result, xout=at)$y

      p.left <- ((1 - blend)/(n.b-1))*(left_interp - 1) + blend
      p.right <- ((blend - 1)/(n.b - 1))*(right_interp - (n - n.b + 1)) + 1
      p.left[p.left < blend] <- blend
      p.left[p.left > 1] <- 1
      p.right[p.right < blend] <- blend
      p.right[p.right > 1] <- 1

      res1[bl_idx_interp] <- res1[bl_idx_interp] * p.left + res2_left * (1 - p.left)
      res1[br_idx_interp] <- res1[br_idx_interp] * p.right + res2_right * (1 - p.right)
            
      # xxx <- x[y_idx]
      # yyy <- y[y_idx]
      # tmp2 <- predict(loess(yyy ~ xxx, deg=0, span=(sp0+0.00000001)/length(yyy), control=loess.control(surface="direct")), newdata=m2)
      # tmp3 <- predict(loess(yyy ~ xxx, deg=degree, span=(span+0.00000001)/length(yyy), control=loess.control(surface="direct")), newdata=m)
      
   }
   
   res1
}

