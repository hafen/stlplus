co2_stl <- stlplus(co2, t = as.vector(time(co2)), n.p = 12,
  l.window = 13, t.window = 19, s.window = 35, s.degree = 1,
  sub.labels = substr(month.name, 1, 3))

plot(co2_stl, ylab = "CO2 Concentration (ppm)", xlab = "Time (years)")
plot_seasonal(co2_stl)
plot_trend(co2_stl)
plot_cycle(co2_stl)
plot_rembycycle(co2_stl)

# post-trend smoothing

co2_stl_pt <- stlplus(co2, t = as.vector(time(co2)), n.p = 12,
  l.window = 13, t.window = 19, s.window = 35, s.degree = 1,
  sub.labels = substr(month.name, 1, 3),
  fc.degree = c(1, 2), fc.window = c(201, 35),
  fc.name = c("long-term", "so. osc."))

plot(co2_stl_pt, scales = list(y = list(relation = "free")),
  ylab = "CO2 Concentration (ppm)", xlab = "Time (years)",
  aspect = 0.25, type = c("l", "g"))

# with NAs

y <- co2
y[201:224] <- NA

y_stl <- stlplus(y, l.window = 13, t.window = 19, s.window = 35,
  s.degree = 1, sub.labels = substr(month.name, 1, 3))

plot(y_stl, ylab = "CO2 Concentration (ppm)", xlab = "Time (years)", type = c("l", "g"))
plot_seasonal(y_stl)
plot_trend(y_stl)
plot_cycle(y_stl)
plot_rembycycle(y_stl)

# if you don't want to use a time series object:
y_stl <- stlplus(y, t = as.vector(time(y)), n.p = 12,
  l.window = 13, t.window = 19, s.window = 35, s.degree = 1,
  sub.labels = substr(month.name, 1, 3))

# with an outlier
y2 <- co2
y2[200] <- 300

y2_stl <- stlplus(y2, t = as.vector(time(y2)), n.p = 12,
  l.window = 13, t.window = 19, s.window = 35, s.degree = 1,
  sub.labels = substr(month.name, 1, 3), outer = 10)

plot(y2_stl, ylab = "CO2 Concentration (ppm)", xlab = "Time (years)")
plot_seasonal(y2_stl)
plot_trend(y2_stl)
plot_cycle(y2_stl)
plot_rembycycle(y2_stl)

# compare to R's stl

x1 <- stlplus(co2, t = as.vector(time(co2)), n.p = 12,
  l.window = 13, t.window = 19, s.window = 11, s.degree = 1,
  sub.labels = substr(month.name, 1, 3))

x2 <- stl(co2, l.window = 13, t.window = 19, s.window = 11, s.degree = 1)

# will be different due to interpolation differences
plot(seasonal(x1) - seasonal(x2))

# but not if all jump parameters are 1
x1 <- stlplus(co2, t = as.vector(time(co2)), n.p = 12,
  l.window = 13, t.window = 19, s.window = 11, s.degree = 1,
  sub.labels = substr(month.name, 1, 3),
  s.jump = 1, t.jump = 1, l.jump = 1)

x2 <- stl(co2, l.window = 13, t.window = 19, s.window = 11, s.degree = 1,
  s.jump = 1, t.jump = 1, l.jump = 1)

plot(seasonal(x1) - seasonal(x2))
