
context("User Input Checks")

input.ok <- rbind(data.table(vals = rnorm(2000, mean = 5.0, sd = 0.75), seg = 1),
               data.table(vals = rnorm(1000, mean = 0.5, sd = 1.0), seg = 1),
              data.table(vals = rnorm(2000, mean = 5.5, sd = 0.75), seg = 2),
              data.table(vals = rnorm(2000, mean = 4.0, sd = 3.0), seg = 2))

parameters.ok <- data.table(
   mean=c(6,4),
   nu=c(0.5),
   scale=c(4),
   shape=c(0.125),
   component.type=c("common","specific"))

input.onesegment <- data.table(vals = rnorm(2000, mean = 5.0, sd = 0.75), seg = 1)
	
input.onecolumn <- data.table(vals = rnorm(2000, mean = 5.0, sd = 0.75))

input.threecolumn <- data.table(vals = rnorm(2000, mean = 5.0, sd = 0.75), seg = 1, extra = "hi")

test_that("Input data is ok", {
	expect_error(fit_comixture(input.onesegment, parameters.ok, algorithm = "VB"), "At least two segments required", fixed=TRUE)
	expect_error(fit_comixture(input.onecolumn, parameters.ok, algorithm = "VB"), "Two columns of input required", fixed=TRUE)
	expect_error(fit_comixture(input.threecolumn, parameters.ok, algorithm = "VB"), "Two columns of input required", fixed=TRUE)
})


parameters.nospecific <- data.table(
   mean=c(6,4),
   nu=c(0.5),
   scale=c(4),
   shape=c(0.125),
   component.type=c("common","common"))

parameters.nocommon <- data.table(
   mean=c(6,4),
   nu=c(0.5),
   scale=c(4),
   shape=c(0.125),
   component.type=c("specific","specific"))

parameters.nocomponenttype <- data.table(
	mean=c(6,4),
	nu=c(0.5),
	scale=c(4),
	shape=c(0.125))

parameters.nomean <- data.table(
	nu=c(0.5),
	scale=c(4),
	shape=c(0.125),
	component.type=c("specific","common"))

parameters.noscale <- data.table(
	mean=c(6,4),
	nu=c(0.5),
	shape=c(0.125),
	component.type=c("specific","common"))

parameters.noshape <- data.table(
	mean=c(6,4),
	nu=c(0.5),
	scale=c(0.125),
	component.type=c("specific","common"))

parameters.nonu <- data.table(
	mean=c(6,4),
	shape=c(0.5),
	scale=c(0.125),
	component.type=c("specific","common"))

parameters_no_w <- data.table(
	mean = c(6,4),
	variance = 1,
	component.type = c("common","specific"))

parameters_no_variance <- data.table(
	w = 0.5,
	mean = c(6,4),
	component.type = c("common","specific"))


test_that("Input parameters are ok", {
	expect_error(fit_comixture(input.ok, parameters.nospecific, algorithm = "VB"), "At least one common and one specific component required", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.nocommon, algorithm = "VB"), "At least one common and one specific component required", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.nocomponenttype, algorithm = "VB"), "At least one common and one specific component required", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.nomean, algorithm = "VB"), "Mean values required", fixed=TRUE)
	
	# VB specific
	expect_error(fit_comixture(input.ok, parameters.nonu, algorithm = "VB"), "Nu values required", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.noshape, algorithm = "VB"), "Shape parameter values required", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.noscale, algorithm = "VB"), "Scale parameter values required", fixed=TRUE)

	# EM specific
	expect_error(fit_comixture(input.ok, parameters_no_w, algorithm = "EM"), "Component weighting values required", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters_no_variance, algorithm = "EM"), "Variance values required", fixed=TRUE)
})

test_that("required function options are ok", {
	expect_error(fit_comixture(input.ok, parameters.ok, algorithm = "NO"), "Algorithm should be either EM or VB", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.ok, algorithm = "VB", rho = "hi"), "rho should be a single numeric value", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.ok, algorithm = "VB", rho = 2), "rho should be between 0 and 1", fixed=TRUE)
	expect_error(fit_comixture(input.ok, parameters.ok, algorithm = "VB", rho = c(0.5,0.5)), "rho should be a single numeric value", fixed=TRUE)
})

