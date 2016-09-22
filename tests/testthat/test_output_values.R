
context("Reasonable output values")

# single common and specific component

input <- rbind(data.table(vals = rnorm(2000, mean = 5.0, sd = 0.75), seg = 1),
               data.table(vals = rnorm(1000, mean = 0.5, sd = 1.0), seg = 1),
              data.table(vals = rnorm(2000, mean = 5.5, sd = 0.75), seg = 2),
              data.table(vals = rnorm(2000, mean = 4.0, sd = 3.0), seg = 2))

initial.parameters.vb <- data.table(
   mean=c(6,4),
   nu=c(0.5),
   scale=c(4),
   shape=c(0.125),
   component.type=c("common","specific"))

initial.parameters.em <- data.table(
	w = 0.5,
	mean = c(6,4),
	variance = 1,
	component.type = c("common","specific"))

output <- fit.model(input, initial.parameters.vb, quiet = T, algorithm = "VB")

output.em <- fit.model(input, initial.parameters.em, quiet = T, algorithm = "EM")

test_that("common parameter estimates are reasonable", {
	expect_equal(output$common_parameters$mean, 5, tolerance = 0.3)
	expect_equal(output.em$common_parameters$mean, 5, tolerance = 0.3)
	
	expect_equal(output$common_parameters$variance, 0.75^2, tolerance = 0.1)
	expect_equal(output.em$common_parameters$variance, 0.75^2, tolerance = 0.1)
	
	expect_equal(output$common_parameters$mix_weights, 1, tolerance = 0.05)
	expect_equal(output.em$common_parameters$mix_weights, 1, tolerance = 0.05)
	
	expect_equal(output$rho, c("1" = 1/3, "2" = 1/2), tolerance = 0.05)
	expect_equal(output.em$rho, c(1/3, 1/2), tolerance = 0.05)
})

test_that("specific parameter estimates are reasonable", {
	expect_equal(output$specific_parameters$mean, matrix(c(0.5,4), ncol = 1), tolerance = 0.3)
	expect_equal(output.em$specific_parameters$mean, matrix(c(0.5,4), ncol = 1, dimnames = list(c("1",2),NULL)), tolerance = 0.3)
	
	expect_equal(output$specific_parameters$variance, matrix(c(1,9), ncol = 1), tolerance = c(0.1,0.4))
	expect_equal(output.em$specific_parameters$variance, matrix(c(1,9), ncol = 1, dimnames = list(c("1",2),NULL)), tolerance = c(0.1,1.1))
	
	expect_equal(output$specific_parameters$mix_weights, matrix(c(1,1), ncol = 1), tolerance = 0.05)
	expect_equal(output.em$specific_parameters$mix_weights, matrix(c(1,1), ncol = 1, dimnames = list(c("1",2),NULL)), tolerance = 0.05)
})

# multuple common and specific components

test.data.basic <- rbind(
	# common components
	data.table(vals = rnorm(2000, mean = 7, sd = 0.5), seg = 1),
	data.table(vals = rnorm(2000, mean = 7, sd = 0.5), seg = 2),
	data.table(vals = rnorm(1500, mean = 11, sd = 0.5), seg = 1),
	data.table(vals = rnorm(1500, mean = 11, sd = 0.5), seg = 2),
	# unique components
	data.table(vals = rnorm(3000, mean = 1, sd = 0.5), seg = 1),
	data.table(vals = rnorm(2000, mean = 3, sd = 0.5), seg = 2),
	data.table(vals = rnorm(1500, mean = 16, sd = 0.5), seg = 1),
	data.table(vals = rnorm(1500, mean = 14, sd = 0.5), seg = 2)
	)

initial.parameters.basic.vb <- data.table(
  mean = c(5,10,2,15),
  nu = c(5),
  scale = c(4),
  shape = c(0.125),
  component.type = c("common","common","specific","specific"))

initial.parameters.basic <- data.table(
  w=c(0.5),
  mean=c(3.5,10,2,12),
  variance=c(0.6^2),
  component.type=c("common","common","specific","specific"))

output.2.em <- fit.model(test.data.basic, initial.parameters.basic, quiet = T, algorithm = "EM")

output.2 <- fit.model(test.data.basic, initial.parameters.basic.vb, max.iterations = 30, quiet = T, algorithm = "VB")

test_that("common parameter estimates are reasonable", {
	expect_equal(output.2$common_parameters$mean, c(7,11), tolerance = 0.1)
  expect_equal(output.2.em$common_parameters$mean, c(7,11), tolerance = 0.1)
  
	expect_equal(output.2$common_parameters$variance, c(0.25,0.25), tolerance = 0.05)
	expect_equal(output.2.em$common_parameters$variance, c(0.25,0.25), tolerance = 0.05)
	
	expect_equal(output.2$common_parameters$mix_weights, c(2/3.5,1.5/3.5), tolerance = 0.05)
	expect_equal(output.2.em$common_parameters$mix_weights, c(2/3.5,1.5/3.5), tolerance = 0.05)
	
	expect_equal(output.2$rho, c("1" = 4.5/(3.5+4.5), "2" = 3.5/(3.5+3.5)), tolerance = 0.05)
	expect_equal(output.2.em$rho, c(4.5/(3.5+4.5), 3.5/(3.5+3.5)), tolerance = 0.05)
})

test_that("specific parameter estimates are reasonable", {
	expect_equal(output.2$specific_parameters$mean, matrix(c(1,3,16,14), ncol = 2), tolerance = 0.1)
  expect_equal(output.2.em$specific_parameters$mean, matrix(c(1,3,16,14), ncol = 2, dimnames = list(c("1",2),NULL)), tolerance = 0.1)
  
	expect_equal(output.2$specific_parameters$variance, matrix(c(0.25,0.25,0.25,0.25), ncol = 2), tolerance = 0.05)
	expect_equal(output.2.em$specific_parameters$variance, matrix(c(0.25,0.25,0.25,0.25), ncol = 2, dimnames = list(c("1",2),NULL)), tolerance = 0.05)
	
	expect_equal(output.2$specific_parameters$mix_weights, matrix(c(3/4.5,2/3.5,1.5/4.5,1.5/3.5), ncol = 2), tolerance = 0.05)
	expect_equal(output.2.em$specific_parameters$mix_weights, matrix(c(3/4.5,2/3.5,1.5/4.5,1.5/3.5), ncol = 2, dimnames = list(c("1",2),NULL)), tolerance = 0.05)
})

