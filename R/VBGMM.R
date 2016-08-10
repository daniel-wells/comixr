y <- c(rnorm(5000, mean = 5.0, sd = 0.75),rnorm(2500, mean = 0.5, sd = 1.0))
hist(y,breaks=100)

n.comp <- 2

##### reponsibilities
gamma <- matrix(nrow=length(y),ncol=n.comp)

##### parameters
# mixing weight
pi <- matrix(rep((1/n.comp),n.comp),nrow=1)

# precision
beta <- matrix(rep(0.5,n.comp),nrow=1)

# mean
mu <- matrix(c(1,2),nrow=1)
  

##### hyperparameters
# mixing prior
lambda <- 100*pi
lambda_0 <- 5

# means prior
m <- mu
m_0 <- mean(y)
nu <- (1/beta)/(length(y)*pi)
nu_0 <- (sum(abs(range(y)))/3)^2

# precisions prior
b <- matrix(rep(0.7,n.comp),nrow=1)
b_0 <- 1000
c <- matrix(rep(0.7,n.comp),nrow=1)
c_0 <- 0.001

##### E step
pi = exp( digamma(lambda) - digamma(sum(lambda)) )
beta = b*c

gamma <- matrix(pi*beta^0.5,nrow=length(y),ncol=n.comp) * exp(-0.5 * matrix(beta,nrow=length(y),ncol=n.comp) * ( outer(y^2, as.vector(m^2 + nu), FUN="+") - 2*y%*%m ) )
gamma <- gamma/rowSums(gamma)

#### M Step
# number of data points in each component
N_s <- colSums(gamma)
stopifnot(sum(N_s)==length(y))

# proportion of data in each component
pi <- N_s/length(y)
stopifnot(sum(pi)==1.0)

y_w <- colSums(gamma*y)/length(y)

y2_w <- colSums(gamma*y^2)/length(y)

lambda <- N_s + lambda_0

sigma2 <- y2_w + pi*(m^2+nu) - 2*m*y_w

b <- 1/( length(y)/2 * sigma2 + 1/b_0 )

c <- N_s/2 + c_0

m_data <- y_w/pi

tau_data <- N_s*beta

tau <- 1/nu_0 + tau_data

nu <- 1/tau

m <- nu_0/tau * m_0 + tau_data/tau * m_data

n <- length(y)
y_infered <- c(rnorm(as.integer(n*pi[1]), mean = m[1], sd = 1/beta[1]),
              rnorm(as.integer(n*pi[2]), mean = m[2], sd = 1/beta[2]))

toplot <- rbind(data.frame(y=y,source="original"),data.frame(y=y_infered,source="Infered"))
ggplot(toplot,aes(y,colour=source))+ geom_freqpoly(bins=50)
