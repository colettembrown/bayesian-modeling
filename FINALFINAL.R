# Load data
setwd("/Users/Colette/anaconda/1CS 146/FINAL/")
raw.data <- read.csv('co2.csv', header = FALSE)

df <- subset(raw.data, select = c('V1', 'V2', "V3", "V4"))
colnames(df) <- c('date', 'ppm', 'raw_date', "scaled_ppm")

# removing weird NA's
df <- df[-c(1, 3041, 3040), ]

plot(df$date, df$ppm, type = "l", 
     xlab="Time", ylab="PPM",
     main = "Mauna Loa CO2 (PPM)")

t_future <- seq(max(df$date), (max(df$date) + (2058-2017)*365), by = 7)
predict <- data.frame(t_future)

N = 3038              # 2038 data points for training
n_future = 2138
t = df$date        # dates/timesteps
t_future = predict$t_future
ppm = df$ppm     # ppm

stan_data <- list(
  N = N,               # 3040 data points
  t = t,               # dates/timesteps
  ppm = ppm,           # ppm
  n_future = n_future,         # number of weeks
  t_future = t_future
)

stan_model = "
data {
int <lower = 0> N;   // the number of timesteps in the data 
int <lower = 0> n_future; // number of future timesteps
vector [N] t;      // timesteps 
vector [n_future] t_future;      // timesteps 
vector [N] ppm;    // ppm values
}

parameters {
real c0;
real c1;
real c2;
real c3;
real c4;
real c5;
}

transformed parameters {
real c0_transf;
real c1_transf;
real c2_transf;
real c3_transf;
real c4_transf;
real c5_transf;

c0_transf = exp(c0);
c1_transf = exp(c1);
c2_transf = exp(c2); // keep this for sure (mapping real numbers to positive real)
// phi_transf = (atan2(phi_1, phi_2)); // phase between 0 and 2pi (think logistic function (2pi()/1+e^-x))
c3_transf = (c3);
c4_transf = exp(c4);
c5_transf = exp(c5);
}

model {
// Priors
c0 ~ normal(0, 1);
c1 ~ normal(0,1);
c2 ~ normal(0,1);
c3 ~ normal(0,1);
c4 ~ normal(0,1);
c5 ~ normal(0,1);

// Likelihood 
for (n in 1:N)
  ppm[n] ~ normal(c0_transf + c1_transf * t[n] + (c5_transf * t[n] * t[n]) +c2_transf*cos(2 * pi() * t[n] / 365.25 + c3_transf), c4_transf); 
}

generated quantities {

  real ppm_future[n_future];
  for (x in 1:(n_future)){
    ppm_future[x] = normal_rng(c0_transf + c1_transf * t_future[x] + (c5_transf * t_future[x] * t_future[x]) + c2_transf*cos(2 * pi() * t_future[x] / 365.25 + c3_transf), c4_transf);
  }
}
"

library(rstan)

fit <- stan(
  model_code = stan_model,
  data = stan_data,       # named list of data
  chains = 4,             # number of Markov chain
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 4,              # number of cores 
  refresh = 1000,         # show progress every 'refresh' iterations
  control = list(adapt_delta = 0.9)
)

samples <- extract(fit)
print(fit, par=c('c0_transf','c1_transf','c2_transf','c3_transf','c4_transf', 'c5_transf'), probs = c(0.05, 0.5, 0.95))

results <- apply(samples$ppm_future, 2, quantile, probs = c(0.025, 0.5, 0.975))

print(all_time[min(which(results[2,1:n_future] * sd(data$CO2) + mean(data$CO2) > 450))])

plot(
  df$date, df$ppm,
  col='black', type = 'l',
  xlim=c(0, t_future[n_future]),
  ylim=c(min(c(results[1,], ppm)), max(c(results[2,], ppm))),
  xlab="Date", ylab="PPM",
  main='Data, future data, and predicted 95% interval')
lines(c(predict$t_future), c(results[1,]), col='blue')
lines(c(predict$t_future), c(results[2,]), col='black')
lines(c(predict$t_future), c(results[3,]), col='blue')
abline(v=t[N], col='red')

plot(df$date,df$ppm,xlim=c(20100,23500),ylim=c(380,430),
     col='black',type="l",
     xlab="Date", ylab="PPM",
     main='Data, future data, and predicted 95% interval')
lines(c(predict$t_future), c(results[1,]), col='blue')
lines(c(predict$t_future), c(results[2,]), col='black')
lines(c(predict$t_future), c(results[3,]), col='blue')
abline(v=t[N], col='red')

projected_ppm <- max(results[2,])
print(cat("CO2 ppm is projected to be", projected_ppm, "in the year 2058"))

plot(df$date,df$ppm,xlim=c(20100, 23500),ylim=c(450,452))
abline(v=t[N])

