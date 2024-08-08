## Introduction
This project estimates the GARCH(1,1) model using Markov Chain Monte Carlo or Gaussian Variational Bayes (GVB) to analyze the volatility in the stock market. In this project, we work with the dataset of the daily return values of the Australia All Ordinaries Index (AORD) from 03-Jan-2012 to 28-June-2022. 

GARCH model is used to analyse the volatility of financial time series data. The model takes in values of the past squared observations and volatility to model the variance at time 𝑡. The GARCH (1,1) model is favoured for its relatively simple implementation which forecasts volatility by fitting one autoregressive lag or ARCH term and one moving average lag (GARCH term).

$$ y_{t} = \sigma_{t} \epsilon_{t},\ \ \  \epsilon_{t} \sim N(0,1),\ \ \ t = 1,2,...,T $$
$$ \sigma_{t}^2 = \omega + \alpha y_{t-1}^2 + \beta\sigma_{t-1}^2,\ \ \ t = 2,...,T $$

I was mainly responsible for estimating the GARCH model by implementing Gaussian Variational Bayes and generating the forecast of the AORD volatility for 29-June-2022 with the estimated model. 

