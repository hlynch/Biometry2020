Week 2 Lab
=============

Confidence intervals
-----------------------

Before getting too far, we need to circle back and make sure we understand what is meant by a confidence interval. 

A 95th percentile confidence interval say “If I repeat this procedure 100 times using 100 different datasets, 95% of the time my confidence intervals will capture the true parameter”. It does NOT say that there is a 95% chance that the parameter is in the interval.

**Quiz time! (Don't worry, not a real quiz)**

*Important note*: This is an area where Aho is WRONG. I will not repeat Aho's interpretation here because I think he's just wrong. Aho is correct on only one point. It is true that ONCE THE 95th CI HAS BEEN CONSTRUCTED, it is no longer possible to assign a $%$ to the probability that that CI contains the true value or not. Because that CI, once created, either DOES or DOES NOT contain the true value. However, we often talk about the interval in the abstract. When we say "There is a 95$%$ chance that the interval contains the true value" what we mean is that there is a 95$%$ probability that a CI created using that methodology would contain the true value.

Do not let Week 2 pass by without fundamentally understanding the interpretation of a confidence interval. 

Testing hypotheses through permutation
------------------------------------

We'll start off by working through two examples from Phillip Good's book "Introduction to Statistics Through Resampling Methods and R/S-PLUS":

**Example #1**: Use permutation methods to test the null hypothesis that the treatment does not increase survival time (in other word: $H_{0}$: No difference in survival between the treated and control groups):

Survival.treated=$\{94,197,16,38,99,141,23 \}$

Survival.control=$\{52,104,146,10,51,30,40,27,46 \}$

(Is this a one-tailed or a two-tailed test?)

Make sure that you understand what is being done here, as this example is very closely related to the problem set.


**Example #2**: Using the same data, provide a 95% confidence interval for the difference in mean survival days based on 1000 bootstrap samples

Note that these two approaches are very closely related. Do you see why either approach can be used to test the null hypothesis? (What is the null hypothesis here?)

Now we will do one slightly more complicated example from Phillip Good's book "Permutation tests: A practical guide to resampling methods and testing hypotheses":

Holmes and Williams (1954) studied tonsil size in children to verify a possible association with the virus \textit{S. pyrogenes}. Test for an association between \textit{S. pyrogenes} status and tonsil size. (Note that you will need to come up with a reasonable test statistic.)

<div class="figure" style="text-align: center">
<img src="Table2categories.png" alt="Data on tonsil size and S. pyrogenes status. Source: Good (1994)" width="40%" />
<p class="caption">(\#fig:unnamed-chunk-1)Data on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

Now lets consider the full dataset, where tonsil size is divided into three categories. How would we do the test now? What is the new test statistic? (There are many options.) What 'labels' do you permute?

<div class="figure" style="text-align: center">
<img src="Table3categories.png" alt="Fill dataset on tonsil size and S. pyrogenes status. Source: Good (1994)" width="50%" />
<p class="caption">(\#fig:unnamed-chunk-2)Fill dataset on tonsil size and S. pyrogenes status. Source: Good (1994)</p>
</div>

Basics of bootstrap and jackknife
------------------------------------

To get started with bootstrap and jackknife techniques, we start by working through a very simple example. First we simulate some data


```r
x<-seq(0,9,by=1)
```

This will constutute our "data". Let's print the result of sampling with replacement to get a sense for it...


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 3 4 5 7 8 
## 1 1 2 3 1 1 1
```

Now we will write a little script to take bootstrap samples and calculate the means of each of these bootstrap samples


```r
xmeans<-vector(length=1000)
for (i in 1:1000)
  {
  xmeans[i]<-mean(sample(x,replace=T))
  }
```

The actual number of bootstrapped samples is arbitrary *at this point* but there are ways of characterizing the precision of the bootstrap (jackknife-after-bootstrap) which might inform the number of bootstrap samples needed. *In practice*, people tend to pick some arbitrary but large number of bootstrap samples because computers are so fast that it is often easy to draw far more samples than are actually needed. When calculation of the statistic is slow (as might be the case if you are using the samples to construct a phylogeny, for example), then you would need to be more concerned with the number of bootstrap samples. 

First, lets just look at a histogram of the bootstrapped means and plot the actual sample mean on the histogram for comparison



```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=2)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-6-1.png" width="672" />

Calculating bias and standard error
-----------------------------------

From these we can calculate the bias and standard deviation for the mean (which is the "statistic"):

$$
\widehat{Bias_{boot}} = \left(\frac{1}{k}\sum^{k}_{i=1}\theta^{*}_{i}\right)-\hat{\theta}
$$


```r
bias.boot<-mean(xmeans)-mean(x)
bias.boot
```

```
## [1] -0.0443
```

```r
hist(xmeans,breaks=30,col="pink")
abline(v=mean(x),lwd=5,col="black")
abline(v=mean(xmeans),lwd=2,col="yellow")
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-7-1.png" width="672" />

$$
\widehat{s.e._{boot}} = \sqrt{\frac{1}{k-1}\sum^{k}_{i=1}(\theta^{*}_{i}-\bar{\theta^{*}})^{2}}
$$


```r
se.boot<-sd(xmeans)
```

We can find the confidence intervals in two ways:

Method #1: Assume the bootstrap statistics are normally distributed


```r
LL.boot<-mean(xmeans)-1.96*se.boot #where did 1.96 come from?
UL.boot<-mean(xmeans)+1.96*se.boot
LL.boot
```

```
## [1] 2.726867
```

```r
UL.boot
```

```
## [1] 6.184533
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7000 6.1025
```

Let's compare this to what we would have gotten if we had used normal distribution theory. First we have to calculate the standard error:


```r
se.normal<-sqrt(var(x)/length(x))
LL.normal<-mean(x)-qt(0.975,length(x)-1)*se.normal
UL.normal<-mean(x)+qt(0.975,length(x)-1)*se.normal
LL.normal
```

```
## [1] 2.334149
```

```r
UL.normal
```

```
## [1] 6.665851
```

In this case, the confidence intervals we got from the normal distribution theory are too wide.

Does it make sense why the normal distribution theory intervals are too wide? Because the original were were uniformly distributed, the data has higher variance than would be expected and therefore the standard error is higher than would be expected.

There are two packages that provide functions for bootstrapping, 'boot' and 'boostrap'. We will start by using the 'bootstrap' package, which was originally designed for Efron and Tibshirani's monograph on the bootstrap. 

To test the main functionality of the 'bootstrap' package, we will use the data we already have. The 'bootstrap' function requires the input of a user-defined function to calculate the statistic of interest. Here I will write a function that calculates the mean of the input values.


```r
library(bootstrap)
theta<-function(x)
  {
    mean(x)
  }
results<-bootstrap(x=x,nboot=1000,theta=theta)
results
```

```
## $thetastar
##    [1] 2.1 3.7 4.3 5.4 7.0 4.9 3.7 6.0 3.3 5.0 4.0 5.4 5.6 3.4 5.3 3.0 5.5 4.4
##   [19] 4.8 7.0 4.4 5.4 3.7 2.9 4.3 4.8 4.3 4.6 4.7 4.6 5.9 4.6 3.2 2.7 3.3 4.3
##   [37] 5.0 5.9 4.8 4.4 4.8 3.6 3.1 6.0 3.2 5.2 3.3 5.4 5.3 4.6 3.2 3.0 2.8 5.2
##   [55] 4.4 4.8 3.4 4.5 5.1 4.3 6.8 5.1 5.4 4.3 6.8 4.3 4.5 4.3 3.9 4.7 6.0 3.9
##   [73] 4.6 1.1 4.5 4.6 4.2 5.3 5.7 5.3 4.0 3.1 5.1 4.0 4.2 4.7 3.4 1.9 4.0 3.8
##   [91] 6.2 2.7 3.5 6.5 5.0 5.3 4.0 4.2 4.7 5.5 5.4 4.5 5.6 4.4 4.9 5.7 3.9 3.6
##  [109] 5.2 4.8 4.8 5.0 5.9 5.0 5.2 3.8 3.9 5.6 5.5 4.6 3.7 5.5 4.2 4.5 5.9 3.7
##  [127] 5.5 1.7 5.2 4.8 4.7 4.5 4.6 4.4 4.0 3.1 4.3 4.1 4.1 4.5 5.3 4.4 4.7 3.8
##  [145] 3.3 5.5 3.2 5.6 4.7 3.5 4.0 5.5 5.0 5.6 3.3 4.9 3.9 4.9 5.0 3.9 3.0 4.4
##  [163] 6.4 4.5 4.5 5.2 5.1 4.3 5.8 4.0 4.6 4.3 4.8 2.7 3.6 5.1 4.3 4.3 3.9 4.5
##  [181] 3.9 5.0 5.0 2.5 4.1 4.4 4.0 5.8 3.3 4.3 4.9 4.2 3.9 4.7 4.1 2.7 6.4 4.5
##  [199] 4.7 5.0 3.1 5.2 3.7 3.6 2.7 5.0 4.8 4.9 5.5 3.9 5.2 4.2 3.5 5.0 4.6 4.0
##  [217] 3.5 5.6 4.2 4.4 4.7 4.7 3.4 4.5 5.8 4.5 3.4 4.7 5.4 4.0 5.4 3.6 2.3 3.6
##  [235] 5.1 4.8 3.6 5.3 4.1 3.8 4.4 4.7 4.9 4.5 4.2 4.5 5.4 4.2 4.9 3.6 5.5 4.7
##  [253] 4.6 4.1 4.1 5.0 4.3 4.3 4.5 6.1 5.2 3.7 2.9 4.2 3.8 4.0 6.6 4.6 4.9 4.9
##  [271] 5.8 5.9 4.9 3.1 4.3 2.4 4.5 3.5 3.6 3.8 4.8 4.4 3.5 4.4 4.0 6.1 6.0 5.7
##  [289] 3.3 4.3 4.0 4.0 4.7 3.9 4.9 4.7 4.6 3.3 5.9 3.9 4.4 4.6 6.4 4.2 5.8 2.9
##  [307] 5.4 5.1 5.1 4.2 3.6 4.7 5.5 3.4 4.0 4.2 5.1 5.5 5.7 5.6 3.9 4.9 4.0 5.2
##  [325] 4.7 2.9 4.4 4.6 4.9 4.6 4.3 5.2 2.9 3.1 4.9 3.7 6.8 6.4 5.4 6.4 5.0 4.2
##  [343] 4.0 1.5 4.5 4.9 3.7 5.4 4.8 4.1 4.7 3.0 5.4 3.3 1.6 4.4 4.8 5.9 4.8 5.7
##  [361] 5.3 4.8 4.6 5.5 4.1 4.3 3.9 4.8 6.6 4.2 5.8 3.6 4.0 4.2 4.0 4.6 4.8 4.3
##  [379] 3.0 5.0 6.8 3.8 5.2 5.4 5.7 3.6 6.1 3.3 4.0 5.5 6.1 3.6 4.6 2.7 4.4 6.0
##  [397] 3.2 4.0 4.3 5.0 4.2 4.3 4.0 5.4 5.0 4.2 4.4 5.3 6.3 5.8 4.9 6.1 4.6 2.1
##  [415] 3.2 2.7 4.7 4.3 6.4 5.5 4.3 5.0 4.7 4.2 3.8 3.4 5.8 3.9 2.5 3.8 3.4 5.5
##  [433] 4.8 3.7 3.0 3.5 5.2 6.4 5.3 4.2 4.1 4.9 2.9 4.0 7.2 4.7 5.3 4.2 4.9 4.2
##  [451] 4.2 2.2 2.8 4.3 5.1 2.3 4.8 6.0 4.8 4.0 5.8 3.9 4.7 4.8 4.0 5.1 2.7 5.3
##  [469] 4.2 4.3 5.2 3.6 5.3 3.8 4.2 5.0 5.2 3.8 5.3 4.7 5.1 5.7 4.4 6.1 4.1 5.5
##  [487] 4.3 5.1 3.9 4.3 5.1 3.3 4.2 2.1 6.1 3.7 4.9 5.3 4.1 3.9 7.2 4.5 4.8 4.9
##  [505] 4.6 5.9 3.1 6.1 4.0 4.7 3.8 4.2 4.6 3.8 3.7 4.7 5.3 4.1 3.7 4.7 3.4 4.2
##  [523] 3.4 5.4 4.4 3.9 5.6 3.6 3.5 3.3 6.8 5.0 3.5 3.0 4.4 4.9 3.8 5.1 3.9 5.9
##  [541] 4.1 3.8 2.9 3.8 4.8 3.6 1.5 5.3 5.7 5.1 5.7 5.0 4.6 6.4 4.1 4.7 3.7 4.7
##  [559] 5.3 4.9 4.3 3.0 3.2 5.6 4.1 6.4 5.5 3.7 4.7 4.9 3.2 6.1 4.7 4.3 6.4 5.3
##  [577] 5.3 5.4 5.6 5.5 4.1 3.7 4.6 5.9 4.0 3.7 5.2 4.5 4.9 5.1 4.3 5.5 5.0 4.7
##  [595] 3.5 4.9 5.0 5.9 4.7 3.9 3.3 3.7 4.3 5.9 4.6 4.4 5.2 5.2 4.5 5.5 4.8 3.9
##  [613] 3.9 4.1 2.6 3.6 4.4 7.4 3.9 3.7 3.4 4.3 4.5 4.4 5.4 4.8 4.3 6.0 5.2 2.6
##  [631] 2.1 4.2 5.7 4.0 5.9 4.7 5.3 4.0 4.3 5.1 4.6 4.4 4.7 5.7 2.9 4.2 4.8 4.3
##  [649] 4.9 4.5 6.0 2.9 5.0 2.8 5.0 6.3 4.9 2.1 4.3 3.2 4.3 4.6 3.6 4.0 4.4 4.7
##  [667] 5.7 4.5 3.8 5.3 3.4 4.5 4.5 5.3 3.0 5.7 4.6 3.8 5.1 4.4 5.5 5.7 3.9 4.9
##  [685] 4.6 2.8 3.2 4.2 4.7 2.8 4.4 3.5 4.0 4.0 5.0 4.1 4.6 4.9 2.9 4.7 5.2 4.8
##  [703] 5.9 3.9 4.3 5.2 4.0 3.8 4.2 4.2 4.8 6.6 4.0 5.5 3.8 4.4 3.2 4.3 3.5 3.7
##  [721] 4.0 3.6 6.3 4.8 5.3 5.3 3.6 4.1 4.1 4.7 4.3 3.6 4.5 4.0 4.8 4.4 4.9 5.2
##  [739] 4.7 6.9 3.6 4.8 4.8 3.8 5.1 5.4 4.3 5.8 5.4 5.0 5.6 5.0 4.3 4.1 4.2 4.7
##  [757] 3.7 5.5 5.6 4.2 4.1 2.8 4.2 4.2 5.8 3.5 5.0 3.8 4.6 4.4 4.6 4.6 5.1 6.4
##  [775] 4.5 3.1 3.4 3.2 3.2 4.4 6.4 2.6 3.5 4.6 4.7 5.2 2.3 4.2 3.7 3.9 4.7 3.4
##  [793] 4.1 5.1 5.4 4.7 5.5 6.5 4.6 5.2 3.2 3.8 4.0 3.8 4.2 6.3 3.7 2.7 5.1 3.9
##  [811] 4.4 3.6 4.0 5.9 4.5 2.9 4.5 5.7 5.1 4.5 5.1 5.8 5.5 5.1 5.0 4.6 4.1 4.4
##  [829] 5.0 4.5 5.4 3.2 4.0 2.6 3.3 4.1 3.7 4.8 4.7 5.5 3.8 5.0 4.4 4.7 3.8 5.8
##  [847] 4.2 2.9 4.7 4.0 5.2 4.6 3.5 5.8 3.9 6.0 3.9 5.3 5.5 2.8 4.7 4.8 3.9 1.1
##  [865] 6.6 4.9 5.8 5.0 2.6 4.0 6.1 3.7 2.4 5.0 6.3 6.7 4.1 5.0 5.4 4.3 4.6 3.7
##  [883] 4.4 4.3 2.2 4.9 4.5 3.5 4.3 4.5 3.5 5.5 4.0 4.4 5.4 3.9 5.2 3.8 4.9 5.9
##  [901] 3.6 4.4 4.7 4.3 5.0 4.0 4.9 3.7 5.1 5.4 3.9 5.7 5.0 3.9 5.1 2.9 5.4 5.4
##  [919] 5.4 3.8 5.4 3.5 4.9 5.2 2.8 3.9 3.5 5.5 5.2 3.4 4.3 3.8 4.0 4.3 5.3 4.2
##  [937] 4.9 5.9 2.4 3.6 4.9 5.3 4.3 5.6 4.4 5.2 5.1 5.1 3.9 4.5 5.2 3.0 4.7 5.9
##  [955] 3.9 4.9 2.9 3.2 4.9 3.2 5.2 4.4 3.9 5.6 3.7 3.3 6.5 4.5 4.6 5.8 3.5 5.3
##  [973] 4.8 3.8 5.1 3.0 5.3 4.5 4.8 3.5 4.6 4.2 2.9 5.8 3.6 5.0 4.2 3.2 3.2 6.6
##  [991] 4.2 3.0 4.3 3.7 3.8 4.1 4.7 3.5 3.0 4.1
## 
## $func.thetastar
## NULL
## 
## $jack.boot.val
## NULL
## 
## $jack.boot.se
## NULL
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta)
```

```r
quantile(results$thetastar,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.6   6.4
```

Notice that we get exactly what we got last time. This illustrates an important point, which is that the bootstrap functions are often no easier to use than something you could write yourself.

You can also define a function of the bootstrapped statistics (we have been calling this theta) to pull out immediately any summary statistics you are interested in from the bootstrapped thetas.

Here I will write a function that calculates the bias of my estimate of the mean (which is 4.5 [i.e. the mean of the number 0,1,2,3,4,5,6,7,8,9])


```r
bias<-function(x)
  {
  mean(x)-4.5
  }
results<-bootstrap(x=x,nboot=1000,theta=theta,func=bias)
results
```

```
## $thetastar
##    [1] 5.2 4.1 3.3 4.2 4.9 3.8 3.9 5.1 4.2 4.0 3.5 4.3 4.8 4.1 4.9 5.8 5.4 4.9
##   [19] 4.2 5.5 5.5 2.6 3.7 4.7 5.0 4.5 4.1 6.0 4.7 4.8 5.5 4.3 5.0 5.6 3.9 4.4
##   [37] 3.9 3.7 4.8 4.3 4.5 3.9 5.4 4.1 3.7 6.3 3.0 5.5 2.5 3.9 3.9 4.3 3.0 4.3
##   [55] 4.3 5.1 4.5 4.7 3.7 3.9 4.2 3.7 3.3 4.6 4.8 5.6 2.9 4.7 3.5 4.9 5.4 4.6
##   [73] 4.5 4.6 7.1 3.5 6.0 3.0 4.0 3.9 4.8 2.3 5.5 5.3 3.4 5.9 5.3 5.3 3.1 4.2
##   [91] 5.8 2.7 3.7 4.3 4.1 2.9 2.9 3.9 5.2 6.1 4.1 5.5 4.8 4.5 4.3 3.5 4.9 6.3
##  [109] 5.2 4.5 4.8 4.5 4.4 5.4 3.8 4.2 5.0 5.5 4.7 4.2 4.5 4.5 3.4 5.3 2.5 4.7
##  [127] 4.5 4.4 5.5 4.4 5.3 5.2 5.4 4.5 3.5 3.3 5.5 5.4 4.9 3.5 4.9 4.3 3.9 4.8
##  [145] 5.8 3.6 4.3 5.2 4.3 5.8 5.4 3.9 4.5 3.8 4.7 3.8 5.4 4.4 4.8 4.8 5.0 5.6
##  [163] 5.2 4.7 5.0 5.2 4.7 3.0 5.1 4.5 3.6 3.2 3.1 3.2 3.9 3.8 5.7 6.2 4.8 5.8
##  [181] 2.6 4.7 6.3 3.6 6.1 4.6 4.9 2.8 4.3 2.4 4.5 4.0 4.9 3.7 4.9 3.9 4.1 6.0
##  [199] 5.2 3.0 4.9 3.9 3.9 4.4 4.6 5.9 5.4 3.5 5.2 3.6 6.2 3.7 4.5 4.8 3.4 4.6
##  [217] 5.1 3.5 5.7 4.3 3.9 3.9 5.9 3.5 5.0 3.6 4.5 5.1 5.2 3.4 5.2 4.0 3.7 5.9
##  [235] 4.2 5.7 4.5 4.8 3.8 2.5 6.3 4.6 4.2 2.6 4.7 5.2 3.6 3.9 5.5 4.1 5.7 5.6
##  [253] 4.1 4.0 5.2 2.6 3.7 3.5 6.0 4.0 3.3 3.2 5.2 4.2 2.5 2.7 5.1 4.6 4.6 4.0
##  [271] 4.9 5.1 2.7 4.0 5.6 5.8 5.3 5.9 4.2 4.0 6.1 4.0 4.7 4.8 5.3 4.7 4.2 4.3
##  [289] 4.8 4.6 5.1 5.2 4.8 3.3 3.3 2.9 4.8 5.6 4.6 4.6 4.9 3.5 4.8 4.3 4.2 4.0
##  [307] 3.9 4.7 5.5 5.4 4.4 5.3 4.1 4.0 6.5 6.0 3.7 4.9 6.1 5.0 5.0 4.1 3.5 6.5
##  [325] 5.9 4.6 4.0 4.9 5.3 3.9 3.3 4.9 4.2 5.6 4.8 3.6 5.2 3.8 3.4 3.3 4.3 5.9
##  [343] 4.1 4.0 5.0 2.7 4.9 3.0 3.3 5.5 6.2 4.0 4.8 4.9 3.9 4.2 5.5 3.8 3.9 4.7
##  [361] 5.7 3.8 6.5 4.0 4.8 3.6 4.9 5.6 4.6 4.4 4.7 3.7 5.0 4.6 4.1 4.9 5.8 4.2
##  [379] 3.8 3.6 4.2 4.8 4.1 4.4 4.4 4.5 4.6 3.2 4.3 3.5 4.6 3.8 4.2 4.8 3.7 4.3
##  [397] 3.7 5.8 4.1 6.0 4.0 5.1 3.9 5.5 5.3 4.8 4.8 3.6 4.8 4.3 4.3 2.9 5.1 4.6
##  [415] 5.8 3.1 4.1 5.6 4.4 3.3 3.8 5.4 4.8 3.7 4.4 6.3 3.3 5.6 3.6 5.0 4.6 5.5
##  [433] 5.2 6.1 5.2 3.5 4.0 4.5 3.8 3.9 5.0 4.2 6.0 4.8 6.0 4.2 2.6 5.4 5.0 3.6
##  [451] 5.1 5.8 2.4 5.5 4.6 4.5 2.8 4.2 5.5 3.7 6.0 2.2 3.9 2.9 4.4 5.2 4.2 2.9
##  [469] 3.3 3.6 4.8 2.7 3.7 4.9 5.7 4.7 4.2 4.5 5.1 4.6 3.4 5.5 4.9 4.8 4.4 7.0
##  [487] 3.3 5.4 4.9 2.5 5.0 5.0 5.2 2.3 5.0 3.8 4.9 4.7 4.6 3.6 4.1 4.9 3.5 3.6
##  [505] 4.7 4.1 4.3 4.6 4.0 4.4 5.2 3.6 5.6 4.9 3.9 4.3 4.3 4.7 2.0 4.4 4.8 3.7
##  [523] 4.1 3.8 4.4 4.9 4.5 4.6 5.7 4.6 5.2 5.1 5.0 3.8 3.8 5.8 3.9 5.4 3.8 5.8
##  [541] 4.1 5.6 3.1 4.3 5.2 5.3 3.4 5.9 3.4 4.2 4.7 5.3 5.9 2.8 5.9 3.9 4.3 5.9
##  [559] 4.9 3.8 5.9 4.2 3.4 4.5 3.8 6.1 6.1 3.7 5.0 5.9 5.6 5.1 6.2 4.7 4.0 4.2
##  [577] 4.0 4.2 4.1 4.2 5.0 3.2 4.1 4.7 5.6 3.5 4.3 4.5 5.6 5.4 5.8 4.0 5.5 4.9
##  [595] 5.0 4.6 3.8 4.3 3.5 4.8 5.3 5.4 4.4 4.7 3.4 4.6 6.0 4.8 5.7 4.2 5.2 4.7
##  [613] 4.3 3.7 5.0 7.7 6.1 4.3 4.9 4.4 4.3 3.6 4.5 4.3 3.3 5.0 3.2 4.2 5.8 4.9
##  [631] 3.7 3.8 6.1 4.3 5.7 5.8 3.6 3.8 4.7 3.4 6.0 4.3 5.9 4.7 3.6 4.6 3.9 5.1
##  [649] 3.5 6.0 3.5 5.3 4.9 4.1 5.1 3.0 3.9 4.4 3.6 4.8 4.7 5.6 4.8 4.8 5.1 4.7
##  [667] 4.9 4.0 2.5 4.2 5.0 5.9 5.2 3.3 3.9 3.9 3.5 5.6 4.7 5.6 2.2 4.0 4.8 4.4
##  [685] 4.6 3.7 4.0 5.2 3.6 4.1 4.5 4.2 4.0 5.0 4.6 3.5 6.0 4.8 4.6 3.6 4.7 3.0
##  [703] 3.8 4.9 5.6 5.5 4.8 5.0 4.8 3.8 3.7 4.4 4.3 5.1 4.6 3.8 4.3 3.8 6.3 3.2
##  [721] 5.3 4.6 3.4 6.4 3.7 3.2 3.6 3.7 4.3 6.1 3.8 3.3 4.9 5.8 3.0 6.0 3.9 6.1
##  [739] 5.3 3.0 6.8 5.0 5.8 6.0 3.7 4.2 3.8 4.6 4.9 5.2 3.4 5.3 4.3 4.1 4.6 3.1
##  [757] 4.8 4.6 4.2 4.1 5.4 4.2 5.1 5.3 5.2 5.4 4.2 3.7 4.5 4.3 5.7 3.8 5.4 4.2
##  [775] 5.1 4.1 4.5 5.2 4.0 6.0 3.7 4.4 2.6 5.0 6.3 4.7 3.1 4.0 2.7 4.8 4.1 1.5
##  [793] 4.7 4.7 4.6 4.2 4.1 4.3 6.0 4.7 4.9 3.9 4.2 5.5 2.8 4.5 5.1 5.2 3.9 6.4
##  [811] 4.5 4.7 4.9 5.2 4.2 3.3 5.8 4.1 4.5 4.9 3.3 5.1 4.1 3.2 6.6 4.7 4.5 5.3
##  [829] 3.9 5.0 5.1 5.0 3.9 3.3 3.3 4.1 5.1 6.0 4.4 3.0 5.2 6.1 4.2 3.0 3.4 6.3
##  [847] 4.7 3.5 5.0 5.6 3.5 3.3 4.9 4.5 3.2 4.9 3.3 5.1 4.6 3.5 3.5 4.1 3.8 5.0
##  [865] 4.7 5.1 2.1 3.8 4.6 5.6 3.1 4.0 6.3 4.4 3.9 5.6 4.3 5.4 4.3 4.1 4.5 4.1
##  [883] 4.1 4.8 5.2 4.8 4.0 3.8 4.5 3.3 3.7 3.4 4.9 5.4 4.5 4.8 3.3 5.2 4.7 5.0
##  [901] 4.7 4.4 3.5 4.9 3.5 3.6 4.5 2.2 3.4 4.0 3.4 5.7 5.5 4.7 5.4 4.4 4.2 4.6
##  [919] 5.0 3.3 4.9 4.2 4.9 3.2 6.2 4.7 3.6 3.9 4.4 3.5 5.5 3.8 3.8 5.9 3.6 4.9
##  [937] 4.9 4.9 4.3 3.9 5.3 4.1 3.9 6.5 3.0 5.3 3.4 3.3 4.4 5.6 3.0 4.1 3.7 5.5
##  [955] 5.5 5.2 4.9 6.4 6.1 3.7 5.3 4.2 5.7 4.5 3.8 6.0 4.4 4.5 3.8 5.1 5.4 4.4
##  [973] 5.6 5.8 4.8 3.4 5.4 5.0 4.2 4.9 6.6 4.1 4.8 3.1 5.4 3.9 4.2 3.9 4.8 5.5
##  [991] 5.4 4.9 3.9 4.5 4.0 4.5 4.7 5.3 4.2 3.9
## 
## $func.thetastar
## [1] 0.003
## 
## $jack.boot.val
##  [1]  0.47391304  0.38689459  0.25314286  0.21913043  0.11265060 -0.04580838
##  [7] -0.08931751 -0.29455587 -0.44985755 -0.44887640
## 
## $jack.boot.se
## [1] 0.9479218
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta, func = bias)
```

Compare this to 'bias.boot' (our result from above). Why might it not be the same? Try running the same section of code several times. See how the value of the bias ($func.thetastar) jumps around? We should not be surprised by this because we can look at the jackknife-after-bootstrap estimate of the standard error of the function (in this case, that function is the bias) and we can see that it is not so small that we wouldn't expect some variation in these values.

Remember, everything we have discussed today are estimates. The statistic as applied to your data will change with new data, as will the standard error, the confidence intervals - everything! All of these values have sampling distributions and are subject to change if you repeated the procedure with new data.

Note that we can calculate any function of $\theta^{*}$. A simple example would be the 72nd percentile:


```r
perc72<-function(x)
  {
  quantile(x,probs=c(0.72))
  }
results<-bootstrap(x=x,nboot=1000,theta=theta,func=perc72)
results
```

```
## $thetastar
##    [1] 6.0 4.3 5.5 5.8 5.4 5.0 5.6 3.4 5.0 3.7 5.0 5.1 4.1 4.6 2.9 4.4 4.2 4.8
##   [19] 4.0 4.2 3.9 5.9 5.0 3.6 4.4 4.5 4.6 3.9 4.3 4.9 5.1 5.4 3.5 3.8 3.9 5.8
##   [37] 4.5 3.4 4.2 6.0 4.5 4.1 3.1 3.3 4.5 5.4 3.5 5.1 3.2 3.7 5.0 5.9 4.9 4.9
##   [55] 3.6 4.7 4.2 3.5 4.3 4.2 5.4 5.4 5.4 4.7 2.5 4.2 5.4 3.9 5.3 4.4 3.1 2.9
##   [73] 4.1 4.2 5.0 3.6 4.4 4.1 6.0 5.6 5.3 4.5 2.3 5.4 3.5 5.7 6.4 4.7 2.1 5.3
##   [91] 4.9 4.2 2.8 5.4 4.4 5.9 5.6 4.6 3.6 5.2 4.4 4.8 3.9 4.6 2.9 3.9 5.1 5.4
##  [109] 3.5 5.4 3.5 4.4 3.2 5.2 4.2 5.9 5.3 3.4 4.6 3.4 4.4 2.9 2.5 4.0 4.9 4.3
##  [127] 4.9 4.3 4.9 3.2 4.5 3.5 3.7 3.3 4.9 4.8 5.1 4.5 5.2 5.1 4.6 3.3 3.7 5.9
##  [145] 4.6 3.4 3.2 5.9 5.9 3.4 3.1 5.4 4.5 5.6 5.1 5.1 3.6 5.2 4.5 4.8 4.1 4.5
##  [163] 4.7 4.1 3.7 4.6 3.2 4.0 4.2 4.7 4.4 3.7 3.0 4.4 6.1 6.2 5.6 2.9 4.3 4.7
##  [181] 4.9 4.9 4.0 4.2 4.6 2.7 4.2 4.1 4.2 3.4 5.1 5.8 3.1 4.0 3.9 4.9 4.6 4.6
##  [199] 4.6 5.4 4.3 3.5 4.0 5.6 3.5 4.5 4.8 4.9 3.8 3.7 5.6 3.7 3.4 5.4 3.8 4.0
##  [217] 3.2 3.5 4.6 7.1 2.6 3.1 3.4 3.2 5.1 5.3 4.5 5.3 4.6 2.9 5.5 3.5 2.9 4.3
##  [235] 5.3 4.5 4.6 4.8 4.3 4.6 5.1 4.0 5.0 3.9 4.8 4.2 4.1 4.0 5.2 2.5 4.9 3.9
##  [253] 3.4 4.5 5.2 4.5 4.5 4.3 5.9 3.8 3.9 6.3 4.3 4.4 3.5 3.5 3.0 5.3 4.0 5.8
##  [271] 5.2 4.7 3.9 5.4 4.2 2.9 4.9 4.1 5.8 4.4 4.5 3.3 4.4 4.9 4.8 4.6 5.6 4.9
##  [289] 5.8 5.5 4.2 3.6 5.0 5.7 5.4 4.2 4.2 4.7 4.2 4.7 2.8 3.7 4.2 5.5 4.5 4.4
##  [307] 5.6 4.3 4.5 4.8 4.7 4.7 5.4 2.9 4.6 3.5 4.9 4.4 5.5 5.3 3.3 3.5 5.1 4.5
##  [325] 3.8 4.5 5.0 4.7 4.7 5.5 5.0 4.6 5.5 3.3 4.7 3.4 4.3 5.2 4.5 5.8 4.3 5.8
##  [343] 3.7 4.7 4.8 5.7 4.4 3.6 4.2 4.6 2.7 4.5 4.8 4.1 4.7 4.2 3.9 5.6 3.4 5.1
##  [361] 3.0 4.7 5.2 3.2 4.6 5.3 3.8 5.3 4.2 3.0 4.9 3.7 3.5 4.0 3.3 2.5 4.6 3.3
##  [379] 4.0 3.3 3.5 3.8 5.8 3.8 4.6 4.0 6.0 5.2 4.9 4.6 4.7 5.7 5.0 4.9 5.3 4.8
##  [397] 4.7 3.7 5.7 5.0 4.8 3.8 5.7 3.3 3.7 3.3 4.3 2.8 3.7 4.3 4.4 3.8 5.2 5.1
##  [415] 3.2 4.1 3.7 5.1 3.1 4.4 3.8 5.0 4.7 5.4 5.5 5.4 5.3 5.1 4.6 4.6 4.0 4.3
##  [433] 3.5 3.2 3.6 5.3 4.5 4.1 4.7 4.4 6.7 5.7 5.1 4.4 4.2 4.7 5.5 6.9 3.1 5.4
##  [451] 5.4 5.1 5.2 4.3 5.9 3.6 5.2 2.7 4.4 4.6 6.4 6.5 5.5 5.3 4.8 4.0 4.4 4.1
##  [469] 5.3 4.4 3.3 6.2 4.8 5.5 4.1 3.0 2.1 4.8 4.0 5.3 4.5 5.4 5.3 4.6 3.9 4.6
##  [487] 3.4 3.9 4.4 4.3 4.7 5.1 6.1 4.5 4.3 3.3 3.8 4.1 3.7 5.0 5.5 4.0 5.0 4.5
##  [505] 4.7 4.1 3.8 3.1 4.2 3.3 5.2 4.7 5.4 5.0 4.4 3.9 4.8 3.8 4.3 5.1 3.2 3.4
##  [523] 5.1 2.9 4.3 4.9 6.1 3.5 3.8 4.9 4.1 4.9 4.7 3.8 5.4 3.8 3.2 3.5 4.5 3.0
##  [541] 4.0 3.1 4.2 4.6 4.2 4.1 3.9 4.4 4.3 4.6 3.0 5.1 4.2 3.5 4.1 5.1 5.3 2.6
##  [559] 4.4 5.8 4.5 4.6 4.1 3.3 4.4 5.3 4.7 3.7 5.3 5.5 5.7 5.3 4.5 5.9 4.2 3.2
##  [577] 6.2 4.1 3.5 3.9 3.2 5.1 5.6 5.6 4.9 5.5 4.8 6.1 4.0 5.4 4.4 6.1 4.4 4.1
##  [595] 6.5 4.1 4.1 4.8 4.0 4.4 5.8 2.4 3.2 5.8 5.2 4.8 3.8 4.8 6.2 5.2 2.2 5.0
##  [613] 4.9 5.1 4.3 5.3 4.2 5.0 5.3 5.3 3.2 5.4 4.9 5.0 4.8 5.1 4.6 5.3 5.6 4.6
##  [631] 4.1 3.7 4.3 3.7 3.2 4.7 4.3 4.4 3.3 3.9 4.5 4.2 4.2 2.4 3.8 4.5 4.0 3.8
##  [649] 4.8 4.7 5.0 6.0 5.3 3.8 3.6 4.4 4.4 4.6 4.7 4.3 4.3 5.1 4.2 3.7 5.3 5.5
##  [667] 4.0 5.9 4.4 4.1 3.2 5.2 3.0 4.7 5.6 4.5 4.8 4.6 5.1 5.1 5.1 4.6 5.0 2.0
##  [685] 4.9 2.5 6.0 6.3 3.2 4.0 4.3 4.4 4.1 5.4 4.8 5.9 4.0 4.8 4.7 5.4 5.7 5.0
##  [703] 4.5 5.2 5.1 4.5 5.9 5.2 4.5 3.7 5.3 4.5 3.6 2.6 5.0 4.0 4.7 4.0 3.6 4.9
##  [721] 5.5 4.8 4.9 5.9 4.5 4.5 4.1 4.5 4.0 4.3 6.7 3.8 5.9 5.1 5.1 4.9 6.4 5.4
##  [739] 4.3 2.5 3.2 5.4 4.1 7.3 4.5 4.3 5.5 4.8 3.7 3.9 5.3 3.7 4.7 5.0 6.0 4.1
##  [757] 5.7 5.7 3.4 4.0 4.6 5.7 2.1 3.9 5.1 4.0 3.9 2.3 2.6 4.4 3.8 4.7 3.9 3.1
##  [775] 4.4 4.2 3.8 4.8 5.0 4.3 4.8 5.7 4.5 4.4 4.8 4.7 3.9 4.7 4.0 3.6 4.1 5.5
##  [793] 4.8 4.7 3.1 2.6 4.0 3.7 4.1 4.3 5.4 4.5 4.8 4.1 4.3 4.1 3.8 5.2 4.5 3.6
##  [811] 4.2 4.0 3.3 4.1 5.9 5.2 4.0 4.2 3.9 4.0 5.5 3.6 4.8 4.6 4.9 3.4 4.8 5.5
##  [829] 5.1 4.7 4.8 6.6 5.6 4.0 4.8 4.4 4.7 5.1 5.7 2.6 5.9 4.7 4.8 3.5 5.0 3.6
##  [847] 4.0 5.5 4.5 5.5 4.7 4.3 2.7 5.5 2.1 3.5 4.5 4.0 5.0 4.4 4.7 6.2 5.5 4.4
##  [865] 2.8 3.6 6.0 6.3 3.7 5.3 4.3 4.1 4.3 3.9 3.1 3.9 4.4 5.0 5.4 4.4 3.4 4.9
##  [883] 6.5 5.3 4.5 5.9 3.8 3.4 3.6 7.1 3.5 4.5 2.8 4.1 6.6 4.4 4.2 4.0 4.7 4.6
##  [901] 3.5 4.4 3.1 4.3 4.9 4.6 4.7 5.2 3.8 4.1 6.4 5.1 3.6 4.7 5.3 3.9 4.9 4.6
##  [919] 4.7 4.4 3.7 3.5 3.1 3.0 5.3 3.4 4.4 4.9 4.9 4.7 4.0 5.4 4.3 6.0 4.3 5.1
##  [937] 5.1 4.0 5.8 2.9 3.7 4.8 3.8 5.6 4.4 5.8 3.9 4.9 3.5 5.6 4.6 4.2 3.8 5.3
##  [955] 5.4 4.0 4.9 5.1 6.0 5.5 2.6 4.5 4.6 3.5 5.0 4.5 4.9 3.5 5.2 5.0 3.8 4.3
##  [973] 4.8 4.8 3.7 5.2 6.1 4.4 3.6 6.0 4.4 5.6 5.3 4.1 5.2 5.0 4.2 3.8 4.9 4.9
##  [991] 4.0 5.2 5.6 5.0 5.9 4.5 4.8 3.8 4.8 4.2
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.3 5.3 5.1 5.0 5.0 4.9 4.7 4.6 4.5
## 
## $jack.boot.se
## [1] 0.8777243
## 
## $call
## bootstrap(x = x, nboot = 1000, theta = theta, func = perc72)
```

On Tuesday we went over an example in which we bootstrapped the correlation coefficient between LSAT scores and GPA. To do that, we sampled pairs of (LSAT,GPA) data with replacement. Here is a little script that would do something like that using (X,Y) data that are independently drawn from the normal distribution


```r
xdata<-matrix(rnorm(30),ncol=2)
```

Everyone's data is going to be different. With such a small sample size, it would be easy to get a positive or negative correlation by random change, but on average across everyone's datasets, there should be zero correlation because the two columns are drawn independently.


```r
n<-15
theta<-function(x,xdata)
  {
  cor(xdata[x,1],xdata[x,2])
  }
results<-bootstrap(x=1:n,nboot=50,theta=theta,xdata=xdata) 
#NB: xdata is passed to the theta function, not needed for bootstrap function itself
```

Notice the parameters that get passed to the 'bootstrap' function are: (1) the indexes which will be sampled with replacement. This is different that the raw data but the end result is the same because both the indices and the raw data get passed to the function 'theta' (2) the number of bootrapped samples (in this case 50) (3) the function to calculate the statistic (4) the raw data.

Lets look at a histogram of the bootstrapped statistics $\theta^{*}$ and draw a vertical line for the statistic as applied to the original data.


```r
hist(results$thetastar,breaks=30,col="pink")
abline(v=cor(xdata[,1],xdata[,2]),lwd=2)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-17-1.png" width="672" />

Parametric bootstrap
---------------------

Let's do one quick example of a parametric bootstrap. We haven't introduced distributions yet (except for the Gaussian, or Normal, distribution, which is the most familiar), so lets spend a few minutes exploring the Gamma distribution, just so we have it to work with for testing out parametric bootstrap. All we need to know is that the Gamma distribution is a continuous, non-negative distribution that takes two parameters, which we call "shape" and "rate". Lets plot a few examples just to see what a Gamma distribution looks like. (Note that the Gamma distribution can be parameterized by "shape" and "rate" OR by "shape" and "scale", where "scale" is just 1/"rate". R will allow you to use either (shape,rate) or (shape,scale) as long as you specify which you are providing.

<img src="Week-2-lab_files/figure-html/unnamed-chunk-18-1.png" width="672" />


Let's generate some fairly sparse data from a Gamma distribution


```r
original.data<-rgamma(10,3,5)
```

and calculate the skew of the data using the R function 'skewness' from the 'moments' package. 


```r
library(moments)
theta<-skewness(original.data)
head(theta)
```

```
## [1] -0.1135385
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##     shape       rate  
##   5.527087   9.584470 
##  (2.400900) (4.358411)
```

Now lets sample with replacement from this new distribution and calculate the skewness at each step:


```r
results<-c()
for (i in 1:1000)
  {
  x.star<-rgamma(length(original.data),shape=fit$estimate[1],rate=fit$estimate[2])
  results<-c(results,skewness(x.star))
  }
head(results)
```

```
## [1] 0.2552544 1.5048163 0.6655179 0.7069533 0.5778461 0.6820755
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-22-1.png" width="672" />

Now we have the bootstrap distribution for skewness (the $\theta^{*}$ s), we can compare that to the equivalent non-parametric bootstrap:


```r
results2<-bootstrap(x=original.data,nboot=1000,theta=skewness)
results2
```

```
## $thetastar
##    [1] -0.4565122266  0.2797629630 -0.7463174997 -0.1407450953  0.0140556591
##    [6] -0.4542025180  0.2778877491  0.2331933854  1.4323266026 -0.1826555213
##   [11] -1.3676504469  1.0045630093  0.0477319324 -0.7153806657 -0.2625512677
##   [16] -0.1558296742 -1.6037619581 -0.3776726792 -0.0829454849  0.2888300721
##   [21]  0.4797761972 -0.4734391660 -0.1972876637  0.4152507332 -0.2618757423
##   [26] -0.3455773724  0.7925306267 -0.2852498474 -0.0419285993 -0.2587608958
##   [31]  0.4298559436 -0.6827798529 -0.0007692209 -0.4305400287  0.1344855019
##   [36] -1.4061287783  0.3208884284 -0.1977560360 -0.1066244493  0.4340642497
##   [41]  0.1660300868  0.3045996705 -1.5206822475 -0.4224280311 -0.7909889121
##   [46] -0.0620125320  0.1881939255 -1.2552311789 -0.0059979808 -0.9716656275
##   [51] -0.4085572047 -0.2193106522 -1.2111673064  0.2942098420 -0.7435022038
##   [56] -0.4571647608 -0.5412922472 -0.2796156444 -0.2500045288 -0.0201711525
##   [61] -0.0747969368 -1.3703651318 -0.6064536368  0.6993089229 -0.1346014918
##   [66]  0.3290451453  0.3127266962  0.1368072984 -0.3106020677  0.3037146098
##   [71] -0.2119555750 -0.0645052566 -0.2032442068 -0.0228120112  0.0980414880
##   [76] -0.4079899101 -0.1821833562 -0.1482157369 -0.0689671068  0.2426185000
##   [81]  0.1213998552 -1.0959864431  0.0015118573  0.0842913623  0.3687036341
##   [86] -0.0521939181 -0.4712472763 -0.5735422781 -0.4691872017 -0.4490190641
##   [91]  0.3193959989  0.1079829621 -0.1003452295  0.0107170428 -0.2090418735
##   [96] -0.0063948339  0.1175574003  0.3987642944  0.3889389203 -0.2568183905
##  [101] -0.0403476626  0.9261992226 -0.4413863704 -0.0968484100 -0.2723473623
##  [106] -1.1815286757  0.1544179876 -0.6991752429  0.0450351432 -0.2785802366
##  [111]  0.1868456641 -0.3201985050 -0.0753339642 -0.1116975537  0.1827983634
##  [116] -0.0716130644 -0.3586115715  0.4993388411 -0.4345638042 -0.6817775414
##  [121]  0.8487074700 -0.0302423458  0.1288884903 -0.3180436965 -0.2098701928
##  [126] -0.1582476503  0.0362649907 -0.4581151608  0.6862967728 -0.8188360447
##  [131] -0.1129619956 -0.5225791646  0.0235730587 -0.4109377222 -0.1918100630
##  [136] -0.3104288826  0.5053564450 -0.1966282969  0.0367533745 -0.5271362665
##  [141]  0.4727887357  0.3419177210 -0.1585012145  0.1294258233  0.0827187707
##  [146] -0.4365964475  0.3755769086  0.1843904327 -0.1785546821 -0.6372660885
##  [151]  0.2242129856 -0.9058015781  0.2890308790 -0.4131519638 -0.6840710848
##  [156] -0.2645935869  0.4540442226  0.1245853411 -0.1044290508 -0.4785818037
##  [161] -0.6605642152 -0.1047271606 -0.4370416711 -0.1291559498 -0.7516739290
##  [166]  0.4562248161 -0.1376309498 -0.9509301733  0.2936771234 -0.6550070141
##  [171]  0.4275219280 -0.6843214726 -0.7747548335 -0.2124858359 -0.1931149519
##  [176] -0.1632452880 -0.0304428240 -0.4010709262  0.1086227072  0.0536891777
##  [181] -0.2263684356 -0.0186170174  0.2811527926 -0.3414916493 -0.3069095328
##  [186]  0.4921381559 -0.0445074498  0.5394030073  0.4261256320 -1.1407696079
##  [191] -0.7804651707  0.3995380288  0.2610280748 -1.0266945780 -0.0148568641
##  [196] -1.0209348537  0.0324322104  0.2403859220 -0.2082721916 -0.7805800331
##  [201] -0.3971634681 -0.5579129727 -0.6228342854 -0.0197184585  0.0713269839
##  [206] -0.2778678966  0.1589824394 -0.9679748612 -0.1116084557 -0.2366759520
##  [211] -1.0483764571  0.2796273729  0.1805925488  0.1403120606 -0.5776934543
##  [216] -0.2446117511 -0.2069444638  0.5556815159 -0.3036726149 -0.8979517381
##  [221]  0.1277422561 -0.4531904595  0.4820320883  0.2268350575 -0.9903858312
##  [226]  0.1176636588 -0.6630552174  0.3556218160  0.5621349137 -0.1274126821
##  [231]  0.1083853259 -0.7136360112  0.3546897394  0.1122920341 -0.2390012507
##  [236] -0.4557735761  0.3766241579 -0.0328780741 -0.2818309271 -0.2803568638
##  [241]  0.5222434974  0.5749547896  0.1590798730 -0.3360370711 -0.9903858312
##  [246]  0.5334844372 -1.0578934466  0.7187555640 -0.3874480439 -0.8291957296
##  [251] -1.0208667315  0.6932183676 -0.4863866363 -1.1130465216  0.0300879788
##  [256]  0.2402025865  0.2937674847 -0.2232565391 -0.0745444684 -0.7825187077
##  [261] -0.1969817093  0.0535968737  0.7719413396  0.1721352852 -0.5699092256
##  [266] -0.3568372106 -0.7876375486 -0.1906763183  0.4001384497 -0.2796989107
##  [271] -0.3375856230 -0.0787044676 -0.5870161431  0.4739359456  0.4857859162
##  [276]  0.1510238096 -0.7258900751 -0.1571037435 -0.0219660881 -0.6090748793
##  [281]  0.2209630039 -0.9629115814  0.2367802300  0.2975849294 -0.4710034399
##  [286]  0.3097956369 -0.2696421440 -0.0425123603 -0.6424844736  0.1134836629
##  [291] -0.1258351733 -0.2886650750 -0.1064232523  0.2316323796  0.5932309101
##  [296] -0.2841300233 -0.5595486644  0.1448937418  0.0047442826 -0.3108839542
##  [301] -0.1173904826  0.5841243154  0.2353901312 -0.1915153515 -0.3167868501
##  [306]  0.3265158359  0.1361936907  0.5663606226 -0.0093763418  0.8029942667
##  [311] -0.6354392972 -0.1711712442 -0.0952772943 -0.4308225382  0.4424007122
##  [316]  0.3208884284 -0.1926187447  0.5892612682  0.3741145601  0.1933080347
##  [321]  0.0758755437 -0.5189139293  0.3516494317 -0.3963716767 -0.2707785543
##  [326]  0.8226503760  0.4581887811  0.7837016639  0.0565366891 -0.7552441823
##  [331]  0.0001405620  0.8962176436 -0.2228493928 -0.7427123951 -0.1112704372
##  [336]  0.1269930030 -1.3229533547  0.2584966335 -0.0309793098  0.2768705543
##  [341]  0.3620626800  0.1368072984  0.9163073770  0.8851667351  0.7767202232
##  [346]  1.5743698119  0.1916820346  0.2016388510 -0.1478001014  0.3029682011
##  [351] -0.4296323560  0.0775220431  0.1911553747 -0.4033710652 -0.2001507592
##  [356] -0.4862235211 -0.0832958271 -1.2672072395 -0.2108574825 -1.3766983935
##  [361] -0.9691101599 -1.1290650860 -0.1280873224 -0.4909517345 -1.4833949468
##  [366]  0.2002826738 -0.4469063810 -1.0569295376  0.8746053778  0.2926866629
##  [371]  0.2191785155  0.5039294431 -0.0446218999  0.4492111837  0.1667582136
##  [376] -0.2281876929 -0.0331136795  0.0628630637  0.9624799410  0.3025939621
##  [381]  0.5617018205 -0.0338108143 -0.4671599236 -0.6427504823 -0.1733090472
##  [386]  0.2012011849 -0.9092050201 -1.1241006917 -0.0750771592  0.5928905234
##  [391]  0.3244352403 -0.5858378635 -0.5921180745  0.1801832926  0.0898529498
##  [396]  0.2555140491 -0.4411289821  0.5049509364 -0.4507124765  0.2440484588
##  [401]  0.4610433951  0.2890787117 -0.5548225957  0.5425192547 -0.1870446454
##  [406] -0.1632448540 -0.0973630596 -0.3271469883 -0.2245246143 -0.5753657675
##  [411]  0.1499878123 -0.4451708391 -0.2074496071  0.1305870988 -1.1730365343
##  [416]  0.0288027794  0.4850874723 -0.6733579517 -0.1494550115  0.3248477862
##  [421]  0.6072217236  0.4199489118 -0.4022288127  0.3110617995 -0.4498305398
##  [426] -0.2997451882 -0.9604997858 -0.0034041263  0.2094583802  0.4912817055
##  [431] -0.7694168166 -0.3418739954 -0.3833825665  0.3928110453 -0.4790386285
##  [436] -0.2785802366 -0.0225834149 -0.2509532426 -0.1050912347 -0.7434208390
##  [441]  0.8388064983 -1.0432334343 -0.0483437779 -0.8473706903 -0.5401910074
##  [446]  0.5359322231 -0.8440451009 -0.1559416807  0.5419140140 -1.4470394839
##  [451]  0.1715344197  0.0928996117 -0.5422211669 -0.0117650927 -0.0863194054
##  [456] -0.8398085716  0.0908733159 -0.1789305317 -0.2830725298 -0.0235492163
##  [461]  0.0061309144 -0.0526322159  0.2556194359  0.0551978942  0.2203804186
##  [466]  0.7314460368  0.2565610299  0.8316402193 -0.5463603415 -0.4776447179
##  [471]  0.1237075221 -0.1446471857 -0.3897643212 -0.7265331527 -0.8858241876
##  [476] -0.5874571320 -0.4435138260  0.8656428062 -0.5453558881 -0.5942466657
##  [481]  1.5519350271 -0.1092410541 -0.4109491502  0.5895354766  0.4856426399
##  [486]  0.9926604180 -0.5378296912  0.5274076784  0.0514383113 -0.0455752292
##  [491]  0.3472144418  0.2235613434 -0.4391465706 -0.1845859250  0.5192484834
##  [496] -0.7082656031 -0.6142907447  0.2127793769  0.3021192198  0.4423921905
##  [501] -0.2736098346 -0.3177579597 -0.3670893012 -0.3531362153  0.2112239903
##  [506]  0.3393501689  0.6352315850 -0.0091291495  0.1971128289 -0.0446708968
##  [511] -1.2134586418  0.4317356406  1.1056556899 -0.4435138260  0.0981197331
##  [516]  0.2948111183  0.7566960697 -1.2820628505  0.1086915272  0.4723996199
##  [521] -1.1090561719 -0.1787502753  0.0159944956 -0.2547592175  0.0012672770
##  [526]  0.3430671474  0.2666149844 -0.6828059729 -0.0840769288  0.1894269110
##  [531] -0.0906797133  0.8384203729  0.0995745788  0.2985573496 -0.0585340445
##  [536]  1.2550212251 -0.0578629686  0.2410271280 -0.1274958113 -0.2158891155
##  [541] -0.7637140619 -0.0252307802  0.0547317090  0.0891908723  0.1817949124
##  [546] -0.1291559498 -0.1761759437 -0.8784488834  0.6765924243 -0.5156619589
##  [551] -0.9104263735 -0.3394562534 -0.4308225382  0.3714872991 -0.3103949062
##  [556]  0.1667582136  0.5851402066  0.9431048962  1.3554323560  0.1294178669
##  [561]  1.3861308546 -0.1723645336  0.4570911460  0.0640888007 -0.4891330717
##  [566]  0.6027456734 -0.2512664455 -1.0909842817  0.1705702515 -0.6392663031
##  [571] -0.0754751681  0.2553270279 -0.2319202004 -0.5143604903  0.6035302144
##  [576]  1.1808430771 -0.3849655648 -0.5004000772 -0.2831089651 -0.0022315492
##  [581] -0.5543656632 -0.2048603748 -0.1657363533  0.1945832849  0.5617958089
##  [586] -0.8901232941  0.0019961233  0.1118380242  0.4599178184 -0.7067250832
##  [591] -0.1462720627  0.3028964103 -0.6095675481 -0.0347402843  0.0476078273
##  [596]  0.4269116111  0.5992481121 -0.0322606111  0.5746475440  0.1052192742
##  [601] -1.8113270469 -0.8355225013 -0.2070539519 -0.6779784805  0.6807429228
##  [606]  0.5903562097 -0.2313030196 -0.4608118382  0.0666957106 -0.0235497653
##  [611] -0.5211221709  0.1035284097  0.2269751775  0.5819622191 -0.4198011690
##  [616] -0.0019124247  0.3093339776  0.2567269824  1.0331144688  0.6039869147
##  [621] -0.5909807425  0.3022126102 -0.9584648090 -0.7513676772 -0.0097898540
##  [626]  0.3175013672  0.2407098613 -0.3532166114  0.0668517157  0.0196717976
##  [631]  0.0447407548 -0.1051683870 -0.5189729875 -0.0303539952 -0.2304632371
##  [636]  0.2822064464  0.1654847755 -0.7552822983 -0.2072225796 -0.5335343516
##  [641]  0.4966637332 -0.3887186699 -0.4399732944 -0.3285556011  0.0128041331
##  [646] -0.4975608968  0.3622496444 -1.0660035380 -1.1105999360 -0.3849655648
##  [651] -0.2202479667 -0.0352414426  0.2338943678 -0.5109739119 -0.5759183908
##  [656] -0.6065547674  0.4779184965  0.1286534578 -0.2484058406  0.1527495760
##  [661] -0.5609382408 -0.5464585190 -0.0452921293 -0.3136968908  0.1592562595
##  [666]  0.3995851405 -0.1573026695 -0.2208248394  0.5889516952 -0.5725744680
##  [671] -0.2398076044 -0.4500354604 -0.2648617779  0.3182689956  0.5605620178
##  [676] -1.0127546064  0.7378933728 -0.0997109453 -0.2955208414  0.1893990546
##  [681] -0.5058035977  0.0006271068  0.0808583996  0.2273424284 -0.5160618997
##  [686] -0.0578629686 -0.1374225771  0.5345870083  0.3921855173 -0.4321431599
##  [691] -0.5695496158 -0.0515632324 -0.1135384779  0.1215227526  0.0425840530
##  [696]  0.1475466347  0.1963914631  0.0028852037 -0.1971159462 -0.3513110892
##  [701] -0.9895774003  0.0122376355 -0.3446773088 -0.8300298613  0.9191469453
##  [706] -1.1722668185 -0.0493102892  0.5546465119  0.6223245732 -0.5508047788
##  [711] -0.4281715069 -0.2516147128  0.2255475719  0.4130332939  0.4076214265
##  [716]  1.0835059286 -0.8094222490  0.6941071405 -0.4033710652 -0.4865744190
##  [721] -0.7993402637  0.0098809764  0.0955667413 -0.2448851437  0.4510847332
##  [726] -0.0948187988 -0.4478862213  0.6614720090 -0.8200202655  1.0049984628
##  [731]  0.1767861104  0.0147383856 -0.4245051569  0.2183801864  0.1753754252
##  [736]  0.0371139738 -0.0492377354 -0.0725432959  0.2219088266  0.2380618330
##  [741]  0.5748203893 -0.2792744298 -0.5575620960  0.3996484012  0.3334204901
##  [746]  0.8771681324 -0.1657497917 -0.3192881141  0.9231008276 -0.2324925965
##  [751] -0.3691148387  0.4176417416 -0.4965638934 -0.0770156940 -0.0175760689
##  [756]  0.0938402002 -0.6378997265  1.2039853607  0.6719109068  0.5319100299
##  [761]  0.4524953014 -1.0941494671  0.5412120395  0.2006803424  0.4508741474
##  [766] -0.1728542086  0.3213100806 -0.4965336899 -0.9767821253 -0.9037399527
##  [771]  0.1112243194  0.4333743135 -0.8993428826 -0.4956301092  1.0597073347
##  [776] -0.5113755608  0.3987642944 -0.7158920897  0.1041076211 -0.4109491502
##  [781]  0.0524975118  0.2612448776  0.5281813575 -0.2762962649 -0.8230172878
##  [786]  0.6261301683  0.0344561916  0.4738401176 -0.0638050377 -0.9503310308
##  [791] -0.6648853292 -0.3882500652 -1.1685444864  0.8212260386 -0.1118809572
##  [796]  0.1115398260 -0.1876188128 -1.0273692730  0.3051201215  0.3657681620
##  [801] -0.1839937471  0.4682614765 -0.4456383548  0.1836039902  0.4579430129
##  [806] -0.5303548356  0.3051201215 -0.5238074735  0.0115048897  0.2928826179
##  [811] -0.4351842243  0.1764769536 -0.6540656042  0.1244484737  0.4757670142
##  [816] -0.1527916641  0.3434786413 -0.3763610649 -0.7232028990  0.2796273729
##  [821] -0.5710606666 -0.2196460348  0.3006906903  0.2510798830 -0.4171982067
##  [826]  0.0625563049  0.4697575980 -0.7230287090  0.0060982124 -0.1932259628
##  [831] -0.0245903504 -0.5016670985  0.5038961915 -0.2016069879 -0.6954022857
##  [836]  0.2111659435  0.5280779636 -0.5876127178 -1.0215542674 -0.3555733974
##  [841]  0.1112799995 -0.2719959548  0.0428945219 -0.2056546249  0.3194933661
##  [846]  0.3679924809  0.8566884578 -0.7773172019 -0.8207646083 -0.3182948647
##  [851] -0.3594976346 -0.4490603581 -0.5999312812 -0.1339015443 -0.2542594862
##  [856] -1.1315292024  0.1286464289  0.5689428239  0.7089204086  0.8987464695
##  [861] -0.4622835832  0.3256482751 -0.3738360219  0.6855963945  0.2854441821
##  [866]  0.1590458131 -0.3088583326 -0.5626183470  0.0422380575 -0.9507411066
##  [871]  0.5410232645 -0.2875135620 -0.4366219061 -0.1829034227  0.4062755938
##  [876] -0.5275982299 -0.4945999726  0.0619223753 -0.0136422420  0.1876154804
##  [881] -0.4412121602 -0.0491888292 -0.1680473213 -0.7668672318  0.2992064657
##  [886] -0.1966150879 -1.3265953840  0.0845552570  0.9869043786  0.5284011717
##  [891] -0.2050577984 -0.5880578621 -0.2779565934 -0.1073719215  0.0720433395
##  [896] -0.0710289479 -0.0653160100 -0.8948854622  0.1576830588 -0.8157593195
##  [901]  0.5707869631  0.5419361797 -0.2215609277  0.2483846309  0.2885385801
##  [906] -0.4037120143 -0.1790642146  0.0763438447 -0.6056746085 -0.4934786239
##  [911]  0.3164048571 -0.2153358491 -0.5944253863 -1.0595727395 -0.2069444638
##  [916]  0.3657616247 -0.2454330148 -0.1790593705  0.3914993289  0.3705536831
##  [921] -0.0172948823 -0.7142658773  0.0511011866 -1.4353655418  0.5935069321
##  [926] -0.1199596236  0.4856345361 -0.3460454468  0.3180346183 -0.1668267094
##  [931]  1.3862163406  0.5014610964  0.1466206431  0.2240016260  0.5429518377
##  [936]  0.3921855173  0.5243532322  0.4988754713 -0.1357916819 -0.3312148077
##  [941] -0.0069496676  0.3265158359  0.0856171817  0.4314728103 -1.9236520042
##  [946] -0.3366037357  0.2441533681  0.4148282403  0.2078159835 -0.2862598466
##  [951] -0.0115724380  0.0514642316 -0.9764294149  0.2135686672  0.1643537176
##  [956] -0.1420430683 -0.1859653458  0.1052192742  0.2158823633 -0.5331153391
##  [961] -0.4296323560  0.2555179485  0.7430488500 -0.5575786947 -0.4317039966
##  [966]  0.5329229938  0.0973966765 -0.3007983359 -0.2472213202  0.3582739962
##  [971] -0.3637261404 -0.2438809803 -0.8666000042 -0.0068715697 -0.2403881310
##  [976] -0.6170559822  0.3086390867  0.3880044971 -0.5636690958 -0.3177579597
##  [981] -0.2764488981 -0.8216071160  0.5319611196  0.6862967728 -0.5318589481
##  [986]  0.4678998811 -0.8349539664 -0.6934023333  0.4825507595 -0.7427123951
##  [991] -0.0511072163 -0.2162659066  0.0689770127 -0.0607373923 -0.2845934479
##  [996] -0.2911754666  0.3308058282 -0.8351753501 -0.1772012802  0.7652554932
## 
## $func.thetastar
## NULL
## 
## $jack.boot.val
## NULL
## 
## $jack.boot.se
## NULL
## 
## $call
## bootstrap(x = original.data, nboot = 1000, theta = skewness)
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-23-1.png" width="672" />

What would have happened if we would have fit a normal distribution instead of a gamma distribution?


```r
fit2<-fitdistr(original.data,dnorm,start=list(mean=1,sd=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced

## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
fit2
```

```
##       mean          sd    
##   0.57666709   0.22382344 
##  (0.07077919) (0.05004521)
```

```r
results.norm<-c()
for (i in 1:1000)
  {
  x.star<-rnorm(length(original.data),mean=fit2$estimate[1],sd=fit2$estimate[2])
  results.norm<-c(results.norm,skewness(x.star))
  }
head(results.norm)
```

```
## [1]  0.7865789 -0.1714547 -0.4559057 -0.3748067 -0.4779693  0.2976969
```

```r
hist(results,breaks=30,col="pink",ylim=c(0,1),freq=F)
hist(results.norm,breaks=30,col="lightgreen",freq=F,add=T)
hist(results2$thetastar,breaks=30,border="purple",add=T,density=20,col="purple",freq=F)
```

<img src="Week-2-lab_files/figure-html/unnamed-chunk-24-1.png" width="672" />

All three methods (two parametric and one non-parametric) really do give different distributions for the bootstrapped statistic, so the choice of which method is best depends a lot on the situation, how much data you have, and what you might already know about the underlying distribution.

Jackknifing is just as easy at bootstrapping. Here we will do a trivial example for illustration. We will write a little function for the mean even though you could put the function in directly with 'jackknife(x,mean)'


```r
theta<-function(x)
  {
  mean(x)
  }
x<-seq(0,9,by=1)
results<-jackknife(x=x,theta=theta)
results
```

```
## $jack.se
## [1] 0.9574271
## 
## $jack.bias
## [1] 0
## 
## $jack.values
##  [1] 5.000000 4.888889 4.777778 4.666667 4.555556 4.444444 4.333333 4.222222
##  [9] 4.111111 4.000000
## 
## $call
## jackknife(x = x, theta = theta)
```

Why do we not have to tell the 'jackknife' function how many replicates to do?

Let's compare this with what we would have obtained from bootstrapping


```r
results2<-bootstrap(x,1000,theta)
mean(results2$thetastar)-mean(x)  #this is the bias
```

```
## [1] -0.0369
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9151863
```


Everything we have done to this point used the R package 'bootstrap' - now lets compare that with the R package 'boot'. To avoid any confusion (a.k.a. masking) between the two packages, I recommend detaching the bootstrap package from the workspace with


```r
detach("package:bootstrap")
```


The 'boot' package is now recommended over the 'bootstrap' package, but they give the same answers and to some extent it is personal preference which one prefers to use.

We will still use the mean as the statistic of interest, but we will have to write a new function for it because the syntax of the 'boot' package is slightly different:


```r
library(boot)
theta<-function(x,index)
  {
  mean(x[index])
  }
boot(x,theta,R=999)
```

```
## 
## ORDINARY NONPARAMETRIC BOOTSTRAP
## 
## 
## Call:
## boot(data = x, statistic = theta, R = 999)
## 
## 
## Bootstrap Statistics :
##     original     bias    std. error
## t1*      4.5 0.04874875   0.8899213
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 4 5 6 7 9 
## 2 2 2 1 1 2
```

```r
xmeans<-vector(length=1000)
for (i in 1:1000)
  {
  xmeans[i]<-mean(sample(x,replace=T))
  }
mean(x)
```

```
## [1] 4.5
```

```r
bias<-mean(xmeans)-mean(x)
se.boot<-sd(xmeans)
bias
```

```
## [1] 2e-04
```

```r
se.boot
```

```
## [1] 0.8984694
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

