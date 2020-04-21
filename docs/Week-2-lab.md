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
## 0 1 2 3 7 8 9 
## 1 1 3 1 1 2 1
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
## [1] 0.0229
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
## [1] 2.74793
```

```r
UL.boot
```

```
## [1] 6.29787
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.9   6.3
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
##    [1] 4.0 5.3 5.5 6.7 2.7 3.7 4.4 3.5 4.9 2.8 5.0 3.7 4.8 3.2 3.0 3.5 4.9 5.3
##   [19] 4.0 4.8 5.1 4.6 3.3 5.3 3.1 4.2 4.4 2.3 5.1 4.5 4.4 3.7 5.5 5.3 3.4 4.6
##   [37] 5.1 4.0 3.4 4.6 4.5 5.4 4.4 6.4 4.0 4.0 4.0 6.0 5.3 4.7 3.9 5.6 4.4 4.8
##   [55] 4.9 4.9 3.9 4.2 3.4 2.7 4.7 3.7 4.6 3.9 3.6 3.6 4.7 6.7 6.0 2.9 5.7 5.1
##   [73] 5.3 5.5 3.6 4.3 4.3 5.2 5.3 4.1 5.2 5.6 3.2 6.0 2.6 5.4 4.5 4.2 4.5 4.1
##   [91] 3.3 4.6 3.0 4.8 3.2 5.4 4.4 4.0 4.9 4.9 3.7 4.7 5.2 4.7 3.4 4.1 4.9 4.1
##  [109] 4.7 3.9 4.4 5.5 6.6 4.0 3.6 3.1 4.3 5.6 4.6 1.8 4.3 5.3 3.8 4.1 4.6 4.3
##  [127] 4.9 4.0 4.8 3.6 4.2 4.1 4.5 6.2 5.0 2.9 4.4 4.6 4.8 6.0 4.3 3.4 4.1 6.3
##  [145] 3.8 3.6 1.9 5.3 3.7 6.6 4.2 5.6 4.7 4.8 3.5 4.6 4.6 3.8 4.1 6.0 4.7 5.8
##  [163] 3.5 5.0 5.6 3.9 5.9 3.7 3.7 3.4 4.3 4.0 4.8 6.1 3.6 4.8 4.9 3.3 4.3 5.3
##  [181] 3.2 3.5 3.6 5.1 4.3 3.8 4.4 5.0 4.8 4.8 4.8 4.6 3.8 4.3 4.5 4.7 4.9 4.2
##  [199] 5.4 3.4 5.6 4.4 4.9 3.6 5.6 5.0 5.3 4.7 4.9 4.9 2.7 5.4 5.0 5.2 5.1 2.8
##  [217] 4.0 4.3 4.5 2.5 5.6 4.0 4.1 6.7 5.0 6.6 5.5 5.1 3.9 4.8 5.1 5.4 4.9 5.1
##  [235] 5.7 4.9 5.1 5.3 4.6 4.8 5.1 3.6 3.3 4.8 4.6 4.5 3.9 5.7 4.1 4.0 4.4 4.2
##  [253] 5.5 2.6 3.9 4.7 5.3 4.4 5.0 4.9 2.3 2.3 4.8 4.3 4.1 4.6 6.0 4.5 5.6 4.1
##  [271] 5.4 2.6 6.4 3.3 4.8 4.8 4.9 5.4 2.0 3.9 4.1 4.9 4.8 5.1 3.9 4.5 3.3 4.4
##  [289] 4.8 5.1 5.1 3.7 4.8 5.0 5.5 4.9 5.2 4.5 5.4 4.6 3.3 2.8 3.4 4.9 5.8 4.4
##  [307] 4.9 4.0 4.7 3.4 5.1 4.9 4.6 5.2 3.9 3.3 3.8 3.9 4.9 4.3 5.2 4.4 5.1 4.9
##  [325] 3.6 5.5 5.1 3.2 6.3 4.0 3.9 4.0 3.4 3.7 3.7 2.6 6.0 7.2 5.4 4.7 6.3 4.0
##  [343] 5.3 5.0 4.6 2.8 5.3 4.6 5.1 3.9 3.8 4.7 5.1 3.0 4.1 3.6 2.9 4.0 5.6 3.7
##  [361] 4.0 4.3 3.4 2.9 4.0 4.7 5.4 5.4 3.9 3.3 4.8 4.5 4.0 4.2 4.1 5.9 4.7 3.6
##  [379] 5.6 5.7 3.6 4.5 4.0 5.7 2.9 5.5 5.8 2.9 3.8 5.7 5.2 3.6 5.2 4.2 4.1 5.3
##  [397] 4.2 4.0 4.4 4.1 6.3 5.0 4.8 2.5 5.0 5.0 3.9 5.0 4.3 4.2 5.1 4.7 5.5 6.3
##  [415] 4.7 3.6 4.8 4.5 4.8 3.4 4.0 3.3 3.7 4.6 4.1 4.0 2.9 4.2 4.0 4.9 5.4 3.6
##  [433] 3.7 4.9 4.0 5.0 3.4 3.1 4.5 3.2 6.0 6.3 3.8 4.9 4.7 3.7 3.9 5.6 5.3 5.2
##  [451] 2.5 4.8 3.3 4.1 5.0 4.0 2.8 3.2 3.8 3.8 5.3 6.2 3.6 4.8 3.7 3.8 4.9 4.7
##  [469] 4.0 5.8 3.9 5.5 3.0 3.3 3.2 4.1 6.0 4.8 4.5 4.8 5.8 3.9 4.4 3.7 4.9 5.3
##  [487] 5.4 5.8 3.9 6.5 4.8 5.2 4.6 5.4 3.9 6.9 3.8 5.0 4.7 3.5 3.2 5.1 4.0 4.2
##  [505] 4.5 2.7 3.9 4.0 5.1 5.5 4.2 4.7 5.1 4.7 4.8 3.0 5.3 3.5 3.4 5.1 3.0 3.0
##  [523] 5.4 3.7 2.1 4.0 4.2 3.0 5.5 6.3 4.4 4.6 4.1 4.1 5.9 4.4 6.1 4.9 5.3 3.5
##  [541] 3.8 4.8 4.5 3.7 4.4 5.7 4.7 5.1 3.4 5.9 4.2 3.7 4.9 3.9 5.8 4.2 5.3 4.7
##  [559] 4.3 6.2 5.6 5.5 3.5 3.4 4.7 3.6 4.0 5.1 4.0 5.1 4.0 4.3 2.6 4.4 5.1 4.1
##  [577] 4.0 3.3 4.7 3.9 5.9 4.9 5.3 4.1 4.2 4.6 5.3 2.5 5.3 6.4 4.0 4.6 5.4 4.5
##  [595] 3.9 2.8 4.2 3.6 4.9 5.6 4.3 3.0 3.5 4.7 3.6 3.3 4.5 5.4 6.4 4.6 3.1 5.6
##  [613] 5.6 4.0 5.4 5.6 4.2 4.7 4.9 3.3 4.4 3.4 3.8 3.7 5.8 4.4 5.5 3.7 2.3 5.2
##  [631] 5.0 4.3 3.8 3.9 5.6 5.3 4.4 3.2 3.9 5.3 5.0 4.8 4.5 4.3 4.2 4.4 4.3 5.2
##  [649] 5.1 4.4 4.1 4.3 4.9 3.4 3.5 4.3 3.3 3.8 5.1 4.2 4.7 4.6 6.6 4.2 5.8 4.2
##  [667] 2.4 4.3 4.2 4.4 4.4 4.9 3.0 2.6 3.4 5.9 5.1 4.7 4.2 5.1 5.4 4.7 4.1 3.8
##  [685] 2.6 3.3 3.9 5.8 4.6 5.0 5.5 6.1 3.0 5.0 4.5 3.2 3.8 6.0 5.8 4.1 4.9 3.9
##  [703] 3.6 2.2 3.7 3.0 5.6 4.5 5.1 3.9 5.7 4.6 3.4 5.5 4.0 3.6 5.0 4.3 2.9 5.1
##  [721] 4.8 5.0 4.8 4.7 5.1 3.9 6.2 4.8 3.5 5.0 4.7 4.3 4.8 4.5 5.3 4.4 6.4 2.8
##  [739] 3.6 5.6 4.4 3.9 4.7 5.3 3.9 4.9 4.4 4.3 4.7 3.6 3.4 4.4 4.8 5.2 5.0 4.1
##  [757] 5.1 5.8 4.6 4.3 4.5 4.7 4.8 3.6 4.5 3.1 5.5 5.1 6.3 4.4 5.2 4.6 2.6 4.8
##  [775] 5.9 4.3 5.3 4.8 4.4 3.4 4.8 3.7 4.8 4.5 5.2 6.5 4.8 3.8 3.9 3.0 3.5 4.2
##  [793] 4.2 3.5 5.2 5.0 4.3 3.2 4.5 5.7 3.7 4.1 4.1 6.4 4.5 4.3 3.8 4.4 4.0 5.0
##  [811] 5.3 5.9 3.9 5.0 3.6 4.9 5.1 3.6 5.6 3.2 4.1 5.8 4.2 4.6 3.8 5.2 5.1 3.7
##  [829] 4.4 4.1 4.6 4.5 2.8 3.7 5.1 2.3 4.1 5.2 4.1 4.2 3.4 5.2 3.7 3.3 4.0 5.1
##  [847] 4.4 5.7 2.8 4.4 4.9 5.5 5.2 4.6 4.2 3.3 5.4 3.9 4.8 5.5 2.7 4.5 4.2 4.3
##  [865] 4.2 5.4 5.3 5.5 4.9 5.0 3.0 4.0 4.8 4.1 3.7 6.0 5.3 4.3 4.4 6.5 5.6 5.4
##  [883] 5.4 4.3 4.7 3.6 4.2 5.9 5.0 4.4 4.1 4.3 4.4 3.9 3.8 3.9 5.3 3.7 4.5 4.5
##  [901] 4.9 4.7 4.0 4.9 3.9 4.8 6.0 4.3 5.2 6.0 4.1 4.9 4.9 5.1 6.1 4.5 4.2 4.3
##  [919] 4.0 4.4 4.1 4.2 4.3 4.6 5.0 3.8 4.6 5.1 3.5 3.8 5.1 3.0 4.4 3.7 6.2 5.3
##  [937] 4.8 5.2 4.2 4.6 2.5 3.8 4.8 4.5 4.5 2.9 3.5 4.1 4.8 4.2 4.4 2.8 4.9 5.3
##  [955] 3.8 4.8 4.0 3.5 4.4 5.2 5.7 4.7 4.9 4.2 4.0 3.0 6.8 6.3 3.7 5.1 3.6 5.1
##  [973] 4.7 4.0 4.7 4.9 4.5 3.2 4.1 4.6 4.4 4.8 5.2 4.3 4.1 5.8 3.7 5.1 4.9 4.8
##  [991] 4.2 3.3 3.0 3.4 3.8 4.2 4.7 4.4 3.9 4.7
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
##   2.7   6.3
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
##    [1] 4.2 5.8 4.3 3.7 2.9 4.5 5.5 4.5 3.2 4.8 4.8 2.0 3.7 4.8 4.4 4.3 3.7 3.5
##   [19] 4.2 4.7 5.1 4.9 4.4 3.5 4.4 6.5 4.7 3.1 3.1 5.0 4.2 5.3 4.0 4.4 3.8 4.1
##   [37] 5.5 3.6 3.3 3.7 3.5 3.7 3.8 5.2 6.7 4.0 6.2 4.7 5.1 3.6 5.0 4.6 4.8 3.5
##   [55] 4.7 3.1 2.2 4.5 5.4 4.5 4.5 5.6 5.0 3.8 3.6 3.9 4.7 4.0 5.3 4.9 4.8 3.2
##   [73] 4.0 5.1 4.9 5.7 5.6 5.6 5.1 5.1 5.1 4.0 5.1 4.8 4.2 3.4 3.1 3.0 4.2 4.2
##   [91] 6.5 3.4 4.7 3.3 4.2 5.8 3.4 4.6 3.3 4.4 3.6 2.7 4.2 3.9 4.6 6.5 4.6 3.2
##  [109] 5.3 4.9 2.9 4.3 2.8 3.5 5.8 4.3 5.3 4.3 4.4 4.5 4.7 4.6 3.7 4.5 4.0 3.6
##  [127] 3.5 4.8 4.5 3.7 2.8 5.3 5.2 4.7 4.8 4.2 6.9 6.2 4.6 3.8 4.8 3.8 5.4 3.7
##  [145] 3.4 4.8 4.7 3.0 4.4 4.4 5.4 4.1 6.3 4.4 4.5 4.8 5.4 2.7 5.3 6.4 3.4 3.5
##  [163] 3.6 6.2 3.5 3.6 4.5 4.5 4.9 3.6 4.8 4.4 5.3 4.5 5.8 4.6 6.2 4.2 4.0 4.6
##  [181] 4.3 3.4 3.1 4.3 3.8 5.1 4.1 4.0 4.6 3.9 3.4 5.7 4.6 3.8 4.8 4.4 5.6 4.1
##  [199] 4.8 4.5 3.8 3.5 4.1 2.9 4.3 5.2 4.2 3.2 5.5 4.2 3.5 4.8 3.7 4.1 4.2 6.4
##  [217] 4.7 4.2 4.2 3.9 5.0 3.4 3.2 3.9 3.4 5.5 4.8 4.5 6.5 3.9 3.0 2.4 4.9 3.7
##  [235] 4.8 3.9 5.0 3.7 3.2 6.3 5.3 3.3 5.4 3.6 4.1 4.0 3.9 3.0 4.7 2.9 3.6 4.1
##  [253] 3.5 5.8 3.3 5.7 3.8 7.1 4.8 5.0 4.7 3.8 4.3 5.4 3.2 3.0 6.1 3.8 3.9 5.5
##  [271] 6.0 5.6 4.8 4.6 3.7 4.8 5.6 3.7 3.3 5.0 4.7 4.7 4.3 4.4 2.9 4.8 3.4 4.0
##  [289] 5.2 3.4 5.0 4.8 4.0 4.9 4.1 4.7 3.7 4.4 4.5 5.9 5.9 5.1 5.3 3.2 6.2 5.2
##  [307] 4.1 6.0 3.6 3.5 4.4 5.2 4.9 3.7 4.8 3.9 4.6 6.3 3.4 4.5 4.5 5.9 5.3 3.0
##  [325] 3.4 3.4 4.1 3.3 2.7 3.1 5.0 4.0 4.5 4.0 4.5 4.1 4.7 5.3 3.0 4.2 5.8 4.3
##  [343] 6.2 3.7 5.1 4.2 4.5 4.7 4.4 4.6 3.8 5.6 5.9 5.4 4.5 4.5 5.0 4.4 4.8 3.1
##  [361] 5.4 4.7 4.8 5.8 4.5 5.8 4.2 5.1 5.2 5.1 4.1 5.3 4.4 4.3 5.4 5.3 5.3 3.0
##  [379] 4.2 4.0 5.7 3.9 3.2 5.3 5.7 5.4 5.1 4.4 5.2 3.0 5.4 3.3 5.5 4.1 4.8 4.2
##  [397] 4.8 4.3 3.4 4.6 4.4 5.0 4.0 4.8 4.5 3.7 5.1 4.8 5.4 4.1 3.7 2.9 5.6 5.1
##  [415] 4.3 4.8 5.1 6.2 4.5 3.2 5.0 4.8 3.6 4.5 4.8 4.6 4.5 5.8 4.4 3.2 4.6 4.4
##  [433] 3.8 3.8 6.0 4.9 4.6 5.3 4.2 4.1 6.5 5.6 4.2 4.8 4.9 5.4 4.7 3.6 3.4 4.8
##  [451] 5.7 4.8 5.3 3.6 4.0 5.9 3.9 4.8 4.1 4.2 5.5 6.1 4.1 5.0 5.8 4.6 5.9 4.0
##  [469] 4.6 5.0 3.7 4.8 2.9 4.3 3.9 6.4 4.0 5.2 4.5 4.6 4.8 4.5 3.5 4.8 4.5 5.0
##  [487] 4.3 5.4 5.8 3.2 3.9 4.5 4.0 3.8 2.6 5.4 4.6 4.9 4.5 4.7 5.0 6.3 5.1 5.0
##  [505] 4.7 3.7 4.7 6.2 4.2 4.9 4.0 5.3 4.9 2.9 3.9 4.0 4.0 3.9 6.3 4.2 4.6 4.2
##  [523] 5.4 5.1 5.5 4.1 5.2 5.3 3.8 4.8 6.3 4.4 2.5 3.5 4.7 5.9 4.4 5.0 4.8 3.3
##  [541] 4.1 3.9 4.5 3.6 6.0 6.7 4.4 5.2 3.3 3.9 4.5 6.1 4.6 4.7 4.8 4.2 4.7 5.3
##  [559] 2.1 4.7 5.9 3.1 4.6 4.0 4.5 5.2 4.9 4.4 6.1 5.4 3.3 6.2 5.9 4.9 3.4 6.1
##  [577] 3.6 6.2 3.6 2.9 5.7 2.3 5.4 4.7 5.6 4.5 5.4 4.5 6.3 3.1 4.1 4.8 5.7 4.4
##  [595] 5.2 4.9 4.9 3.7 4.8 4.3 4.4 3.7 3.2 4.3 5.2 6.4 6.4 4.7 3.7 2.0 4.5 3.4
##  [613] 2.9 5.3 5.5 6.7 3.7 5.7 6.0 5.7 5.7 4.4 2.9 4.8 4.6 2.5 5.6 4.1 4.4 5.1
##  [631] 4.7 3.7 2.4 4.2 3.8 5.7 5.4 3.3 3.0 4.9 4.1 3.7 4.0 5.4 5.5 3.8 6.3 3.3
##  [649] 4.9 3.8 3.1 4.3 5.1 4.3 4.5 5.2 4.1 3.7 4.7 4.8 3.2 3.1 4.5 6.3 5.0 4.6
##  [667] 4.6 4.1 5.4 5.6 6.1 4.1 6.3 4.8 2.6 4.5 5.2 5.3 4.5 4.1 5.6 3.4 4.3 4.0
##  [685] 3.6 3.6 3.3 6.0 2.9 5.5 2.9 5.6 4.8 5.0 4.2 2.7 3.8 3.6 2.7 4.5 4.9 4.7
##  [703] 4.5 3.3 4.8 5.3 5.3 4.2 5.5 4.3 4.4 5.0 5.0 4.2 3.8 4.3 5.3 4.6 4.8 4.2
##  [721] 5.6 3.2 5.1 3.8 4.2 6.9 7.0 4.2 5.9 6.1 5.9 3.5 2.6 5.0 3.4 3.8 4.0 5.5
##  [739] 4.5 5.6 4.4 3.1 4.6 4.1 4.8 4.8 4.9 4.2 5.0 4.9 5.0 3.6 3.4 4.8 3.8 5.0
##  [757] 3.8 5.7 4.0 5.2 5.0 4.3 4.2 4.2 6.1 3.9 4.9 4.8 6.0 3.8 5.7 3.7 5.2 2.1
##  [775] 3.7 4.9 5.5 3.5 4.8 3.9 4.9 4.0 3.9 5.2 3.0 5.0 3.2 2.7 4.3 5.4 6.3 5.4
##  [793] 3.1 4.3 4.7 4.5 4.2 2.7 3.8 5.7 4.2 4.4 5.6 4.7 4.9 5.4 4.5 3.1 4.5 4.4
##  [811] 5.2 3.1 4.0 4.2 4.9 2.5 4.3 5.2 4.7 4.9 3.8 4.4 5.5 6.0 5.3 2.7 5.3 6.4
##  [829] 3.7 5.6 4.5 5.3 4.2 4.7 5.8 5.2 5.0 3.9 4.9 4.7 4.3 4.4 4.7 5.8 5.3 3.5
##  [847] 4.3 4.1 6.2 3.5 3.9 3.2 3.6 4.4 4.2 5.1 7.1 5.8 2.6 5.5 3.9 3.5 4.3 6.2
##  [865] 3.3 4.1 2.5 4.1 5.5 5.7 5.5 5.7 4.9 4.5 4.5 5.0 2.8 5.2 3.9 5.5 3.9 5.0
##  [883] 5.2 5.9 5.2 3.8 4.3 6.3 4.5 4.6 7.0 5.2 4.5 5.5 3.1 4.9 5.8 4.6 6.8 5.6
##  [901] 4.8 4.7 5.3 5.3 4.5 5.1 5.0 3.0 4.0 5.2 3.3 4.9 5.1 4.7 4.6 6.6 4.7 3.2
##  [919] 4.9 4.4 3.7 4.7 4.7 6.3 5.2 4.6 4.3 4.5 5.4 4.2 2.3 4.1 6.1 5.7 2.9 4.6
##  [937] 3.2 4.0 5.3 5.1 6.0 5.2 5.5 6.1 3.9 4.8 5.0 4.5 5.4 5.5 4.2 3.8 5.2 4.8
##  [955] 5.9 4.5 4.6 5.5 5.6 4.9 4.0 4.8 6.2 3.4 3.7 4.8 2.7 4.3 5.0 5.5 5.4 3.0
##  [973] 4.5 5.7 5.4 5.1 4.0 5.1 4.2 3.2 3.4 4.5 4.6 4.0 5.0 4.7 5.1 2.8 2.4 3.1
##  [991] 4.6 3.9 3.7 5.3 3.8 5.8 4.6 5.8 5.8 4.3
## 
## $func.thetastar
## [1] 0.0283
## 
## $jack.boot.val
##  [1]  0.51805930  0.45212121  0.30311526  0.19318801  0.06657825 -0.10773994
##  [7] -0.18534031 -0.29532164 -0.34849398 -0.47710145
## 
## $jack.boot.se
## [1] 0.9892098
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
##    [1] 4.6 4.2 5.0 5.0 4.7 4.8 3.4 4.3 4.9 2.8 5.9 4.0 6.8 5.0 4.2 4.3 4.2 4.7
##   [19] 4.8 3.6 5.1 5.8 4.0 4.6 4.6 6.7 6.0 5.0 4.1 5.6 5.2 6.2 4.6 3.5 4.0 5.4
##   [37] 4.8 5.1 5.3 3.4 4.0 6.0 5.0 2.3 3.0 5.4 4.3 5.9 3.6 4.9 4.6 3.9 3.0 3.0
##   [55] 5.9 5.2 5.0 4.2 4.8 5.7 3.4 4.9 2.9 4.6 4.4 3.4 5.5 4.3 4.7 3.8 3.2 4.7
##   [73] 5.0 5.2 4.2 4.1 3.2 5.3 5.2 4.2 5.0 3.9 4.2 4.3 3.0 3.5 3.9 4.5 4.8 5.9
##   [91] 5.7 3.8 6.8 5.4 4.1 5.5 6.3 3.0 3.2 4.5 4.0 6.2 4.1 5.2 3.9 3.7 4.2 4.5
##  [109] 3.7 6.4 3.6 4.2 5.2 4.8 4.0 4.5 3.8 4.8 3.7 4.5 5.9 4.1 5.1 3.7 2.2 6.2
##  [127] 5.2 4.0 3.9 4.6 3.3 2.2 3.9 4.7 4.6 5.7 4.5 5.4 5.2 5.0 5.7 4.9 5.2 4.7
##  [145] 5.7 4.9 4.7 4.1 5.1 5.8 2.7 6.0 5.3 5.6 5.9 4.5 5.5 4.1 4.8 4.9 3.0 4.1
##  [163] 4.6 4.7 3.5 3.4 4.7 3.3 5.1 5.2 4.6 5.0 4.1 4.5 4.2 5.4 3.9 4.5 4.2 6.4
##  [181] 5.3 5.0 3.9 5.6 5.5 2.7 4.7 3.6 3.7 5.5 3.9 4.5 5.1 3.1 4.7 4.3 6.1 4.4
##  [199] 3.9 4.0 3.7 3.5 5.5 3.8 4.3 3.4 3.8 5.4 3.4 5.5 3.3 3.7 3.5 4.0 4.3 4.5
##  [217] 4.2 3.4 4.4 5.2 4.0 5.5 3.4 4.0 2.2 5.6 4.9 3.5 4.7 5.8 5.3 5.0 4.6 4.1
##  [235] 4.4 4.8 4.7 4.2 5.6 4.8 4.3 3.4 4.1 3.9 5.6 5.0 6.7 4.6 4.2 3.3 3.8 5.2
##  [253] 3.8 5.1 5.9 5.0 3.4 4.9 4.7 5.3 4.4 5.9 4.5 6.1 4.7 4.9 4.7 4.6 2.5 3.5
##  [271] 4.5 5.7 3.8 5.1 3.7 4.9 3.4 7.4 4.2 4.0 4.0 5.1 4.1 6.0 4.0 5.2 3.9 2.8
##  [289] 3.7 3.0 3.7 4.7 4.7 3.4 3.8 3.3 4.8 5.2 4.6 4.6 4.8 5.0 3.2 4.5 4.5 3.2
##  [307] 4.0 4.4 5.3 5.2 5.2 4.3 4.7 4.4 3.3 5.2 2.8 6.7 3.4 5.5 4.7 5.7 3.5 4.1
##  [325] 2.0 4.6 3.4 3.5 3.5 4.5 5.0 4.4 4.2 4.1 4.8 5.1 4.8 5.2 4.4 3.6 3.5 5.2
##  [343] 4.3 4.1 4.3 4.1 5.6 4.1 3.5 4.4 5.5 2.8 4.4 5.2 3.1 4.4 4.1 3.2 3.8 5.1
##  [361] 2.8 3.9 4.1 4.6 5.1 4.0 5.0 3.7 5.3 4.2 5.8 6.6 5.8 5.7 4.7 3.2 5.1 4.2
##  [379] 6.0 4.1 5.7 4.5 5.6 5.1 4.7 3.4 4.7 6.5 3.3 4.3 2.9 4.0 4.7 4.2 4.8 5.0
##  [397] 5.4 3.3 5.1 4.7 4.7 5.4 5.0 5.3 4.6 5.1 5.9 4.2 5.0 5.9 3.3 4.4 4.1 3.5
##  [415] 3.7 2.9 4.7 3.1 5.2 4.3 3.7 4.5 5.2 6.4 5.1 4.5 5.3 4.6 5.5 3.3 5.1 5.1
##  [433] 5.4 3.1 4.8 4.6 2.7 3.4 4.4 6.3 4.5 3.4 4.2 4.6 3.5 3.5 3.2 2.9 5.7 3.3
##  [451] 5.8 6.4 4.5 5.0 4.9 6.0 5.2 4.4 3.2 5.3 4.5 3.5 4.2 6.2 5.2 5.3 4.8 3.9
##  [469] 4.7 4.5 4.0 6.0 5.0 5.2 5.3 4.7 4.0 3.7 3.4 4.5 4.3 5.1 5.9 4.4 3.9 3.0
##  [487] 4.3 2.5 5.2 5.9 5.0 4.1 6.1 3.3 5.5 4.1 3.3 3.9 4.9 2.9 1.9 3.6 1.9 4.2
##  [505] 5.3 4.4 5.6 4.5 2.7 3.3 3.9 5.1 5.5 3.6 3.6 4.1 3.4 4.4 4.2 4.8 4.4 4.3
##  [523] 4.1 4.7 4.4 3.9 4.8 5.4 5.0 4.5 4.6 5.5 3.1 5.5 3.2 3.2 4.8 4.0 5.9 4.4
##  [541] 5.3 3.4 3.9 3.3 4.2 4.2 4.7 4.1 5.1 4.3 4.6 5.3 5.7 4.8 4.1 2.5 3.4 4.9
##  [559] 4.0 3.6 4.3 2.3 4.0 4.3 6.8 4.5 6.0 4.3 5.2 5.1 4.0 3.6 4.2 4.4 4.1 4.6
##  [577] 5.1 3.8 5.0 6.2 3.6 4.0 4.7 5.4 3.2 3.9 3.9 5.0 4.2 2.2 4.3 5.6 4.1 3.1
##  [595] 4.0 4.3 3.7 2.4 3.7 4.2 6.1 5.1 5.8 5.0 4.5 3.7 4.4 4.4 5.7 3.8 4.9 4.0
##  [613] 4.1 4.1 3.8 4.8 3.2 4.2 4.5 4.1 4.3 3.7 4.9 5.2 3.0 5.2 3.9 4.8 5.1 5.9
##  [631] 5.4 3.3 5.2 4.2 4.9 5.0 5.2 3.9 4.4 6.5 5.2 4.7 5.6 3.5 3.8 4.0 4.7 5.0
##  [649] 4.5 4.5 3.8 5.5 4.1 4.8 6.1 5.6 2.9 4.6 4.4 3.1 4.0 5.7 3.7 5.6 4.8 6.7
##  [667] 3.2 5.8 4.2 4.0 4.5 5.6 4.1 4.3 4.6 4.3 5.3 2.5 3.8 5.0 5.3 4.4 3.0 4.1
##  [685] 4.0 4.6 4.1 4.3 3.6 4.5 4.1 3.6 5.3 3.2 5.9 3.4 4.0 3.9 3.7 4.9 6.5 3.7
##  [703] 5.0 4.4 4.7 4.0 5.4 3.7 6.0 5.8 5.0 4.9 3.6 4.6 4.6 4.6 4.5 3.7 4.3 4.3
##  [721] 4.7 4.8 3.3 4.4 4.8 5.2 5.4 4.4 3.4 4.6 6.7 3.3 4.3 3.5 3.0 3.7 4.2 4.0
##  [739] 4.5 2.6 3.1 3.3 5.4 3.7 4.8 4.1 4.9 2.8 3.6 4.0 3.5 4.8 4.1 3.9 5.5 3.4
##  [757] 3.0 4.9 4.3 4.6 3.1 4.0 6.7 4.2 4.4 4.7 5.5 3.4 5.4 5.7 4.0 5.8 3.2 5.4
##  [775] 5.7 4.6 4.2 6.4 4.8 4.8 4.4 3.6 4.9 2.7 3.7 5.5 4.0 4.2 4.0 4.6 5.1 5.5
##  [793] 1.8 5.1 5.6 5.3 4.3 5.7 5.6 4.1 3.5 4.8 3.8 2.9 5.9 4.8 5.3 6.1 4.3 5.3
##  [811] 4.8 3.9 5.1 5.1 4.5 4.6 3.7 3.9 5.3 5.3 6.0 3.9 4.0 5.1 6.1 3.5 5.4 3.4
##  [829] 3.9 5.1 4.6 5.2 4.8 4.4 2.8 3.7 4.3 5.5 6.1 5.3 5.9 5.6 3.0 4.5 5.5 4.3
##  [847] 4.5 5.4 4.5 4.9 4.8 3.6 5.7 4.4 4.1 4.0 5.3 5.7 4.4 4.2 4.0 2.9 4.4 2.7
##  [865] 4.6 3.4 4.4 4.6 4.1 5.5 5.3 6.6 3.8 3.5 4.2 5.3 2.1 3.4 5.9 4.9 3.8 3.6
##  [883] 4.3 3.9 5.0 4.7 4.6 4.8 3.5 6.2 5.7 4.1 3.4 5.0 5.8 5.6 5.1 5.9 4.6 3.8
##  [901] 4.8 5.8 5.4 5.1 3.6 4.3 2.1 4.5 3.8 4.8 4.9 3.6 5.0 4.7 4.7 6.2 4.9 4.4
##  [919] 4.6 5.7 4.3 6.9 3.6 4.0 4.5 4.6 2.2 3.5 5.6 4.7 4.4 5.7 5.7 4.3 5.4 4.5
##  [937] 4.7 4.8 4.5 4.9 4.7 5.2 4.8 5.3 5.1 4.1 5.0 3.6 4.1 3.9 5.8 4.1 4.3 4.9
##  [955] 3.5 3.6 3.6 5.2 4.9 4.3 5.4 5.5 4.7 3.0 4.6 5.1 5.7 3.9 4.7 4.5 3.0 4.6
##  [973] 5.1 4.0 6.4 4.7 7.0 4.1 3.3 3.3 6.6 6.4 5.2 5.6 4.1 3.9 4.0 5.0 4.4 5.1
##  [991] 3.7 5.4 3.5 4.7 4.2 3.9 4.7 4.8 5.7 3.9
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.5 5.3 5.3 5.2 5.1 5.0 4.9 4.7 4.7 4.5
## 
## $jack.boot.se
## [1] 0.9079648
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
## [1] 1.741258
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
##   3.228422   4.686353 
##  (1.375735) (2.160625)
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
## [1]  0.7799487 -0.2035393  0.5602398  0.2753324  0.6859434 -0.4306753
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
##    [1]  0.6388187208  0.4195018083 -0.0183413250  1.1885552918  1.8671190782
##    [6]  2.0011133317  1.9912914447  1.7292943393 -0.3300519078 -0.0829291132
##   [11]  1.2203485940  0.2704584968 -0.6373152544  1.5878722439  1.6675519738
##   [16]  0.9392238012  2.2114876129  2.2605174027 -0.8955398981  0.5826981294
##   [21] -0.5902674639  1.0932928771  0.3404430236  1.9161570862  2.3945432964
##   [26] -0.4391570211 -0.1632780168  1.1722621218  2.0261044354  1.0618802189
##   [31] -0.7802346317 -0.8336441689  2.0852879657  2.3687148062 -0.5403388625
##   [36]  1.0573839665 -0.0560335644  0.1797586475  1.8014108437  1.3945309740
##   [41]  1.3485247501  0.3797163651  0.8660793400 -0.8433555298 -0.4157140709
##   [46]  2.2374315924  1.4453013162  0.3688566482 -0.3763317156  1.5044459348
##   [51] -0.3774376312 -0.2498375400  0.7212307489  1.2122450319 -2.0000160439
##   [56] -0.2699542341  1.1799496223  1.7293760350  0.2366084109  0.0865094656
##   [61]  1.1973361383  0.0969595947 -0.2870425567  0.4069135210  1.4127586773
##   [66]  1.3567140727 -0.6394466256  1.9839907847  1.5248175947  1.3417638972
##   [71]  1.6801583469  0.3144198730  1.1148070777  0.0898002432  1.2970871606
##   [76]  1.1466307139  0.7080708414  0.4316587503 -0.0159387478  1.1650739193
##   [81]  2.2730609902  0.7929323031 -0.8433555298  2.2900301322  0.9888364829
##   [86]  2.0507055934  2.0712967122  2.2763186791  0.7125383649  1.2586531173
##   [91] -0.5056887187  1.2073788807 -0.8406866289  2.1156959363 -0.7621946121
##   [96]  1.3606783704  1.7662213931  0.1365626741  1.6460454121  1.6213088004
##  [101]  1.2267511079  1.9956040706  1.1705236220 -0.4994649809  1.6870907750
##  [106]  1.2138409036  1.9082346416  1.5216432758  0.4958542865 -1.2519770980
##  [111]  0.5329038088  0.6842533202  0.6900340235  1.7290180678  1.0731880049
##  [116]  1.8715397947 -0.0973799883  1.6034161650 -0.8280394112  0.3144198730
##  [121]  1.2852021298  0.7189886729  1.3371963859  1.8523670480 -0.5920823256
##  [126] -0.1098915591  1.2661576554  1.9862804681  0.6638354506  2.1607951839
##  [131]  1.4162849725  2.0834350673  0.1356846836 -0.4958933067 -0.5727564590
##  [136]  2.0920979750  1.1659615444  1.7288235635  1.8995816785  1.9643989132
##  [141] -1.2482332010  1.8655551577 -0.5307251215  1.2178618257  1.2105300132
##  [146]  0.9700654111 -1.7455769581  1.1808558831  0.7530335029  1.2440181988
##  [151]  1.6376755685  1.1694541516  0.7438457705 -0.4581992566  1.9408844008
##  [156]  2.1101668535  0.9401432029  1.0186366257 -0.0443647640  0.1205772141
##  [161]  1.6223475888  1.5555797167  1.2612673677  0.7275563111  1.8826036818
##  [166]  1.3076981205  1.0540740628  0.3266662334 -0.6037537499  1.0461161650
##  [171] -0.5016261568  1.1607808052  1.9463783378  0.1832263635  2.0368036361
##  [176] -0.1005479796  0.2243771381  0.4196226417  0.2550501587 -0.6702902687
##  [181]  2.1970022226  0.3589737264  1.8957554782  1.7425874089  1.6204517388
##  [186]  1.4591196428  1.5501875976  0.9942072732  0.8984132178  0.4308075959
##  [191]  1.8769799007 -0.8065474937  0.9872358097 -0.9314757813 -0.1079641188
##  [196]  1.2321363603  1.1921130800  1.3177507825 -0.8509671047  1.5191455269
##  [201]  1.0060891826  1.4177097940  2.0430793482  1.5733820838 -0.2527454216
##  [206]  1.4321942934  0.0618675863  1.5437414730  0.8180643996  1.1014106264
##  [211] -0.6287204648  0.5615629681  1.0730852061  1.8801091868  0.2408991119
##  [216]  1.6449034288  1.3936753056  0.1924801156  1.7256131321 -0.0974009631
##  [221]  1.2128783215  1.3737518847  2.0847572520 -0.3148004560 -0.5170646987
##  [226]  2.0270775657 -0.1086961455  0.1405052055  0.1029880744  0.7696780348
##  [231]  0.4009418503  1.2469289964  0.3023667287  1.5339050117  0.4966692826
##  [236] -0.0734279694 -0.4870638246 -0.3868082312  1.7644353291  1.0842293360
##  [241] -0.9777387441  0.2940514189  1.9473195514 -0.1937211509  1.5539081943
##  [246]  1.7561157513  1.4363605402 -0.6455546504  0.0037267152  0.1989209307
##  [251]  0.1914193025  1.6381943126  0.7775787590  0.0004463083  0.8731518053
##  [256]  1.0975733428 -0.0795617653 -0.9918459796  0.7664968364  0.9560456473
##  [261] -0.3915504789 -0.3190020252  0.6923498711 -0.3685927232  1.2702829459
##  [266]  0.8485768711  1.1400523208  1.3377868826  1.5807560429  0.4202397292
##  [271]  0.6992081152  1.6566259508  1.7376823959 -0.2200712834  1.9488087510
##  [276]  1.3300109473  1.3457274460 -0.5458292669  2.1289684343  1.0633681949
##  [281]  0.1779401481  0.4786585176  0.3569006271  1.0851737793  2.0649033396
##  [286]  1.8005683864  2.0826194821  1.2927976740  1.3282836288  1.3753361184
##  [291]  1.7886977248  2.2250000609  1.7360741546 -0.4342447169 -0.7186660384
##  [296]  1.2548546373  0.5856252092 -0.0779797449  0.0398591647 -0.0560786253
##  [301] -0.2215377546  0.4294712661 -0.9727027343  0.9560912069  2.3399746748
##  [306]  1.1645103185  1.4460467876 -0.7777392624  1.9060697552  0.3039993060
##  [311] -0.8870527532  1.7229181907  1.1291251486  0.0948952901  1.9204341363
##  [316] -0.3641312712 -0.1836631752  1.6408771903  0.9645374079  1.8069300097
##  [321]  1.6884545709  1.1898334290  1.2552786180  0.3182570025 -0.3072938647
##  [326]  1.8641777220  1.0924025227  1.3810482798  1.2944114020  1.1590342036
##  [331]  0.6928067958  1.2029845871 -0.3032559790  1.3341687865  0.0457427111
##  [336]  1.6151208187  0.8480771790  2.0583591732  2.2156005156 -0.3623347966
##  [341]  0.5605510522  0.4070322920  1.2141613151  1.2291141357  0.1920037714
##  [346] -0.4861461995 -0.5082417173  1.3002821689  1.1346272895  0.5943877352
##  [351]  1.7207287252  0.9840625036  0.3323937600  1.3268155851  1.8299489078
##  [356]  1.7946866313  0.1112915511  1.2648947146  0.0850791780 -0.8266237246
##  [361] -1.3420073773 -0.3508277478 -0.5286600870 -0.9331087227  1.7281665953
##  [366]  1.1789465451  1.9890687574  1.2982926246 -1.1468315151  1.6788911016
##  [371]  1.6684182198  1.6544337933  0.3144198730  1.1073109542  1.3202476998
##  [376]  0.7287113945  2.0918462729 -0.4874056200  0.1812705305  0.0960429803
##  [381]  0.0966859443  0.0847905424  0.8001418203 -0.2225954294 -0.0314413710
##  [386]  1.4407229728  0.7217372320  0.7076415959  1.4149669644  1.9430304117
##  [391] -0.4317088528  0.4911441809 -0.3961607636 -0.1547739241  1.8973287767
##  [396] -1.0080210737  1.1051776181  0.0696263900  1.3746189206  2.3456045449
##  [401] -0.6080490923 -0.4229824828  0.0887935355  1.5668225507  1.3141322981
##  [406]  0.2476615542  1.9801541087  1.4660394385  0.1205772141  1.6624880222
##  [411]  2.3454882414  1.5740232465  1.5579301101 -0.3282410571  0.5804031311
##  [416]  1.1304863147  1.8218584756  2.2589696662  0.9624706164  1.3850463961
##  [421]  1.2119516951 -0.1844518870 -0.5996089285  1.8250174028  0.6863412672
##  [426]  2.0002797924  1.8777344680  1.5337366594 -0.1056732246 -0.8047395969
##  [431] -0.6129942178 -0.4440565205  1.2691093202  1.1005536029  1.9779096994
##  [436]  1.7804722277  0.7211365667  1.1922014566 -1.0742317363 -0.3740830701
##  [441]  1.5239143162 -1.0872931134 -1.2087407764  1.9201227149  1.2521722792
##  [446]  1.1350925718  1.3425691803  1.7710790123  1.1878694280  1.9712878597
##  [451]  1.9174384564  0.9329682590  0.2010904535  2.0572991447  1.2768806353
##  [456]  1.2548956663 -0.2167358666 -0.4818251513  1.1346158484 -0.1186846855
##  [461] -0.4028518193  1.0488964340  2.2699383535 -0.3977384766 -0.4187361936
##  [466]  1.0517260239  2.1010797914  1.5341446123  2.0249878954  0.9718987580
##  [471]  2.4078247729 -0.3843578248  0.0122933023 -0.4783842247 -0.6479988371
##  [476]  0.2692374320 -0.3017032328 -0.5727564590  1.1634158577 -1.3256651703
##  [481] -0.2759142108  2.0087065705  1.6927409108 -0.2105733934  2.1159190876
##  [486]  1.2034065832  1.2287546636 -0.4189908533  1.2408052739 -0.0002769493
##  [491]  0.6824489788 -0.1769999493  2.1783643794 -0.4950908505  1.1609620882
##  [496]  1.7407949496  1.3363922431  2.1243128799  0.9843600767  1.7654315978
##  [501]  2.3241232746  1.3133868977  1.2402571427  1.0041798777  1.0547125699
##  [506]  2.4014664758  1.9964828125  1.1918505272  1.3974083149 -0.2401654507
##  [511]  0.9375030429  0.2032359430 -0.3058558178 -0.1072683297  1.7306545712
##  [516]  1.7965780586  0.9905294817  1.0150133919  0.6629942846 -0.4371272088
##  [521]  1.7705190900  1.6930523370  1.3594980967  1.8064077546 -0.3756671087
##  [526] -0.3496667042  1.3400882399  1.7863156686  1.8312210057  1.1700940070
##  [531]  1.9581205122  1.3965664364  1.2762758263  1.1031033832  1.2171325104
##  [536]  0.9312593775  1.0573480472  0.5388150164  1.3099136609  0.5338889799
##  [541]  0.6578588901  1.6732086542  0.2756754843  2.3453812980 -0.2658301024
##  [546]  1.0430536413  1.1130977065  1.9423654063  1.1882333051  1.9219020568
##  [551]  0.8796998095  0.7168071196  1.9838812382  1.1804967158  1.7965780586
##  [556]  2.3614888204  0.9713421612  1.5460080492  0.0557772169  1.2738006188
##  [561]  2.0134695691  1.5327749425 -0.4580804625  1.1987319005  0.6442402191
##  [566]  0.4509149624 -0.6581749960  1.4024420891  1.6667259580  0.6721657809
##  [571]  0.1899483277  1.1728771204 -0.9320853156  1.1365798755  0.6222063966
##  [576]  2.2838185680  1.8103103803  2.1885949248  2.0889527557  1.1675242590
##  [581]  0.7887925759  2.2748781887  1.2677915812  1.8987210085  2.0765296615
##  [586]  1.8299489078  0.3197225607  1.9309656469  2.2937815222  0.2525800284
##  [591] -0.2656086237  1.8783367394  0.7759499611  1.1010227924 -0.0622119992
##  [596]  1.3605816772  2.0638690945  0.9532424505  0.3147431565  0.6778528262
##  [601]  2.0370446098  1.3372487755  0.9715734266  2.2919951783 -0.2831849593
##  [606]  1.8131050565  1.8897187921  1.8076124910  2.0048065997  1.7998847441
##  [611]  2.2062179474  1.1249243155  1.3397495561 -0.7185552763  1.5549056158
##  [616] -0.3448590307 -0.2561326184  1.1258772291 -0.1986924429  0.6898928626
##  [621]  1.0388626442  1.7581628869  0.4875158665  1.5803539261 -0.6268779973
##  [626]  1.8752243174  1.9123424744  2.2443211904 -1.1468315151  1.7938284401
##  [631]  1.1519693234  0.3707835373 -0.3078386873  1.8846782358 -0.4829665376
##  [636]  1.9234307280  0.6422577586  1.5929316812  1.1945115146 -0.5044248024
##  [641] -0.2523656597  1.4962186033  2.3498924120  1.1159521859  1.4898229563
##  [646]  1.9414720357  1.5272535609  0.2301079919  2.2829875816  0.0926415060
##  [651]  0.3930744192 -0.3986730876  1.4770023695 -0.6723422699 -0.2646031268
##  [656]  0.2073683259  0.9932815727  0.6470279663  1.1604660756  1.5905674008
##  [661]  1.8096761480  0.7289842603  0.6237233331  0.5169819361  1.1732884399
##  [666]  0.2483343249  2.1425165788 -0.1387455466  1.1602121289  1.5268523934
##  [671]  1.3251557181  0.9868297497  2.0392769660  0.2961694763  2.4084637116
##  [676] -0.0502526397  2.0389798427  0.9023316872  1.2374432898  0.4098142692
##  [681] -1.0583489844  1.9717547402  0.6449420702 -0.6373152544  1.2151587047
##  [686] -0.2330334180  0.7112978115  1.4831194669  1.3648625676  2.0412916559
##  [691]  0.9257713995  2.0814745280  1.3309903614 -0.6186890318  1.1010865927
##  [696]  0.3404514662  1.1786541720  1.9642096805  0.0469612034  1.3339723318
##  [701] -0.2979902467  1.0610225796  2.0494136434  1.0629711193 -0.6481715822
##  [706]  0.3007300410  1.2212652637  1.9402388635 -0.2743469107  0.1168128634
##  [711]  1.6653540399  1.6097479371  0.1497028153  1.9082345925 -0.0452083782
##  [716]  0.2503114130  1.4036831279  1.5148891684 -0.7399582055  0.9862369894
##  [721]  1.3160954091 -0.3393518077  1.9892141072  1.9254538901  1.0655056125
##  [726]  1.4935094530  0.5573548861  1.8013800824  0.4638926296  1.6650547932
##  [731]  2.0832535490  1.2996048768  1.8348476314  1.1504542374  1.2466235331
##  [736] -0.3420161844  2.1160206663  1.9406070507  0.6590861710  1.2010005527
##  [741]  1.2732722028  1.3488468329 -0.2833781012  1.8852850024  1.5737425258
##  [746]  0.7321714726  2.0350505426  0.5281675740  1.9303250168  1.8887194990
##  [751]  1.8119750121 -0.1751659815  0.1623624156 -1.1481893794  2.0328831449
##  [756] -0.2496248143  1.1769436885  1.5299252904 -0.7342934234  0.9758254142
##  [761]  1.2636217217  1.7491152186  1.4498192852  1.9130583905  0.9964130043
##  [766]  1.2475760746  1.2831018563  2.0426332872 -0.5928396337 -0.1757306326
##  [771]  2.1829307561  2.0062091800  2.0469247302  2.0328455372  1.2570992471
##  [776] -0.4428006670  0.0043078804  1.6343753682 -0.5193692040 -0.2208770356
##  [781] -0.4647314899 -0.5645412879 -1.3246835075 -0.3249807905  1.9518116637
##  [786] -0.1804044650  0.4857374789  2.3077147173  0.5979784517  1.6716497651
##  [791]  0.0664957106 -0.5022308085 -1.5057620334  1.1644030960 -0.5230797127
##  [796]  1.1049226957  1.3986225545  2.0684149600  0.5422256491  1.9285936141
##  [801]  1.7396200313  1.3220459042  0.3504130243 -0.4874056200  2.0177880042
##  [806]  0.7802071659  2.3091681148  1.6706138744  2.1062056617  1.4053693651
##  [811] -0.2003486457  1.2213022140  1.0714873713  2.1133800555  2.3093276266
##  [816]  1.8773817444  1.1912311415  0.1689087948  1.6292202585  1.1853602862
##  [821]  0.0568330243  2.1742046713  1.0766776303  2.0692317010  0.4009699474
##  [826]  1.2858346548 -0.2752042330  0.4814316736  1.7626501228 -0.6899851065
##  [831]  2.1434582617  1.1098405803  1.7662213931  2.0459717107  0.5874803340
##  [836]  1.4838743645  1.7220893107  1.7708111485 -0.3608185060 -1.1135187682
##  [841]  1.7924819012 -0.5029937555  1.2045129423 -0.4083177555  1.7070983532
##  [846]  1.4807541471  1.8292204349  0.7911772081  1.0698724268  1.1798984503
##  [851]  1.3032734372  1.6122210134  0.9616665664 -0.0829608786  1.5893636384
##  [856]  0.6741773817  1.3021106362  0.9316964577  1.4329429913 -0.6301258072
##  [861]  0.6722626115  1.3163906313  1.0772778385  1.7851033557  1.9899665125
##  [866]  2.1387834267  1.6487780897  0.6928269552  1.8036057207  1.6895295740
##  [871]  2.2697481045  0.7299026868 -0.2029796365 -0.5479794381 -0.4953574360
##  [876]  0.0361320003 -1.2960400736  2.0415857278  1.8076671423  1.6316503763
##  [881]  1.8649122525  0.1652133234  2.2960486357  0.2322998144  0.2716614846
##  [886]  2.0068728260  2.0799085790  1.1853602862  1.6613611291  1.7844979964
##  [891]  1.3339807340  1.9235365949 -0.8695728692  0.9291123896  1.9524961405
##  [896]  0.2905511980  0.2254492922  1.1519087808 -0.3311350258 -0.7074487137
##  [901] -1.3245908737  1.3387222515 -0.3739802555  1.5025412757  1.4475064037
##  [906]  1.0023259406  1.6407756967  0.6097760827  0.9812183353  0.1505869837
##  [911]  1.1374291720 -0.4653463178  0.7538498208  1.6665131302  1.0484805353
##  [916]  1.2290747272  1.8024403871  0.5968140346  1.2571226798  1.1060413367
##  [921] -0.1565493630  1.0461161650 -0.0648202916  0.1002781480  1.6619514719
##  [926] -0.1158142228 -0.3808478629  0.5372643858  1.1584723771  1.2549304234
##  [931]  1.3099906396 -0.1123925089  1.5730797631 -0.6382864176  2.2542588176
##  [936]  1.7626501228  2.3022481330  1.5835938948  0.4343101019  0.0222670423
##  [941]  1.2602661130 -0.2860462455  0.9284959643  1.2106414377  1.5353031472
##  [946]  0.7325229425  1.2200635023  1.7377840260  0.6549398974  1.0457670721
##  [951]  1.5643716309  0.0941146792 -0.0380089552  1.9896307952  2.0771641651
##  [956]  0.3299979903 -0.1448205692  1.1518749053  1.2094715137  1.8169351786
##  [961]  2.3720721890  1.2139239905  1.2235472668  0.0556364829 -0.3606099284
##  [966] -0.2926057972  1.3347131452 -0.1380991394  2.2983589075  0.9926635668
##  [971]  0.2558515693  1.2998789173  0.7298951711  2.0808395329  1.2393921001
##  [976] -0.2178889629  1.5448798180  2.2848611572  0.5315447924  1.7996808519
##  [981]  1.1289661614 -0.3186707577  0.9960060766  1.1544411242 -0.6604006060
##  [986] -0.6067190534  1.1322038575  1.6641469387  1.7962089033  0.9909827590
##  [991] -0.7772800250  0.9802866788  1.6322163809  0.5320388692  1.2287546636
##  [996]  1.8836977134  1.9156524760  1.8818836642  1.8001762408  2.5142547371
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
```

```r
fit2
```

```
##       mean          sd    
##   0.68886163   0.43375781 
##  (0.13716626) (0.09698722)
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
## [1]  0.86429334  0.18001840  0.77617123 -0.14553885 -0.44152923 -0.05603099
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
## [1] 0.0168
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.870411
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
##     original      bias    std. error
## t1*      4.5 -0.04404404   0.8953241
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 1 2 4 5 6 8 9 
## 1 1 2 2 1 1 1 1
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
## [1] 0.006
```

```r
se.boot
```

```
## [1] 0.8759631
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

