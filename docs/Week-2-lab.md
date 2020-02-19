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
## 0 1 2 4 8 9 
## 1 4 2 1 1 1
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
## [1] 0.0185
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
## [1] 2.693199
```

```r
UL.boot
```

```
## [1] 6.343801
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.8000 6.3025
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
##    [1] 4.1 4.4 5.4 4.2 3.8 3.8 5.0 3.1 4.4 3.7 4.5 5.2 3.9 4.2 4.8 6.1 3.9 4.0
##   [19] 4.4 6.2 4.4 5.1 4.3 4.6 4.3 5.6 2.8 4.5 3.6 3.0 5.1 3.8 4.4 4.7 5.5 4.4
##   [37] 4.4 5.8 3.6 4.9 3.2 5.5 4.4 4.5 5.4 5.9 3.8 5.3 4.2 4.6 5.1 3.5 5.4 5.0
##   [55] 3.6 3.2 4.5 3.7 5.5 3.2 3.6 3.8 5.2 5.8 4.7 4.3 3.6 5.6 4.2 3.5 3.5 5.8
##   [73] 4.8 5.8 4.4 5.5 3.5 3.5 4.5 5.0 5.5 4.7 3.5 5.3 5.0 3.1 4.3 3.5 4.8 5.2
##   [91] 4.0 4.0 5.9 4.9 6.1 4.7 5.0 3.4 6.1 4.3 4.6 4.0 4.0 4.1 4.2 4.3 2.8 5.7
##  [109] 4.2 3.3 3.1 5.8 3.9 5.5 4.2 2.8 4.1 3.3 5.0 4.5 4.8 4.9 4.1 4.4 4.9 3.1
##  [127] 3.3 5.5 3.8 4.9 3.6 6.0 2.7 3.4 4.7 4.0 4.4 3.7 3.3 4.6 4.1 3.5 4.1 5.7
##  [145] 4.8 5.1 4.7 3.5 4.6 4.1 5.1 4.7 3.5 2.9 3.5 3.4 4.8 5.0 3.7 3.6 4.3 5.3
##  [163] 6.5 4.9 4.6 5.0 4.6 4.8 4.4 4.6 5.1 4.0 4.6 4.8 4.8 5.6 5.0 4.7 5.0 2.3
##  [181] 4.8 3.5 6.8 5.6 4.9 4.6 5.0 5.0 5.1 4.4 3.9 3.1 3.3 5.3 4.6 3.7 3.8 3.6
##  [199] 3.9 4.4 5.2 3.9 3.6 5.2 4.6 3.6 3.0 2.8 4.9 4.2 5.5 3.6 3.9 2.6 4.6 4.6
##  [217] 6.2 4.1 4.5 4.0 3.3 5.4 5.0 5.0 4.3 5.2 3.2 4.5 3.9 3.5 5.4 3.2 4.3 3.6
##  [235] 4.6 5.4 3.8 4.9 4.3 4.3 3.6 2.9 4.0 4.9 4.8 5.9 4.2 2.8 4.1 4.1 3.5 3.7
##  [253] 3.4 4.7 4.9 4.3 4.6 5.4 4.9 2.9 3.7 3.5 4.6 4.7 5.0 4.9 4.5 5.1 4.5 3.4
##  [271] 4.9 4.4 4.0 3.6 6.0 5.0 2.5 4.2 4.9 4.8 5.9 4.4 4.4 3.8 3.3 4.7 3.9 3.9
##  [289] 3.8 4.3 5.1 6.1 3.5 4.8 4.0 6.1 3.6 4.4 3.3 4.3 3.1 3.8 4.4 4.4 5.1 3.8
##  [307] 5.4 5.3 4.0 5.7 3.5 3.0 3.1 3.8 6.0 4.6 4.5 4.2 4.1 3.0 5.8 5.4 6.7 3.7
##  [325] 3.9 4.0 5.2 5.2 5.5 4.6 3.1 3.2 4.6 4.0 4.9 5.1 3.1 3.9 4.7 4.6 3.3 3.2
##  [343] 3.1 4.1 3.3 4.8 4.5 5.7 4.3 4.4 4.4 4.8 5.5 4.6 5.0 4.9 5.6 4.5 4.5 3.6
##  [361] 3.7 4.5 5.9 4.9 4.8 3.0 4.8 4.3 6.2 4.2 5.5 5.5 4.8 4.9 5.2 4.1 5.4 3.6
##  [379] 3.8 5.5 3.9 3.4 5.2 3.4 2.9 4.2 2.9 3.2 5.4 4.3 5.0 4.5 4.0 3.9 4.6 3.7
##  [397] 5.6 3.9 6.8 3.2 3.6 5.5 4.9 5.2 4.2 3.9 4.2 4.6 5.3 5.7 4.9 4.0 5.4 4.8
##  [415] 3.7 5.0 4.1 3.6 3.4 4.6 5.3 3.8 4.4 3.8 2.6 4.6 5.7 4.5 4.3 4.4 3.7 5.5
##  [433] 5.4 4.8 4.1 4.9 4.4 3.1 5.3 3.5 4.4 4.2 3.4 6.1 4.8 3.3 4.5 6.4 5.2 4.8
##  [451] 4.6 4.6 5.7 5.1 5.7 3.5 5.1 4.6 5.8 5.5 3.5 4.8 3.4 3.7 3.3 2.7 5.8 6.0
##  [469] 6.9 5.3 5.0 4.4 5.1 5.1 3.9 5.1 6.2 4.8 4.2 4.5 5.0 3.1 4.0 3.7 5.3 4.8
##  [487] 5.5 4.0 3.3 5.4 5.0 5.0 4.7 4.6 5.2 4.8 3.9 4.6 5.8 3.8 3.7 4.2 4.4 4.7
##  [505] 3.6 3.4 3.6 2.6 4.9 3.1 3.7 5.0 4.4 4.8 4.2 5.1 4.5 5.3 2.9 4.9 4.1 5.7
##  [523] 6.7 3.5 6.2 4.7 5.2 5.2 4.0 4.6 6.4 5.1 4.0 5.5 3.7 6.6 3.4 4.1 3.1 3.7
##  [541] 4.5 4.7 5.9 3.6 4.9 3.9 4.2 4.2 3.7 6.6 5.8 3.9 5.4 5.3 4.5 5.6 2.7 6.0
##  [559] 3.9 3.5 5.0 6.1 6.3 4.2 4.7 5.1 4.5 5.2 4.7 5.4 4.9 5.3 3.3 3.3 6.0 4.9
##  [577] 3.9 2.5 4.5 4.8 3.5 6.5 3.6 5.2 3.9 4.9 5.6 7.2 6.4 5.8 4.4 4.6 3.2 3.8
##  [595] 4.9 5.2 3.8 4.9 5.3 5.2 4.0 4.7 4.2 3.9 4.0 5.3 4.2 4.9 6.3 5.0 5.3 5.5
##  [613] 5.0 3.2 4.2 5.7 5.6 5.0 5.2 4.4 4.4 5.0 4.1 5.0 5.6 5.9 4.9 7.2 5.9 4.2
##  [631] 3.4 5.5 5.2 5.3 2.4 6.3 4.9 5.3 4.0 4.8 3.2 3.0 5.9 6.0 3.5 4.3 4.3 4.6
##  [649] 3.9 3.4 5.5 5.1 3.0 5.6 4.5 3.0 5.2 5.2 4.8 5.2 5.2 3.6 4.4 3.6 4.6 3.8
##  [667] 4.1 5.3 3.9 3.5 5.4 4.3 2.4 5.2 5.4 5.1 3.5 3.8 3.9 4.8 2.8 6.0 6.6 4.4
##  [685] 5.1 3.8 5.8 3.1 5.5 3.9 5.3 6.3 6.2 5.1 6.3 5.9 4.4 3.9 3.5 5.2 4.1 4.6
##  [703] 5.7 3.9 4.9 3.6 5.7 3.7 3.7 5.5 5.5 3.5 4.9 4.1 5.6 3.8 4.7 3.6 2.9 4.2
##  [721] 3.4 5.7 3.9 3.4 4.7 5.3 5.6 6.9 6.5 5.5 4.9 3.5 3.1 4.6 4.1 4.6 3.5 4.7
##  [739] 3.8 4.1 5.0 4.0 5.1 5.0 5.2 4.2 4.6 4.8 4.6 3.4 2.0 4.1 4.2 5.8 4.6 4.3
##  [757] 5.1 2.1 2.8 4.2 5.4 2.6 4.7 3.9 3.1 5.1 5.7 3.5 5.0 4.2 2.8 4.5 2.9 4.5
##  [775] 1.7 6.3 6.5 4.0 3.6 3.2 4.1 2.2 3.7 4.2 5.5 4.8 4.2 3.4 3.9 3.0 2.9 6.6
##  [793] 4.6 3.1 4.2 3.6 4.6 4.0 5.3 3.2 5.2 5.3 3.9 3.6 3.5 4.8 5.4 5.8 3.6 5.5
##  [811] 4.6 5.5 6.1 4.9 4.3 5.0 3.0 4.7 5.3 5.0 3.1 4.4 5.0 5.2 2.4 4.6 6.0 3.7
##  [829] 3.4 3.6 4.7 3.3 4.1 3.6 5.2 3.8 4.3 5.1 3.0 3.9 3.2 6.4 3.8 5.1 4.0 6.1
##  [847] 4.2 3.8 4.3 3.3 4.1 4.3 3.6 4.6 4.7 5.9 5.6 3.9 4.3 5.1 5.7 3.3 6.9 4.1
##  [865] 4.3 4.4 3.5 3.6 5.2 4.2 3.9 4.5 2.6 4.3 4.3 5.4 5.1 5.3 4.0 3.7 5.7 4.1
##  [883] 4.9 5.9 6.6 3.4 3.9 4.3 5.5 3.1 4.0 3.9 3.3 3.2 4.5 6.2 5.2 3.5 3.0 4.8
##  [901] 5.9 5.8 3.7 4.5 5.6 5.5 5.2 6.1 3.0 2.7 4.5 2.9 5.1 4.2 5.3 4.7 6.0 4.1
##  [919] 4.1 4.5 4.1 3.4 5.2 5.5 4.2 4.8 3.9 4.0 4.1 5.1 4.7 4.6 4.7 4.5 4.5 5.1
##  [937] 4.4 4.7 6.0 4.0 3.9 3.8 5.0 6.0 5.2 3.6 5.3 6.2 6.2 4.8 5.5 5.0 3.6 4.4
##  [955] 5.5 4.0 4.5 3.1 2.8 5.9 5.1 3.5 4.7 3.8 6.6 3.7 3.9 3.4 2.6 6.2 4.4 5.8
##  [973] 5.6 4.0 5.0 4.7 5.1 3.9 4.8 4.6 4.4 3.7 3.7 5.9 4.0 4.6 5.5 3.3 5.2 6.0
##  [991] 5.6 4.0 4.4 5.3 5.5 5.3 3.8 4.6 4.5 4.4
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
##   2.8   6.3
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
##    [1] 3.9 5.0 3.7 5.4 4.5 3.4 5.1 5.1 3.9 4.0 4.2 3.9 5.0 4.0 4.3 3.5 3.4 3.6
##   [19] 5.1 3.5 5.8 2.3 4.6 3.7 4.4 4.4 4.2 2.1 4.8 5.9 6.2 4.7 4.7 3.6 3.1 3.2
##   [37] 6.7 3.5 5.0 5.4 3.5 3.7 5.5 5.4 4.9 5.1 4.8 4.7 3.3 4.2 2.8 4.4 4.9 5.0
##   [55] 5.0 2.9 4.4 2.8 4.2 4.2 3.8 3.1 4.2 3.9 5.7 5.8 5.3 6.0 4.9 4.2 4.7 6.4
##   [73] 4.0 3.0 5.6 3.3 5.2 5.1 4.8 5.5 5.4 5.3 4.6 4.7 5.9 3.6 2.6 2.7 4.9 6.5
##   [91] 3.6 5.5 4.0 3.8 3.5 3.8 6.0 5.1 5.2 5.4 3.6 4.3 3.8 5.2 4.8 4.9 2.7 4.6
##  [109] 4.0 3.4 4.3 4.2 3.5 4.8 5.2 4.4 3.8 4.1 4.5 3.8 4.5 3.9 7.6 3.9 3.6 5.7
##  [127] 4.3 5.5 4.4 3.9 3.7 4.0 4.0 4.4 3.1 2.7 4.3 3.8 4.7 3.1 4.5 5.2 5.1 4.9
##  [145] 4.2 6.2 3.3 4.6 5.2 4.5 3.2 3.8 5.7 6.4 4.4 5.8 3.6 2.8 3.7 4.0 4.7 4.5
##  [163] 5.4 5.3 4.1 4.6 4.9 3.6 3.9 5.1 4.6 3.9 5.3 2.8 5.8 3.7 5.4 3.5 3.8 5.7
##  [181] 3.2 3.7 4.6 4.9 5.2 1.7 5.4 5.1 5.0 3.9 5.7 4.8 3.0 5.9 3.4 3.0 5.4 4.3
##  [199] 5.1 3.9 4.3 4.0 4.3 5.0 3.5 5.1 3.6 4.1 5.0 4.8 5.9 3.2 5.8 4.7 4.9 3.5
##  [217] 5.6 5.0 4.1 5.9 2.5 5.4 4.7 5.5 5.8 6.0 4.8 7.0 4.8 4.0 3.8 3.7 5.0 5.8
##  [235] 3.6 5.4 4.5 4.8 4.2 4.5 4.0 5.1 4.2 4.2 5.2 4.1 3.9 4.4 4.7 4.9 2.8 5.6
##  [253] 6.1 4.1 3.8 4.1 5.4 3.8 3.9 4.0 4.2 3.8 4.8 5.3 3.8 2.9 5.1 4.1 5.0 4.4
##  [271] 2.9 3.0 5.7 4.0 2.5 3.1 4.8 5.0 4.3 5.1 5.3 5.4 4.1 4.4 4.5 4.3 5.3 2.6
##  [289] 5.3 3.9 5.8 4.0 4.6 5.7 3.2 3.8 5.1 4.1 4.4 2.5 2.5 4.3 5.1 5.2 5.0 4.6
##  [307] 5.9 4.6 4.2 4.7 5.6 4.4 5.4 4.9 4.0 5.0 5.1 5.1 4.1 5.1 4.1 4.4 3.7 4.5
##  [325] 6.1 5.4 6.2 4.4 4.8 3.7 4.4 4.5 4.9 4.1 4.6 4.4 6.1 4.1 4.5 5.5 3.6 5.2
##  [343] 5.7 4.9 3.8 5.4 4.1 4.8 5.3 5.1 3.9 2.9 3.8 4.6 6.1 3.1 4.5 5.2 5.2 5.1
##  [361] 4.2 3.1 3.4 3.9 4.2 4.0 5.8 4.8 3.4 4.9 4.5 5.3 5.2 5.5 3.8 4.6 4.1 4.8
##  [379] 4.4 4.8 3.9 4.0 4.7 5.3 3.5 5.3 4.0 4.5 4.2 4.9 3.6 3.3 3.0 3.5 4.1 5.6
##  [397] 5.9 4.9 4.9 3.4 5.4 5.4 5.9 6.0 4.1 6.0 5.6 4.7 4.4 3.6 5.4 4.4 4.3 3.9
##  [415] 3.7 5.0 5.3 4.0 4.4 3.9 5.4 6.4 3.2 3.6 6.8 3.0 4.1 4.3 4.6 4.2 4.4 6.5
##  [433] 5.8 4.1 4.9 6.7 3.6 4.8 5.4 4.9 5.0 4.5 4.0 5.0 5.2 4.3 5.2 5.4 3.6 4.5
##  [451] 4.5 4.3 5.2 4.8 5.4 4.5 4.1 3.3 5.2 4.9 3.9 3.9 5.0 5.0 4.8 4.7 3.3 4.5
##  [469] 5.5 4.3 5.4 3.6 4.8 4.8 4.1 4.3 5.0 5.4 5.0 3.9 4.1 4.7 4.2 3.9 4.3 4.7
##  [487] 4.1 5.1 4.8 4.0 4.6 3.8 2.9 4.1 6.7 4.0 6.0 6.1 4.4 3.8 3.4 4.7 4.2 5.0
##  [505] 4.4 4.3 5.8 5.3 4.3 5.2 4.8 2.8 3.8 3.0 4.3 5.6 4.8 5.1 5.0 5.2 5.0 5.3
##  [523] 5.0 5.3 4.3 2.4 4.9 3.8 3.8 4.5 4.5 4.6 5.8 3.9 4.3 4.7 4.9 3.9 5.6 7.3
##  [541] 5.4 5.3 2.9 5.2 4.4 3.1 4.4 5.2 3.9 4.3 4.3 4.7 5.5 5.3 5.7 4.3 4.0 4.6
##  [559] 3.2 4.2 4.4 4.1 3.7 4.5 5.0 3.5 5.5 4.6 5.0 6.8 3.5 4.3 2.9 4.4 5.3 3.1
##  [577] 4.9 4.6 4.2 3.9 6.9 3.8 4.0 4.8 6.4 4.2 4.7 3.6 5.2 5.1 3.5 5.2 5.0 3.9
##  [595] 3.1 4.8 5.2 3.8 3.9 5.5 2.8 4.3 4.5 4.9 3.7 5.2 5.3 4.6 4.3 5.2 5.8 3.2
##  [613] 5.0 4.8 5.2 3.9 5.3 3.8 3.6 4.4 3.6 5.9 5.5 3.3 4.6 3.5 4.4 7.0 3.7 4.7
##  [631] 5.0 5.7 3.6 5.4 3.7 5.2 3.9 4.4 4.8 3.0 4.5 4.4 4.5 2.4 5.7 5.2 4.6 3.3
##  [649] 5.0 2.9 4.8 5.1 4.9 4.6 4.5 6.2 4.1 3.3 4.9 2.4 4.3 4.7 3.7 4.3 5.8 4.7
##  [667] 5.2 3.1 3.8 3.7 3.5 5.0 5.6 3.9 4.4 4.5 3.7 3.7 4.5 4.1 4.4 3.2 4.3 5.2
##  [685] 4.6 4.4 4.4 5.8 4.5 5.4 4.5 5.0 3.5 4.3 4.0 5.9 3.5 5.5 4.0 4.2 4.5 4.3
##  [703] 4.4 4.1 3.4 3.2 3.9 3.8 5.1 4.3 3.9 3.4 5.5 4.6 2.9 4.5 5.5 6.3 3.6 4.0
##  [721] 3.4 3.6 3.8 2.8 5.2 2.9 3.6 5.1 5.9 4.0 4.7 3.3 5.0 5.7 4.3 5.9 4.4 3.8
##  [739] 4.1 4.4 5.2 5.1 4.3 6.0 4.6 3.5 2.7 5.3 4.1 3.2 6.4 4.6 4.4 5.6 3.8 3.7
##  [757] 3.3 4.4 4.5 4.4 3.6 5.0 5.4 3.8 5.0 2.9 3.7 4.8 4.8 4.1 2.8 2.3 3.8 5.4
##  [775] 6.6 3.6 3.8 4.0 5.2 3.7 2.9 3.5 5.2 4.1 4.9 6.5 4.8 4.6 4.5 5.5 4.4 4.3
##  [793] 3.9 4.5 4.1 4.6 5.3 4.3 6.3 4.7 4.5 2.5 3.7 3.9 4.5 5.0 4.8 2.7 5.4 3.3
##  [811] 4.6 5.7 5.3 4.8 4.7 4.9 4.3 3.2 4.5 5.4 4.5 4.3 4.6 3.9 5.8 3.8 3.6 6.1
##  [829] 4.2 5.0 2.7 6.1 7.2 5.4 6.0 5.6 4.0 5.6 4.7 4.9 4.5 3.9 4.8 4.4 3.6 4.5
##  [847] 2.4 4.8 5.0 5.2 4.3 4.9 4.0 3.5 3.8 4.9 4.7 4.3 5.8 5.5 5.5 3.1 5.1 4.5
##  [865] 5.8 5.6 5.8 2.5 3.0 5.0 3.7 4.9 4.3 3.1 4.5 5.1 4.7 5.4 2.9 5.5 5.8 4.0
##  [883] 6.2 3.8 6.7 4.4 5.3 4.0 3.1 4.4 5.4 4.3 3.9 4.2 3.0 4.6 3.7 5.0 3.7 5.5
##  [901] 4.3 5.2 5.0 3.9 3.9 4.2 4.4 6.0 4.3 7.4 5.1 6.7 4.6 6.5 3.8 3.3 4.2 3.1
##  [919] 4.9 4.7 5.4 4.3 4.2 4.5 4.4 4.9 3.7 4.1 4.8 2.7 4.8 5.5 4.4 3.0 5.8 6.2
##  [937] 5.7 4.4 4.8 5.8 5.0 3.6 5.2 5.7 4.4 6.3 4.7 3.6 5.4 5.5 5.1 4.5 3.4 6.3
##  [955] 4.7 4.8 4.4 5.0 4.1 5.2 3.4 3.9 3.1 4.9 5.0 5.2 5.9 4.8 4.0 5.3 4.8 6.0
##  [973] 3.0 3.3 4.8 2.3 5.2 5.3 5.3 4.9 5.0 4.5 3.3 5.4 4.9 5.2 4.3 4.9 5.7 5.2
##  [991] 3.2 5.4 3.1 5.6 4.8 4.2 4.0 5.0 6.0 4.9
## 
## $func.thetastar
## [1] 0.0187
## 
## $jack.boot.val
##  [1]  0.48850575  0.40946746  0.31077348  0.21416431  0.20804598 -0.06772334
##  [7] -0.15058480 -0.20083333 -0.41260745 -0.44424779
## 
## $jack.boot.se
## [1] 0.9555971
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
##    [1] 4.5 3.5 4.9 2.9 4.8 5.5 4.5 3.9 2.6 5.7 5.0 5.0 3.0 4.4 5.0 4.5 4.6 2.8
##   [19] 4.0 2.8 5.4 3.6 4.6 4.0 5.1 4.8 3.0 4.5 4.1 4.4 3.6 4.1 4.7 3.1 6.1 3.5
##   [37] 3.8 5.3 5.1 5.2 3.7 4.9 5.1 5.5 5.0 5.8 5.0 5.0 4.1 5.3 2.5 3.8 3.3 6.2
##   [55] 5.6 5.7 3.4 4.6 4.1 4.3 3.7 4.2 5.5 4.7 3.8 4.1 4.0 3.8 4.4 6.6 4.9 5.1
##   [73] 3.3 3.8 3.7 5.8 3.2 2.4 3.9 3.9 5.2 4.7 4.0 5.5 2.9 3.9 6.0 4.5 5.4 4.9
##   [91] 5.7 4.2 4.4 4.4 3.4 6.0 4.8 3.4 6.4 5.2 4.7 4.5 5.0 3.8 5.9 3.4 3.6 5.4
##  [109] 7.0 4.0 4.8 3.6 5.5 5.0 5.3 6.3 3.6 3.8 4.1 5.8 5.6 4.2 5.3 4.6 4.1 4.5
##  [127] 5.1 3.9 5.4 6.1 5.4 5.4 5.2 5.1 4.7 4.6 5.7 4.1 6.0 5.5 5.1 4.3 5.0 5.2
##  [145] 4.3 4.0 4.1 3.1 3.9 4.6 3.9 3.6 6.2 4.0 2.9 4.9 5.5 4.0 5.3 3.4 3.4 5.6
##  [163] 4.7 3.6 2.9 5.4 4.4 4.1 4.4 3.9 5.1 3.6 5.4 5.2 4.8 4.6 5.1 4.5 5.4 3.9
##  [181] 5.3 4.5 3.6 5.3 4.3 3.6 4.4 4.8 4.4 3.1 4.0 3.7 4.9 5.2 4.9 4.1 5.2 3.8
##  [199] 4.2 5.8 3.4 4.6 3.4 5.0 3.5 3.3 5.2 4.7 3.9 3.1 5.0 3.2 6.2 5.2 5.2 2.9
##  [217] 3.8 3.7 5.0 5.9 3.3 5.8 4.6 4.8 4.2 4.0 3.8 2.5 5.9 4.4 3.7 4.7 5.3 5.1
##  [235] 4.5 5.0 6.3 3.8 3.8 4.6 2.7 4.4 2.6 4.6 4.7 5.1 3.6 5.0 4.9 4.3 4.7 4.5
##  [253] 4.5 4.0 3.5 4.2 2.0 4.4 4.9 4.0 6.3 4.3 4.3 4.8 5.4 3.0 4.7 2.8 4.4 3.5
##  [271] 5.7 5.7 4.6 6.9 3.3 4.3 4.4 5.3 4.8 3.8 5.1 2.7 5.4 3.8 4.7 4.5 3.6 4.0
##  [289] 5.8 5.4 4.7 4.4 5.8 2.9 4.6 4.4 3.2 5.3 6.0 5.6 4.0 6.1 4.2 4.8 5.0 5.4
##  [307] 4.0 3.2 4.1 3.8 4.7 4.1 5.3 5.1 4.5 2.9 2.8 3.8 4.0 5.5 4.0 5.1 3.2 4.1
##  [325] 4.7 3.1 3.6 5.8 4.2 3.4 4.4 4.9 5.0 4.8 4.7 3.5 3.2 4.0 4.9 4.4 4.4 4.8
##  [343] 6.2 4.1 4.0 5.0 4.7 5.1 4.0 4.3 4.2 3.3 4.4 5.2 4.9 5.2 5.4 4.2 4.1 4.4
##  [361] 5.5 5.0 5.3 4.1 3.6 5.0 6.4 3.8 5.7 5.3 4.4 3.4 5.6 5.2 4.4 4.7 4.1 5.6
##  [379] 5.1 4.2 4.4 7.1 3.4 4.1 4.6 3.4 4.6 4.6 4.3 3.8 4.0 5.3 4.2 3.5 5.1 5.0
##  [397] 3.2 5.5 3.8 3.0 4.7 5.2 4.1 2.6 5.2 2.8 3.5 4.9 4.4 5.4 4.1 4.9 4.0 5.6
##  [415] 4.7 4.6 5.3 4.6 3.7 4.2 3.8 5.0 4.4 5.6 3.8 5.3 3.1 3.8 5.0 4.2 3.7 5.7
##  [433] 4.9 4.9 4.8 4.5 5.1 5.1 3.3 4.4 6.3 5.7 3.2 4.3 5.3 4.0 4.5 4.7 5.3 5.4
##  [451] 4.4 5.5 5.0 4.5 4.8 3.6 6.3 3.6 5.4 4.9 4.0 4.0 4.8 4.4 3.7 3.0 4.4 4.1
##  [469] 4.4 5.8 4.8 4.6 5.1 3.2 5.2 3.9 3.2 4.6 2.6 2.5 3.7 4.7 5.8 4.9 3.8 4.4
##  [487] 5.8 4.7 3.1 4.0 4.0 4.1 4.5 5.0 4.5 5.5 3.7 4.7 4.7 5.3 4.2 3.5 3.4 5.6
##  [505] 4.3 6.1 4.8 6.4 3.6 4.3 5.5 4.3 3.9 4.2 4.0 3.6 3.5 5.0 4.7 3.6 4.7 4.2
##  [523] 3.9 4.8 4.2 2.3 2.9 3.8 4.6 4.0 5.9 5.2 5.4 4.2 4.3 3.6 4.8 3.9 4.5 5.7
##  [541] 4.8 5.2 3.9 4.1 4.8 3.6 4.8 3.4 3.4 4.7 3.3 4.2 5.8 5.8 4.3 6.0 5.7 4.0
##  [559] 5.2 3.1 5.0 5.8 5.4 4.8 5.1 5.6 5.6 5.1 3.5 4.3 2.9 4.2 6.4 5.5 5.8 3.3
##  [577] 4.3 7.1 4.0 3.8 5.1 3.5 5.3 5.4 6.2 3.7 4.6 4.3 4.5 5.8 4.1 4.0 3.3 3.0
##  [595] 3.9 3.9 4.7 5.2 3.4 5.1 3.4 5.1 3.7 4.0 3.8 3.4 4.0 4.9 3.1 3.3 3.4 4.4
##  [613] 4.8 3.9 4.7 3.9 3.9 4.4 4.3 3.5 3.9 4.6 3.3 3.5 6.1 4.7 2.9 4.0 4.4 3.7
##  [631] 2.8 4.6 5.1 4.5 4.2 4.3 4.4 4.2 3.6 4.8 5.5 5.0 5.8 5.5 5.5 5.4 6.1 4.6
##  [649] 3.4 3.1 3.6 4.0 4.6 5.6 3.7 5.7 3.8 4.7 4.7 3.8 5.6 3.3 4.0 5.1 5.0 5.1
##  [667] 3.7 5.0 4.8 4.5 2.6 5.0 6.3 5.2 6.3 3.6 4.1 6.3 2.8 4.2 3.0 5.8 4.7 5.1
##  [685] 5.0 4.1 4.9 4.9 4.3 4.7 3.9 4.3 5.8 2.6 3.5 4.8 5.0 4.3 4.4 5.9 4.0 5.4
##  [703] 5.7 4.6 3.4 3.6 3.5 3.8 5.4 4.2 5.7 3.8 3.7 5.1 5.0 3.7 5.1 4.9 6.1 5.7
##  [721] 3.8 4.3 4.5 5.7 4.0 4.4 6.1 4.7 4.0 4.1 2.5 5.1 4.4 4.5 4.6 5.8 5.0 3.4
##  [739] 5.2 4.5 6.0 5.1 3.4 5.6 4.3 4.1 4.4 4.7 4.1 3.3 2.8 4.0 4.2 4.9 4.7 3.8
##  [757] 5.4 3.4 5.5 4.5 5.7 5.0 3.9 5.7 4.7 5.4 5.1 4.4 4.3 5.2 5.3 4.2 2.8 3.2
##  [775] 5.2 3.2 2.1 4.6 3.8 3.3 3.6 5.1 3.6 5.4 4.5 4.0 4.4 4.9 4.9 5.3 4.7 4.5
##  [793] 4.7 3.9 5.1 4.0 5.4 5.0 5.1 4.7 5.9 6.2 4.8 3.8 5.2 3.7 3.6 3.3 4.1 2.6
##  [811] 4.6 4.8 4.9 4.8 5.3 4.2 4.8 4.8 5.1 5.6 5.6 4.3 4.3 3.5 4.0 3.8 4.0 5.9
##  [829] 4.1 4.9 4.3 4.3 3.4 5.2 3.8 4.7 4.1 4.4 4.7 3.0 3.0 3.1 3.7 4.4 5.1 5.4
##  [847] 5.9 3.6 3.9 5.6 3.7 3.8 5.5 5.5 2.9 4.4 2.4 3.4 4.9 2.8 4.9 3.7 5.0 4.1
##  [865] 5.8 4.0 5.3 3.4 4.5 4.8 4.4 4.5 3.6 4.8 4.6 3.1 3.2 4.8 5.0 7.6 3.6 6.6
##  [883] 4.3 6.9 3.3 5.1 4.1 4.5 3.6 4.2 4.0 6.0 5.6 4.6 5.1 4.3 5.4 4.5 2.9 4.9
##  [901] 4.4 4.9 4.6 5.0 4.1 3.8 5.0 4.4 6.9 5.4 2.8 5.4 5.7 5.0 4.9 4.1 4.2 5.0
##  [919] 6.5 3.0 3.8 5.1 2.8 5.4 3.4 5.9 4.8 4.5 6.3 5.0 4.5 5.0 5.4 3.4 4.4 4.8
##  [937] 5.0 5.5 3.3 5.7 6.1 4.3 5.1 5.4 3.2 4.7 2.9 3.4 4.2 3.5 3.8 5.2 4.5 4.3
##  [955] 4.9 4.2 4.7 4.4 3.9 3.8 4.5 4.5 6.0 4.5 5.6 4.0 4.9 6.4 5.5 3.6 3.5 3.1
##  [973] 3.6 4.3 5.2 4.8 4.9 4.3 5.2 6.1 5.8 3.4 4.8 4.2 3.4 4.0 5.5 4.6 5.4 4.6
##  [991] 4.3 4.0 3.7 5.1 6.3 5.2 3.7 3.7 6.0 4.9
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.500 5.364 5.200 5.200 5.100 5.000 4.948 4.800 4.500 4.400
## 
## $jack.boot.se
## [1] 1.006314
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
## [1] 1.364109
```

What is skew? Skew describes how assymetric a distribution is. A distribution with a positive skew is a distribution that is "slumped over" to the right, with a right tail that is longer than the left tail. Alternatively, a distribution with negative skew has a longer left tail. Here we are just using it for illustration, as a property of a distribution that you may want to estimate using your data.

Lets use 'fitdistr' to fit a gamma distribution to these data. This function is an extremely handy function that takes in your data, the name of the distribution you are fitting, and some starting values (for the estimation optimizer under the hood), and it will return the parameter values (and their standard errors). We will learn in a couple weeks how R is doing this, but for now we will just use it out of the box. (Because we generated the data, we happen to know that the data are gamma distributed. In general we wouldn't know that, and we will see in a second that our assumption about the shape of the data really does make a difference.)


```r
library(MASS)
fit<-fitdistr(original.data,dgamma,list(shape=1,rate=1))
```

```
## Warning in densfun(x, parm[1], parm[2], ...): NaNs produced
```

```r
# fit<-fitdistr(original.data,"gamma")
# The second version would also work.
fit
```

```
##     shape       rate  
##   4.896109   9.686837 
##  (2.119220) (4.415464)
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
## [1]  0.1307933  0.7464766  0.7736431 -0.2782426  0.4023117  0.8073835
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
##    [1]  1.1050565331  1.1919406630  0.0715199808  0.5333275942  1.3781649499
##    [6]  2.0763126616  0.5829089657  1.3892239184  0.7265251180  1.2755686703
##   [11]  0.1803270159  2.0606530875  0.8032846205  1.1992533373  1.9694295609
##   [16]  1.8199485434  1.4346508478  1.0347233655  0.6349868695  0.1219058596
##   [21]  0.7560889790  0.9764402618  1.6290573753  0.8314551766  0.1588231310
##   [26]  1.1460248646  2.1967880361  1.6135259778  1.8267086025  0.4685807712
##   [31]  1.1257326805  1.3433788798  0.8985850829  0.6172609831  1.2561660567
##   [36]  1.8740550376  0.7852060279  1.6549479848  0.4122114601  0.1345330787
##   [41]  2.1811594032  0.5192579060  1.4368886855  0.8595689096  0.5364249590
##   [46]  1.1452873669  0.2287405661  0.1051882493 -0.2756757439  1.6094325678
##   [51]  0.5394572809  0.0130101656  1.0234124282  1.3233265512  1.4201749154
##   [56]  0.9496019539  1.0396311603  0.9007307825  1.3735667896  1.8842639367
##   [61]  0.5627619268  0.9247321158  0.5815425904  0.8906225998  1.5332415386
##   [66]  2.2701739264  1.4380297050  1.1374323809  1.6986567176  2.3521568552
##   [71]  1.9850163926  0.4738253243  2.1101625479  0.1453813218  0.3663590514
##   [76]  0.4114852282  0.2038036724  0.8481131813  1.1890313382  0.7870531812
##   [81]  0.1258400495  0.8319628368  1.8357167965  1.8137381039  2.3072610654
##   [86]  1.1839489249  1.2846090981  0.4755857174  0.2931742233  0.5075557020
##   [91]  1.3207927193  0.8034295414  1.3769498491  1.5857490266  0.4566378656
##   [96]  1.9128297794  2.0395320484  0.1318004537  1.3570075858  1.3287639014
##  [101]  0.8206750429  1.4940325172  0.4309876027 -0.2664062563  0.9870124179
##  [106]  2.2295810658 -0.1042173649  1.2695559279  0.6265606179  1.3633886800
##  [111]  1.2725657624  0.9449981327  1.0416865642  0.8069030997  2.2148120139
##  [116]  0.8566200057  0.3876490726  1.0478892519  0.5041019982  0.7017242086
##  [121]  0.9462667911  0.9718737004  0.7518397370  0.4278326900  1.7239439756
##  [126]  0.7203561026  0.6161104318  0.9079012210  1.1775716479  0.9101788257
##  [131]  1.4447908964  1.4745922722  1.3089017260  0.4236889070  0.7445630255
##  [136]  0.9563542099  0.7671770920  0.8497489923  0.3221896180 -0.1802864337
##  [141]  0.7060059667  0.8483117590  1.8537951438  0.5929681714  0.5280364329
##  [146]  0.9918163663 -0.6305375550  0.4000930252  1.5255987219  1.3639389661
##  [151]  1.4977649898  1.6027934659  0.7535179153  0.9368715337  1.4756321600
##  [156]  1.4454514349  1.7489167527  1.8206570432  0.4367807053  0.7617417706
##  [161]  0.8550625110  1.3433743633  1.2559814566  1.5060872659  1.0040830080
##  [166]  1.2306171075  1.2094781752  1.5133430375  1.6430826934 -0.2121781035
##  [171]  1.3683956582  1.0313505012  2.3238965822  2.2540845345  0.3433678685
##  [176]  1.4514101215  0.0510683103  1.2517691512  0.4203925071  0.8978794678
##  [181]  1.5337491668  1.3017079869  0.9245793677  0.4244346573  0.8138616459
##  [186]  0.1338859870  1.0884527220 -0.0055538023  0.7000514058  1.5619685085
##  [191]  0.8599074330  1.5728678602  1.4250882883  0.3981701244 -0.2078794105
##  [196]  2.4568437651  0.7587364602  1.2880269345  1.2849545908  0.0182065293
##  [201]  0.6993077195  0.5026916682  0.7680562795  0.9463986995  1.0923377913
##  [206]  0.3886855215  0.6752243629  0.7320734155  0.3883736599  0.8705087667
##  [211]  0.4312461993  0.9841988300  1.3869906934  0.7404535083  0.8460927038
##  [216]  0.8012270139  1.3978653820  0.9059382122  2.0666713870  0.5984800403
##  [221]  1.2793777033  2.0009915742  0.8370059728  1.1258470946  0.9125978742
##  [226]  1.0076335901  0.8441774086  1.4607044557  1.2268262573  1.4857322805
##  [231]  1.2375347827  1.0408858674 -0.8883503278  0.9484164372  0.3392999954
##  [236]  1.0222635747  0.8283073828  1.3231874784  1.4584468695  0.8273707793
##  [241] -0.3193548036  1.6525253146  1.1371693585  0.7023378240  1.0933334432
##  [246]  0.6026660382 -0.5943231902  0.4474576992  1.7589899897  0.6273869733
##  [251]  1.5014272285  0.2879442569  0.9029921530  0.5977149805  0.8627053328
##  [256]  0.5291285165  0.8103077707  1.8831742528  1.6896105519  0.4953561700
##  [261]  1.4607044557  0.4019527973  0.9566399213  1.7246478233  1.0057388371
##  [266]  0.8742026006  0.8547060618  1.4328438624  1.1170702810  0.9954322561
##  [271]  0.7435786260  0.0230065756  0.8631216401  1.5798455816  1.5325932267
##  [276] -0.0128012184  1.2732293969  1.0164229799  0.4110982987  2.2943721212
##  [281]  1.2492972026  1.6334396047  0.4420681272  1.6690740186  1.4421161619
##  [286]  2.1210565923  1.1890664535  0.3387547556  1.1464221616  1.0497411004
##  [291]  1.1135087862  1.2442807463  0.7233579078  0.4529833448  0.8004138555
##  [296]  1.0284278534  1.3456799740  0.1740951940  0.8049204859  0.5020315931
##  [301]  1.2429692775  2.5451455881  0.0869994096  0.7101975948  0.0901353891
##  [306]  1.0720597524  1.4070117432  1.3126185383  1.3449448482  0.4854654005
##  [311]  2.1127759804  2.4070415183  2.0842701931  1.5055148374  1.6617155603
##  [316]  1.2941291931  1.3523917074  1.5296945536  2.3521564479  1.7632548010
##  [321]  0.9256098162  0.3964299354  1.5607924794  0.2292452808  2.1666296459
##  [326]  0.7702666869  1.4834977105  1.3255521347 -0.2410542065  1.4331700174
##  [331]  1.1981298581  1.8721817176  0.4752364085  0.4834937215  1.0679726855
##  [336]  1.1214642728  0.8422285052  2.0752812749  2.2243748012  1.0148033336
##  [341]  1.0557308534  1.1524126680  0.9285667709  0.8339057925  0.9294891855
##  [346]  0.8560872201  0.8690853933  0.9695868042  1.9074469910  0.4963743707
##  [351]  0.9300935887  1.5602967834  0.4624766913  0.7690079327  0.6834882106
##  [356]  1.0425166084  1.3217856501  0.0277188701  0.9407507051  0.3537809371
##  [361]  0.8741709906  0.5493613687  1.3284743340  0.9614619411  0.8637191773
##  [366]  1.3595990133 -0.1877319569  2.5990689585  0.5685947591  1.5868478988
##  [371]  0.9592124515  1.5497551427  1.2810638678  1.2878060780  0.7309475000
##  [376]  1.2509725066  2.1107695174  2.5435170359 -0.0002820642  0.9191545107
##  [381]  0.4476505593  1.2475728092  0.7170438191  0.2752582876  1.2995155273
##  [386]  0.7246686199  0.4897689325  1.3369870274  1.3369870274  0.8634487776
##  [391]  0.4912464180  1.2956009698  1.2338664738  1.2812256532  0.7267076775
##  [396]  1.1106730560  1.8264731223 -0.6578385798  1.6602100486  0.4252316246
##  [401]  0.7644501292  0.9773556515  0.8939657885  0.9723337680  1.0508172105
##  [406]  1.8511888749  0.9743613641  1.3697129005  0.7856586310  1.5034598537
##  [411]  1.2788277731  1.8660446747  0.8264321241  0.2169619455  0.7621838792
##  [416]  0.4087870361  0.6533975952  1.8304087468  0.4045370583 -0.6068449340
##  [421]  1.9087702101  1.7018109751  1.0508624960  0.9808755736  1.4452101841
##  [426]  1.4003022623  0.0900913813  1.5248120799  1.0257868721  2.0938988031
##  [431]  1.6028341997  1.0384620291  0.8407130278  1.0376602256  0.8947221531
##  [436]  0.4118136055  1.9621245036  1.2800007180  0.7671860207  0.9045711072
##  [441]  1.4181842854  1.1372057045  1.1500674920  1.0339971864  0.7620288157
##  [446]  0.5393292468  0.9412562136  0.5825256474  1.3411175031  1.5582627824
##  [451]  1.2687029073  0.3732744900  1.2531812909  2.1424660057  0.8358440801
##  [456]  0.6855952445  1.3474964144  1.3116357837  0.5972158173  0.8526317014
##  [461]  0.5614712897  2.0421809550  0.7094486912  0.6130115692  1.7233254454
##  [466]  1.7188323618  0.5099797893  1.1079368314  1.1627230562  1.4687322174
##  [471]  1.1037806369  1.9030637169  1.8124707622  0.7371895855  0.5918976656
##  [476]  2.6353222629 -0.2719821956  1.1450180013  0.9704389687  1.4470668293
##  [481]  1.0625428413 -0.2412904915  0.9613904557  1.2851436644  1.0969143937
##  [486]  1.9130211315  0.4756511459  1.1491995226  1.1151973634 -0.6198859513
##  [491]  1.3678071828  0.0761811359  1.6697435887  0.7717371496  1.5549564033
##  [496]  0.9755364981  1.4071391874  0.6010285452  1.6291764814  0.7624991419
##  [501]  0.9157758388  2.1235760094  0.1497992377  1.6538014365  0.2474114410
##  [506]  0.3844621688  0.7961197615  0.9596384814 -0.2427230953  0.8622584254
##  [511]  0.1642304928  2.1177024026  1.2575311770  1.2357156285  0.3998896273
##  [516] -0.2978118699  0.8711847752  1.3683412914  0.7803784705  1.4535462547
##  [521]  1.1257326805  0.9357519981  1.8484030946  1.2073073792  1.8586182113
##  [526]  0.4069111906  0.1323323897  1.2900111384  0.8985850829 -0.2531518979
##  [531]  1.4862635629  2.0123876734  0.6895554947  0.7088749804  1.3936906509
##  [536]  1.2928506018  0.8664417884  1.1180517905  1.8908773940  0.9395612418
##  [541]  1.5258723598  0.7085565383  1.4839582645  1.0940256903 -0.1106154173
##  [546]  1.5919044956  1.4239115067  0.4806549574  1.2638306695  0.3902369863
##  [551]  0.5284618555  0.0320831266  1.1296636623  1.2807401606  1.1970904862
##  [556]  1.0917614565  0.8395152807  1.1799134266  1.9571656671  1.1036039249
##  [561]  0.4145869608  2.0692457809  1.2123213468  0.7995096063  1.8154478843
##  [566]  0.9287642592  1.2681167139  1.1353681067  1.6290573753  1.3084075728
##  [571] -0.0496926087  2.3207891033  0.3151291348  0.3602159882  1.5244338454
##  [576]  0.4595834056  1.1368461369  0.5266983470  1.4437871358  0.8609206093
##  [581]  1.6913106307  0.3952967028  0.9519702114  0.5972158173  2.1977686087
##  [586]  0.7019146760  0.7779920078  1.1571421293  1.4188904457  1.0632898016
##  [591]  1.6397851984  0.1377172570  1.2499452649  0.8052590948  0.0208939526
##  [596] -0.2598760430  0.4281780990  0.8733931972  1.4710106375  0.9602293233
##  [601]  2.2298286931  0.8051963143  0.8027004530  1.1535006902  0.8552446302
##  [606]  0.7467571161  1.2768212271  0.6944046942  1.1548342884  1.2288215611
##  [611]  0.4542474388  0.3783607497  1.5656012886  1.9936741429  0.4993827602
##  [616]  2.5787642939  0.3882034890  0.8351010715  0.8913214073  1.2090676162
##  [621]  0.6076381248  1.1470943707  1.4869134505  1.6690740186 -0.1831375085
##  [626]  0.7135042924  1.3257422065  0.8665888863  0.3192970245  1.0951111991
##  [631]  1.1372057045  0.8307033392  0.7239133632  1.4247928245  0.0150403843
##  [636]  0.7440881074  1.1358152293 -0.4070263174  2.6320618005 -0.2722309648
##  [641]  1.3641094684  1.3391808618  0.7269945434  2.2181649111  0.8719931547
##  [646]  1.1723773597 -0.1501015330  1.4869134505  1.0200160482  1.1233220149
##  [651]  0.1652134769  1.3629254092  1.0166670721  0.9298990696  0.1542084468
##  [656] -0.0069544577  0.6711879501  0.8022108780  0.8639512597  1.4614323268
##  [661]  0.7428253488  0.4925974711  0.8315265581  1.3774207194  0.8075452557
##  [666]  0.7444523481  1.8125441445  1.8705579515  1.5781155843  1.1112084981
##  [671]  1.2648917387 -0.0002323138  0.7668268704  2.2455281429  0.8125866233
##  [676]  1.4940207508  0.9474358549  0.9241671660  1.2805015086  1.0661724463
##  [681]  1.1017156485  1.4180219721  0.6145769062  1.8661613452  1.1691630073
##  [686]  2.4435331029  1.2905772499  1.0745626574  0.6526211496  2.0210134572
##  [691]  1.3397015491  0.4804646817  0.1327086560  0.6047573532  0.9710191310
##  [696]  0.8606788012  1.2938273778  1.2475728092  1.0875964185  1.0219016108
##  [701]  2.1712180444  0.8606788012  0.9539719128  2.1143317582  1.0189774675
##  [706]  2.2920117991  0.5884075794  0.5918696330  1.4569573140  0.9658011679
##  [711]  0.4840533794  1.1949159180  2.2067814799  0.1781387310  0.7588113324
##  [716]  1.0082013878  0.7449262832  0.4795258228 -0.0638759366  0.4231658926
##  [721]  0.5195939860  0.8439518545  0.8495447562  0.3940156647  1.0960625994
##  [726]  0.4410886491  0.8861443479  1.3339903627  1.8327015715  1.1196348710
##  [731]  0.8638264869  1.2083379604  0.8900331661  1.6774231289  0.7132647217
##  [736]  0.2701885376  0.6118476593  2.2461054099  0.6263239743  0.3687312074
##  [741]  1.0461964077  1.4639501323  0.7356968859  0.5146710899  0.9870634212
##  [746]  2.1939358100  1.5316982737  0.7156890223  0.4848723389  1.1341727499
##  [751]  0.3859811128  0.5630893467  0.2688485096  0.9252202138  1.0720993259
##  [756]  1.2968264554  1.8167818695  0.6893868960  1.5740883452  1.0931047622
##  [761]  0.8889033230  1.7097238779  0.9374762875  1.1310902072  1.5060854518
##  [766]  2.4450004862 -0.8235173829  0.8885218684  0.5725362734  1.3469837250
##  [771]  0.9522594167  1.4436194394  1.2466617916  1.5656012886  1.2389965916
##  [776]  2.3121165118  0.4491078670  2.1154231019  0.2169619455  0.6969796609
##  [781]  1.6288717431  2.0771131501  1.5657080148 -0.1873741990 -0.1543300204
##  [786]  1.3740717517  0.7165019944  1.0854254859  0.8186665677 -0.0099072663
##  [791]  0.7469013048  0.9311293673  0.4675903936  1.8484342057  0.1436280865
##  [796]  1.2307893866  0.1122022420  0.6789113803  0.8306407414  1.2961370000
##  [801]  0.7165019944  0.1815456033  0.9758881833  1.4719925357  0.9898180587
##  [806]  0.9199755512  1.3029730782  1.1560623541  1.4802881586  1.4452083907
##  [811]  0.5007364272  0.8860890092  1.8120638019  1.0381565133  1.8859598859
##  [816]  0.4663832907  1.3250084999  0.8713200779  2.6163220014  1.4679056320
##  [821]  0.4466184506 -0.6457074613  1.2334394886  0.4725388715  0.9206058373
##  [826]  0.2017443037  0.1995417650  1.9574412770  1.1920728328  0.9931624353
##  [831]  0.7478054095  1.6270883991  1.4118111094  1.3376467457  0.7095249943
##  [836]  1.8419055259  0.4052723683  0.6483466772  0.2786533306  0.5048160596
##  [841]  1.5582627824  0.8832972890  1.4234428564  1.2629364965  0.6986062154
##  [846]  1.4236687748  1.1266878665  0.6089861529  1.6037543186  0.9574123817
##  [851]  1.1151299957  1.2083782383  1.2256690271  1.0400651871  1.1502255009
##  [856]  0.8155639850  0.8293426274  0.5070059747  0.8417985850  1.6276962374
##  [861]  0.7843047562  0.4066445094  0.9218700318  1.4643458584  0.8905265534
##  [866]  0.7018155243  0.8143035858  1.1560623541  0.9795308955  0.8589295660
##  [871]  0.4545460344  0.3740289360  0.2147744282  0.7385986882  0.4404526385
##  [876]  0.8128681797  0.8242057903  0.6205908612  0.9195073522  0.0094205835
##  [881]  1.1315024830  1.5707585650  0.8158956495  1.1663682775  2.3232799405
##  [886]  1.5746295888  1.0616829524  0.0774368352  0.5108285188  0.9174137204
##  [891]  1.6791918371  0.1607695185  1.5488319671  1.2517900611  2.5280343096
##  [896]  0.8675410900  1.0862090442  1.2010721171  1.2452252776  0.9099004618
##  [901]  0.9273505032  1.4656846361  1.1839483593  1.1964356973  0.0150403843
##  [906]  1.1410597548  0.9171200915  0.9986188859  0.8520230559  0.3889022325
##  [911]  0.1558404116  1.1149463208  0.3981701244  1.2983698081  0.2653763005
##  [916]  1.6580769310  0.7974809297  0.5758205554 -0.2763716477  0.3794034592
##  [921]  0.4055325355  1.8666236905  0.4557828803  1.2888581654  1.3123464368
##  [926]  0.8520499913  2.4626397930  1.7931153179  0.5286283599  1.2835319205
##  [931]  0.9904208433  0.4066445094  1.4263636490  1.3355471402  1.7347639871
##  [936]  1.5113512531  1.1228599873  0.9083055349  0.9641845561  1.3015366622
##  [941]  1.0037035669  1.0258184738  0.4069649240  2.5634491607  1.1619846538
##  [946]  0.8247106226  0.0304841001  1.4114228308  0.8515102248  1.2823283990
##  [951]  0.9956495040  1.8065183718  1.4862635629  0.7767530879  1.5860528139
##  [956]  1.2082547634  1.4499872107  0.3859121544  0.9433517514  1.5065522939
##  [961]  0.9218287572  0.1799701816  1.3848881239  1.9906581680  0.4721815203
##  [966]  0.5668713062  0.4351139962  1.3109540191  2.0080181881  0.5546988319
##  [971]  0.6568218925  1.1479121850  0.8583699163  1.1334216787  1.0595405594
##  [976]  0.8961064888  1.7331307989  0.8326054305  2.3324540485  1.2755475223
##  [981]  0.9466894082  0.2964834679  1.0298063397  0.8524930844  1.3575195086
##  [986]  0.9057562477  2.2074069760  0.2247266871  1.1081058568  1.9864520041
##  [991]  0.3673099665  0.8238124290  0.4695133729 -0.2708085504 -0.2961781691
##  [996]  0.8441443316  0.4920923002  2.1213758295  0.5816608565  0.3125537617
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
```

```r
fit2
```

```
##       mean          sd    
##   0.50544198   0.25672484 
##  (0.08118352) (0.05740228)
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
## [1]  0.46835475 -0.81463707 -0.06509541  0.10837276  0.73572925  0.64236828
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
## [1] 0.0302
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9071112
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
## t1*      4.5 0.05645646   0.9179658
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 1 2 3 5 6 7 
## 2 3 2 1 1 1
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
## [1] 0.044
```

```r
se.boot
```

```
## [1] 0.9276554
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

