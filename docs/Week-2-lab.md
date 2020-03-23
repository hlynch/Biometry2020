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
## 0 2 3 4 5 8 9 
## 1 1 1 2 1 1 3
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
## [1] -0.0054
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
## [1] 2.766831
```

```r
UL.boot
```

```
## [1] 6.222369
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##   2.5%  97.5% 
## 2.7975 6.2000
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
##    [1] 3.9 4.0 4.4 5.6 3.3 4.8 3.8 4.6 5.4 4.5 3.2 4.3 4.4 4.0 3.7 4.6 6.2 3.1
##   [19] 4.1 4.3 4.5 3.8 4.3 3.1 4.3 5.2 4.5 3.1 4.7 4.8 4.6 5.2 6.0 4.3 4.3 4.8
##   [37] 4.8 5.2 4.1 3.3 4.6 5.7 4.9 4.5 4.4 4.9 4.4 5.0 4.5 5.6 4.5 3.4 5.1 3.7
##   [55] 5.2 5.3 5.2 4.1 3.8 5.3 4.1 4.3 4.8 4.3 4.8 2.4 5.0 5.2 4.6 5.3 4.8 5.9
##   [73] 3.5 4.8 5.0 4.9 6.0 4.6 5.2 3.9 4.0 5.2 4.6 5.2 4.7 5.1 4.1 4.3 5.5 5.4
##   [91] 4.8 2.8 5.3 3.3 4.3 3.1 4.9 5.1 4.8 5.9 5.6 4.8 3.8 3.6 3.5 4.5 5.8 3.8
##  [109] 3.8 4.3 3.6 5.8 3.7 4.5 5.9 4.7 5.9 5.1 5.0 4.6 4.1 3.3 4.4 4.3 4.2 3.9
##  [127] 4.2 2.3 5.3 5.5 6.3 6.1 5.7 5.4 3.5 5.5 4.6 3.2 4.8 2.9 4.9 3.8 3.9 6.1
##  [145] 4.6 4.2 4.5 5.7 6.1 6.0 4.8 5.1 4.9 6.3 5.6 4.4 5.1 4.9 4.2 5.1 3.7 5.5
##  [163] 6.0 5.7 5.2 3.9 6.1 6.2 3.6 4.4 4.5 3.4 3.5 5.5 3.7 4.4 4.3 4.5 4.8 3.0
##  [181] 5.6 3.4 5.3 5.8 4.0 5.0 5.3 4.1 4.3 4.6 4.5 3.4 4.5 5.7 4.7 2.2 6.4 4.4
##  [199] 3.4 3.9 4.5 4.2 3.9 4.6 4.1 4.2 3.6 3.8 3.5 5.6 3.7 6.0 4.8 4.3 4.8 2.7
##  [217] 4.4 5.4 4.7 4.0 6.5 3.9 5.1 5.2 4.7 3.1 2.9 4.6 4.6 3.5 4.3 3.7 5.5 4.8
##  [235] 6.5 5.2 4.0 4.5 4.9 4.1 6.4 4.2 5.5 2.1 5.5 4.6 3.9 3.4 4.4 5.5 3.3 4.7
##  [253] 3.6 3.1 3.6 4.9 5.0 4.8 6.2 5.0 4.5 4.1 4.4 3.8 3.9 4.3 4.9 4.7 4.3 4.9
##  [271] 4.7 4.6 5.9 6.4 4.1 5.1 3.4 3.9 3.3 3.9 4.4 5.4 6.1 4.2 3.6 6.1 4.8 4.2
##  [289] 2.7 3.7 5.0 4.8 3.9 4.4 4.5 4.1 4.8 4.6 6.0 4.9 4.5 2.4 5.2 4.2 5.3 3.9
##  [307] 6.0 5.1 3.9 2.7 5.3 4.8 3.4 4.7 4.9 5.0 4.3 4.4 4.8 5.7 3.7 5.0 4.4 5.3
##  [325] 3.9 5.6 6.1 3.6 5.5 4.9 4.5 4.0 3.1 3.9 4.2 4.0 2.9 4.9 4.4 5.0 3.9 5.1
##  [343] 4.6 4.0 3.7 5.2 3.5 3.3 4.1 2.1 3.8 5.1 5.0 4.9 2.8 4.0 5.0 4.2 3.0 4.6
##  [361] 3.9 5.4 4.8 3.9 3.8 5.3 5.2 3.9 4.4 5.1 3.7 4.8 2.3 3.6 3.7 3.3 5.6 3.1
##  [379] 3.8 4.5 3.3 4.2 5.5 3.5 5.6 4.3 5.0 4.2 4.3 3.6 3.4 3.3 3.7 2.4 3.9 3.0
##  [397] 4.3 4.4 3.7 4.6 5.3 4.5 6.7 4.8 5.5 4.3 6.3 5.7 5.0 4.1 3.3 5.2 3.6 5.0
##  [415] 3.3 3.3 5.0 3.6 5.6 3.6 4.0 4.3 3.4 4.4 3.6 6.0 5.1 4.0 3.7 4.1 2.9 5.7
##  [433] 4.6 5.2 4.2 3.3 5.1 4.1 5.4 4.2 5.1 3.7 3.5 3.8 4.1 5.7 4.1 4.3 5.1 5.0
##  [451] 6.4 5.9 4.1 3.2 3.8 3.9 5.9 3.0 4.8 3.6 3.9 4.3 3.0 4.5 3.7 4.3 3.7 5.2
##  [469] 5.0 3.8 3.5 4.8 3.5 5.6 3.8 4.9 3.3 2.7 4.5 4.3 4.6 4.5 4.7 3.9 6.4 5.6
##  [487] 4.8 4.3 4.2 4.4 5.9 5.5 5.2 4.4 5.2 4.7 2.2 4.2 4.7 3.6 4.3 4.0 4.6 4.4
##  [505] 2.8 3.9 5.0 5.2 5.6 4.1 4.3 3.3 3.5 4.6 5.8 5.3 4.0 3.9 4.1 5.3 6.1 5.7
##  [523] 4.1 3.3 3.9 4.6 3.9 4.4 4.6 4.5 5.3 4.6 6.1 4.3 4.3 2.4 4.3 4.1 5.6 5.3
##  [541] 6.2 3.4 4.3 3.4 4.2 3.9 2.9 4.9 5.0 4.6 4.0 5.1 5.3 4.0 4.2 5.6 4.0 5.1
##  [559] 5.0 4.8 3.8 5.6 4.9 3.5 2.3 6.2 2.5 2.8 5.8 4.5 5.0 5.5 4.3 4.4 5.1 4.9
##  [577] 4.7 3.8 5.0 3.3 4.3 4.2 5.0 4.7 4.2 3.5 3.9 4.4 5.5 5.0 3.2 5.1 4.6 3.6
##  [595] 4.0 4.4 4.3 5.0 3.9 4.7 3.3 4.0 3.7 4.3 5.6 5.3 5.3 5.6 3.4 4.3 3.9 4.5
##  [613] 4.1 4.7 5.8 4.1 3.8 3.4 3.7 5.4 4.6 3.4 3.6 4.1 5.7 2.9 2.9 4.5 5.3 4.6
##  [631] 4.6 3.3 4.2 4.7 4.7 4.7 4.3 5.2 5.9 4.2 3.9 4.9 3.4 6.3 3.6 4.7 3.8 3.1
##  [649] 4.7 4.2 5.8 3.8 4.4 4.7 2.8 4.6 5.0 4.5 4.7 3.0 5.3 3.2 6.5 5.0 5.5 4.1
##  [667] 4.7 2.6 5.3 4.9 3.8 2.7 4.6 5.0 4.6 4.1 3.7 5.4 5.1 5.0 8.0 5.2 3.7 6.4
##  [685] 6.6 4.7 4.0 3.8 3.6 3.2 3.7 5.0 4.9 4.5 5.1 3.8 4.9 3.9 5.4 5.6 5.2 5.3
##  [703] 5.1 3.8 3.9 5.0 4.4 3.5 3.9 4.8 5.5 5.5 5.9 5.3 5.9 5.8 4.2 5.0 4.9 4.3
##  [721] 4.9 3.2 6.4 5.4 4.3 3.1 3.7 5.9 5.5 5.0 5.5 4.1 4.7 3.8 5.4 5.0 5.4 3.2
##  [739] 3.7 5.5 5.1 7.0 6.7 3.8 4.2 3.9 4.0 5.7 5.2 4.2 4.9 4.9 3.7 4.6 5.2 5.3
##  [757] 5.2 5.6 3.9 7.2 4.4 3.7 4.4 5.1 3.6 6.6 5.9 4.5 3.8 5.2 5.4 3.8 4.5 4.9
##  [775] 3.3 5.8 3.3 3.4 5.7 4.2 2.7 2.9 2.4 4.5 4.1 3.6 4.1 5.2 4.1 5.3 4.1 3.8
##  [793] 5.1 5.4 4.8 4.3 4.7 5.5 4.0 4.2 4.5 4.3 6.1 4.7 4.5 3.7 4.4 3.7 4.4 5.1
##  [811] 5.1 5.8 4.6 3.8 2.8 4.7 5.7 6.0 4.4 5.2 3.3 4.0 5.0 3.7 2.9 4.6 4.9 2.9
##  [829] 5.6 3.9 4.5 4.2 4.3 3.9 3.6 3.4 4.8 5.8 3.3 4.1 4.4 5.7 3.8 3.6 4.3 2.8
##  [847] 5.4 2.7 4.7 3.2 4.9 5.4 4.7 3.1 5.6 3.6 5.0 2.9 4.7 4.6 4.0 3.3 5.8 6.4
##  [865] 3.9 4.4 4.8 5.9 4.5 4.0 5.9 5.9 5.0 4.9 4.1 4.7 4.1 3.3 4.8 3.3 4.7 4.7
##  [883] 5.7 4.1 4.6 4.3 3.7 4.8 5.1 5.0 3.9 5.0 3.3 4.7 4.9 4.9 4.3 4.8 4.9 2.9
##  [901] 5.1 5.7 5.6 3.5 4.8 3.9 3.3 3.1 4.7 5.2 4.4 5.6 4.4 4.4 5.4 2.0 3.3 5.8
##  [919] 4.9 6.6 5.2 5.0 6.0 6.5 5.0 4.7 5.6 3.6 4.0 3.6 4.5 3.4 5.1 4.0 5.7 5.8
##  [937] 6.6 3.5 4.1 4.5 3.3 5.0 3.7 3.1 3.0 4.4 5.3 3.2 5.3 2.9 5.5 5.1 4.3 4.7
##  [955] 4.1 4.9 6.0 5.2 5.1 3.6 4.3 4.5 3.4 4.4 4.7 2.8 4.7 5.8 4.3 4.8 3.2 6.5
##  [973] 5.5 3.6 5.1 4.9 4.8 4.4 5.3 4.5 4.1 4.5 3.7 3.3 3.6 4.2 5.3 5.5 5.0 4.5
##  [991] 4.7 4.6 5.7 5.4 4.8 4.5 3.4 3.1 5.6 4.6
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
##    [1] 3.1 5.5 2.3 4.9 4.4 3.4 2.6 4.3 4.7 4.6 4.0 3.5 5.8 4.6 5.1 5.0 2.7 3.7
##   [19] 4.1 4.2 4.7 4.9 4.0 4.4 5.2 5.4 3.7 3.0 3.3 4.8 3.8 4.8 6.2 4.6 3.9 4.7
##   [37] 3.9 2.9 5.3 4.4 4.4 5.1 4.2 6.4 3.1 4.0 4.0 3.3 2.7 6.4 4.1 5.2 3.3 4.8
##   [55] 5.6 3.7 5.0 3.8 4.7 3.9 4.0 4.3 4.4 4.9 4.0 5.6 4.1 4.8 3.2 4.5 3.2 3.7
##   [73] 3.3 4.5 4.1 6.0 3.7 4.4 4.5 5.6 5.5 3.1 3.4 5.9 3.5 5.6 4.0 4.0 4.0 5.2
##   [91] 4.8 4.9 5.3 5.5 4.0 4.6 5.5 5.5 3.0 5.6 3.8 3.4 3.8 5.2 3.6 3.3 4.7 3.5
##  [109] 3.7 3.6 4.6 5.0 4.2 5.3 4.4 6.3 4.2 2.8 4.2 4.3 5.3 4.4 5.1 4.4 3.8 4.3
##  [127] 3.6 5.0 4.5 4.1 4.3 4.8 4.4 3.1 4.7 4.3 3.6 3.3 4.3 2.8 3.7 3.6 4.5 3.2
##  [145] 5.1 6.9 3.9 6.2 3.9 5.3 4.9 5.1 3.6 3.5 2.7 5.3 3.8 5.3 5.4 2.7 5.6 4.6
##  [163] 5.1 5.5 4.4 3.9 4.2 5.7 4.1 4.4 5.1 4.6 5.2 6.1 2.7 5.5 5.4 5.2 4.4 4.0
##  [181] 3.9 2.9 6.5 5.7 4.0 6.1 3.6 3.4 4.1 2.7 4.2 3.6 3.9 4.4 5.7 3.5 4.0 4.0
##  [199] 3.8 4.8 3.9 4.7 4.9 5.2 3.7 5.9 5.0 3.7 5.2 5.1 4.8 3.6 3.5 4.5 5.2 4.1
##  [217] 4.5 4.9 4.1 5.2 4.0 4.3 5.0 5.2 5.0 4.8 4.1 4.5 6.1 4.4 4.4 4.5 4.7 5.4
##  [235] 5.4 3.8 3.8 6.0 5.5 4.8 3.7 4.2 5.8 5.2 4.5 3.8 5.0 5.1 2.7 5.4 4.0 3.7
##  [253] 5.7 5.2 4.3 4.1 4.1 4.8 3.6 2.9 4.7 4.2 4.7 4.0 4.1 4.1 3.8 5.2 3.6 3.2
##  [271] 4.6 4.4 4.6 3.9 3.5 4.5 4.3 4.2 4.1 5.2 2.9 5.4 3.9 2.8 3.8 3.9 3.3 4.6
##  [289] 4.2 5.7 4.8 4.5 4.9 5.4 4.7 4.4 4.5 2.7 3.5 4.9 3.2 3.3 4.1 4.4 4.8 4.4
##  [307] 6.8 4.9 3.8 4.4 5.2 5.3 4.6 4.4 3.7 4.4 5.5 4.8 4.3 4.5 5.8 5.3 4.1 3.9
##  [325] 4.0 4.5 4.1 5.0 4.4 4.0 4.9 2.4 4.8 3.4 4.8 3.8 4.9 2.9 4.0 5.5 5.9 6.6
##  [343] 3.3 4.6 3.4 5.5 5.6 5.8 4.5 4.3 4.5 4.8 4.7 4.0 3.0 4.3 5.0 4.8 4.6 4.6
##  [361] 4.6 5.4 2.3 5.6 3.6 5.3 6.4 3.6 4.4 4.0 5.6 3.7 3.3 5.1 4.5 5.2 4.8 4.8
##  [379] 5.2 5.8 3.5 4.5 4.1 3.6 3.9 3.1 4.4 7.0 3.7 7.1 5.0 5.5 3.9 4.2 5.1 3.7
##  [397] 3.5 5.1 5.3 4.3 5.6 3.6 4.6 4.8 4.5 5.5 5.2 3.9 5.5 4.2 4.2 4.6 5.5 4.8
##  [415] 4.9 4.2 4.8 5.6 5.3 4.5 4.3 5.0 4.6 4.4 4.8 2.7 5.1 5.2 4.1 4.2 4.6 3.7
##  [433] 2.4 4.7 3.8 3.2 6.4 6.3 4.7 4.0 4.6 4.0 4.1 3.2 3.5 5.2 3.9 5.2 5.0 4.3
##  [451] 5.4 3.2 4.9 3.1 3.8 5.2 4.2 4.2 6.0 5.5 3.6 4.6 6.5 5.1 6.4 5.3 4.1 4.7
##  [469] 5.8 4.0 5.6 4.9 5.6 5.4 4.5 4.0 4.6 5.8 2.9 4.7 5.2 5.1 3.5 4.3 5.0 4.2
##  [487] 3.9 3.9 4.6 2.3 3.7 3.5 4.0 4.2 4.5 7.2 3.9 5.3 6.7 3.9 5.9 3.5 3.5 4.5
##  [505] 6.1 5.1 5.6 4.6 4.8 3.5 6.2 4.7 4.1 4.5 4.1 5.7 5.1 4.0 4.3 4.4 2.8 4.5
##  [523] 4.7 4.8 5.1 3.8 3.6 3.8 4.3 3.7 3.7 4.9 4.1 3.4 3.3 4.6 4.7 5.3 4.9 4.2
##  [541] 4.7 5.9 5.4 3.7 5.3 5.7 2.8 3.8 3.4 5.1 4.7 3.9 4.0 4.2 4.6 4.2 3.1 5.0
##  [559] 4.1 5.2 5.6 5.6 4.6 4.3 6.3 5.0 4.0 3.8 5.3 5.3 3.0 3.9 3.0 4.7 3.8 4.1
##  [577] 5.0 3.8 5.4 4.2 5.5 5.2 5.9 4.3 5.4 4.8 5.2 5.6 5.0 6.6 3.8 3.7 4.2 4.1
##  [595] 4.5 5.9 5.6 6.8 5.2 4.1 4.1 4.6 4.4 5.3 5.0 5.0 4.3 5.5 4.8 5.6 3.2 6.0
##  [613] 4.1 3.7 4.9 3.0 3.8 3.5 5.3 4.9 3.7 3.1 5.7 4.3 5.0 3.0 5.6 5.4 3.6 4.8
##  [631] 3.8 4.9 4.4 2.6 4.2 4.4 4.3 4.3 6.0 4.6 5.9 3.5 4.8 3.8 5.1 4.8 4.4 5.0
##  [649] 4.7 5.1 4.0 4.3 6.0 4.4 5.3 3.0 4.5 5.9 4.5 4.5 1.7 5.5 4.5 4.9 4.3 6.8
##  [667] 5.3 4.7 5.8 4.8 4.3 4.6 4.6 4.7 3.8 4.8 5.2 4.2 2.7 5.6 4.9 3.3 5.4 4.5
##  [685] 4.4 5.1 4.1 3.4 3.0 4.7 4.4 4.3 3.8 4.3 4.6 4.5 3.3 5.9 3.7 3.8 3.3 6.4
##  [703] 5.4 6.2 5.3 5.1 2.9 4.2 3.9 3.4 5.5 3.0 5.2 3.0 3.5 5.1 5.1 5.7 2.3 5.3
##  [721] 4.7 6.1 4.0 6.0 5.0 4.8 4.3 4.2 4.5 4.9 5.3 3.8 4.7 4.8 4.5 3.1 4.0 4.6
##  [739] 4.5 4.6 4.0 3.9 6.3 3.6 5.0 3.5 5.0 3.3 4.6 3.6 4.9 5.0 3.6 5.3 3.8 3.8
##  [757] 5.4 4.2 3.4 4.6 3.6 3.1 4.6 4.8 5.5 4.3 5.6 6.1 4.0 4.5 3.0 4.7 3.7 4.9
##  [775] 5.1 5.0 4.4 4.4 4.1 4.7 3.7 5.5 5.1 5.4 6.4 4.4 5.7 5.1 4.5 4.7 3.3 4.7
##  [793] 4.5 5.1 5.4 4.5 5.9 6.8 5.3 5.5 5.1 5.0 5.3 4.0 5.5 4.4 4.0 5.3 6.1 2.8
##  [811] 5.0 5.4 5.1 4.9 5.0 4.1 4.2 3.8 3.9 4.5 4.3 4.7 4.4 5.8 3.0 3.0 3.5 3.9
##  [829] 5.3 4.4 5.2 6.3 3.4 4.1 2.9 5.6 3.9 4.4 4.9 6.0 5.2 4.4 4.2 4.7 4.2 4.6
##  [847] 4.1 4.1 4.3 3.8 5.7 2.2 4.6 4.3 4.6 4.4 5.1 3.7 4.2 5.0 5.2 3.7 5.1 6.3
##  [865] 4.7 5.1 4.1 4.7 4.6 6.3 5.1 4.3 3.5 4.4 4.7 6.2 4.6 5.3 4.9 4.4 5.3 3.5
##  [883] 5.5 3.5 3.8 3.9 4.8 4.9 2.7 4.8 3.2 3.0 4.2 4.4 4.6 4.7 4.2 4.9 4.5 4.9
##  [901] 5.5 5.8 3.9 4.4 4.5 2.5 3.7 3.7 4.8 5.5 3.0 4.4 4.9 3.1 5.2 3.8 5.1 3.6
##  [919] 4.3 3.1 4.4 5.4 5.4 3.5 6.0 3.9 5.1 5.8 4.6 5.1 4.1 3.4 5.4 3.3 4.2 4.9
##  [937] 2.5 6.4 4.8 3.2 5.3 5.3 5.3 4.8 3.9 4.6 5.0 4.7 4.6 5.0 4.8 5.0 5.0 5.9
##  [955] 5.2 4.0 5.7 5.3 4.9 3.0 6.3 5.4 4.0 3.9 3.8 5.1 4.0 5.0 3.4 3.7 5.5 4.4
##  [973] 6.7 3.7 3.5 6.7 5.4 4.9 4.2 4.3 4.9 5.0 3.6 4.2 5.1 3.9 5.0 5.4 4.1 5.6
##  [991] 5.9 4.5 5.5 3.7 3.8 5.1 4.7 4.1 2.0 5.8
## 
## $func.thetastar
## [1] 0.016
## 
## $jack.boot.val
##  [1]  0.4756906  0.3764205  0.2281065  0.2323615  0.1002732 -0.0657971
##  [7] -0.1177285 -0.3052941 -0.3433628 -0.4521739
## 
## $jack.boot.se
## [1] 0.9080416
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
##    [1] 5.8 4.8 3.7 2.7 5.4 5.8 4.7 4.7 4.2 2.3 4.2 3.9 4.2 5.7 4.7 4.1 6.1 4.2
##   [19] 5.3 3.2 3.9 4.6 3.9 3.1 4.4 5.6 6.5 6.4 4.2 3.4 5.2 5.3 3.6 3.9 6.0 4.3
##   [37] 4.1 2.6 5.3 4.7 4.3 4.7 5.6 4.8 5.7 5.6 4.3 2.6 5.2 4.6 3.3 3.5 6.2 3.6
##   [55] 3.7 4.0 3.9 5.6 2.7 5.9 3.7 4.5 4.5 3.7 5.5 5.2 3.9 4.8 3.5 5.8 5.3 5.4
##   [73] 5.4 4.6 4.8 3.9 5.2 5.4 4.9 5.6 4.6 5.6 3.8 4.4 5.2 3.8 3.7 3.6 3.8 4.9
##   [91] 6.0 3.1 4.3 5.5 5.9 4.4 5.7 3.9 4.2 6.0 4.4 5.3 5.4 4.5 2.6 3.8 2.8 4.9
##  [109] 4.7 3.8 4.9 2.7 4.2 4.6 4.0 2.7 4.6 3.8 3.4 4.6 4.7 4.7 3.4 3.6 6.1 4.2
##  [127] 5.7 3.7 6.2 4.8 4.1 4.7 3.1 4.9 5.1 4.9 4.0 5.9 4.7 5.0 5.1 4.5 6.2 5.2
##  [145] 3.6 4.7 6.2 3.8 4.0 3.1 5.2 6.4 3.9 5.2 3.5 6.3 4.9 5.5 5.9 5.3 2.4 5.3
##  [163] 2.1 4.2 3.9 4.8 4.7 6.4 5.2 5.1 3.9 3.2 4.0 4.6 3.3 3.7 3.1 3.7 4.0 2.9
##  [181] 4.6 5.7 3.4 4.0 5.3 5.0 3.9 5.5 5.7 3.3 6.0 4.0 4.4 3.6 4.0 4.6 3.9 4.1
##  [199] 4.4 3.3 3.4 5.1 4.2 4.2 3.6 3.9 3.4 3.8 3.5 4.6 4.2 4.6 3.4 5.2 4.6 5.1
##  [217] 4.5 5.9 5.1 3.6 3.6 5.5 4.1 3.4 5.2 5.3 5.1 6.0 4.6 4.6 3.6 4.9 4.0 5.5
##  [235] 3.4 3.9 4.3 5.8 4.3 5.6 6.4 3.9 2.7 3.9 4.4 5.1 7.2 4.8 2.8 4.9 5.5 6.1
##  [253] 4.9 3.6 4.3 4.1 5.6 3.2 4.6 5.3 5.6 3.4 3.4 4.5 3.8 3.7 4.7 4.4 4.7 5.4
##  [271] 4.0 5.1 3.6 5.7 4.2 3.8 4.9 4.9 6.0 4.5 6.1 4.7 5.1 4.3 4.2 5.1 3.3 4.4
##  [289] 5.9 4.8 5.0 4.2 4.2 5.8 4.8 4.6 6.3 3.2 3.8 4.6 4.2 4.8 4.7 4.0 4.6 5.9
##  [307] 5.4 4.3 5.5 4.3 4.2 6.0 4.3 4.7 3.8 3.6 3.3 4.3 3.9 6.4 5.5 4.4 4.6 5.7
##  [325] 4.8 4.3 4.8 4.6 5.2 5.5 4.4 4.5 4.5 5.4 4.5 5.1 4.2 4.7 4.3 4.9 5.1 5.4
##  [343] 4.0 5.9 4.4 4.7 4.6 4.3 4.9 4.3 3.2 4.7 6.0 3.8 4.1 5.4 4.1 5.1 5.7 4.2
##  [361] 3.4 5.3 3.3 3.6 5.3 5.9 5.9 4.0 3.7 5.2 3.8 3.5 4.7 2.6 4.5 4.6 4.1 3.2
##  [379] 4.4 6.2 3.7 5.3 4.2 4.3 4.8 4.3 4.3 4.5 6.5 5.0 4.2 6.2 3.3 4.7 6.3 5.2
##  [397] 4.6 4.5 4.7 4.8 3.7 6.1 3.9 4.5 3.3 3.2 3.7 5.2 5.0 3.5 4.2 4.1 4.2 4.4
##  [415] 3.8 4.0 5.8 3.6 4.9 3.7 2.0 4.1 3.7 4.1 4.5 5.1 3.5 5.4 4.2 3.7 4.1 2.3
##  [433] 4.6 4.1 5.0 3.8 3.9 5.9 5.5 3.9 3.6 5.1 4.2 5.9 5.8 3.7 4.8 5.8 5.0 4.9
##  [451] 4.1 4.1 3.9 4.7 5.7 4.8 5.1 3.9 4.5 5.6 4.0 4.7 4.9 4.3 5.9 4.3 4.7 6.0
##  [469] 3.9 4.1 7.0 5.0 3.2 4.5 5.7 4.4 3.8 6.5 4.5 3.0 4.5 6.2 5.4 3.3 4.0 3.2
##  [487] 3.5 3.8 4.8 4.6 5.1 4.7 3.8 5.1 4.4 4.6 4.1 4.5 4.7 5.0 4.7 4.3 6.0 4.5
##  [505] 6.1 3.6 5.1 4.5 4.8 3.6 4.3 3.9 4.9 4.7 5.0 4.2 4.5 4.4 6.4 5.1 5.0 4.2
##  [523] 4.8 5.1 3.6 4.2 4.4 3.5 4.1 3.9 3.2 3.6 5.8 3.1 5.9 4.2 4.8 5.7 4.6 3.3
##  [541] 5.0 5.4 4.1 4.1 5.3 3.8 4.6 4.4 5.1 4.2 4.8 5.5 3.2 4.4 3.9 4.7 4.3 5.2
##  [559] 5.0 4.3 3.5 4.6 3.6 3.5 3.1 3.9 2.7 4.4 3.5 4.5 5.2 4.1 5.3 4.6 5.2 3.5
##  [577] 4.5 4.4 5.7 4.1 4.2 4.4 5.2 4.9 6.6 4.0 6.3 4.1 3.0 4.1 5.7 5.3 5.1 5.6
##  [595] 5.0 4.8 2.9 5.1 4.8 4.6 4.2 5.0 5.3 4.9 5.7 3.5 5.6 3.2 3.8 5.0 3.5 3.0
##  [613] 5.1 2.4 6.4 4.6 3.1 3.2 5.1 3.7 5.6 5.9 5.9 4.3 5.9 5.3 5.1 4.6 4.0 5.6
##  [631] 4.2 2.4 5.1 4.4 4.4 4.4 2.7 1.9 5.0 3.7 4.5 4.7 4.5 4.3 3.9 4.6 5.5 4.2
##  [649] 4.6 3.8 4.2 4.6 3.3 4.3 2.8 3.3 3.6 4.1 3.4 3.2 5.6 5.2 4.4 6.4 5.7 5.2
##  [667] 3.7 5.7 4.0 4.3 3.9 5.6 3.0 3.6 4.6 5.6 4.8 4.0 3.9 3.4 4.3 4.8 4.6 5.4
##  [685] 5.3 3.4 3.8 2.9 5.2 5.9 3.7 2.4 4.3 5.6 3.1 5.1 4.7 3.6 3.8 5.8 2.6 3.8
##  [703] 4.3 4.3 4.9 3.1 3.9 3.6 5.5 3.7 4.0 4.6 4.7 7.5 3.8 4.7 4.4 4.6 4.5 3.2
##  [721] 4.7 4.8 5.9 4.5 4.1 4.9 5.5 4.2 3.9 3.9 6.3 4.2 2.9 5.6 4.6 4.6 4.7 3.9
##  [739] 3.5 5.4 4.5 3.8 6.0 3.7 3.4 3.7 4.2 4.7 4.3 6.1 3.3 4.8 6.0 3.7 5.0 5.6
##  [757] 5.0 4.4 4.5 5.7 4.9 3.6 4.3 4.7 2.9 2.9 4.9 3.7 4.4 3.2 5.2 5.0 3.5 6.3
##  [775] 6.1 3.7 3.6 2.1 4.4 4.5 3.2 5.7 4.8 3.1 5.0 6.0 4.1 5.1 4.6 5.3 3.7 3.5
##  [793] 4.7 4.4 5.1 6.5 5.1 4.6 4.3 3.3 6.6 3.9 5.1 3.3 5.1 5.9 5.3 3.9 4.4 4.1
##  [811] 4.1 3.9 4.6 3.9 4.7 3.6 4.3 5.0 2.2 6.1 4.2 4.1 5.4 5.4 6.9 4.8 4.2 4.9
##  [829] 4.9 5.4 5.6 6.3 2.9 3.2 3.9 4.4 4.4 4.2 3.8 4.2 4.5 4.2 3.8 5.3 3.8 5.9
##  [847] 4.1 3.7 4.3 5.6 3.7 3.2 5.9 4.2 3.8 4.5 2.6 4.8 4.9 4.2 3.3 3.8 4.5 4.1
##  [865] 5.6 5.9 5.5 6.0 4.3 5.6 4.2 4.0 4.8 3.1 3.5 4.0 5.1 5.4 5.9 4.7 5.2 3.4
##  [883] 2.5 4.2 4.2 3.0 5.4 4.4 4.8 4.7 4.7 4.1 5.8 3.3 2.2 4.4 4.4 5.6 4.9 4.4
##  [901] 4.6 7.1 4.0 3.7 3.5 4.5 5.7 3.7 5.6 2.6 5.1 4.9 4.7 4.6 4.5 3.6 4.1 2.9
##  [919] 4.0 5.0 4.8 4.6 4.9 4.7 3.9 3.9 5.0 2.9 4.0 4.8 4.5 6.5 4.0 2.9 5.1 4.7
##  [937] 4.5 4.7 3.8 3.6 6.0 4.7 4.3 5.6 6.2 4.5 2.8 4.9 5.0 5.2 4.4 5.2 4.2 4.6
##  [955] 5.5 6.8 4.4 3.0 5.7 4.3 3.2 5.7 5.3 3.1 4.0 4.3 4.4 4.5 4.7 4.9 5.6 3.6
##  [973] 3.9 4.7 4.3 3.7 5.7 5.0 4.9 5.4 6.0 6.2 4.7 3.3 5.3 5.7 5.2 5.1 3.7 3.2
##  [991] 3.5 4.2 1.9 3.7 3.5 3.3 4.1 5.4 4.8 4.3
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.6 5.5 5.3 5.2 5.0 5.0 4.8 4.7 4.6 4.5
## 
## $jack.boot.se
## [1] 1.071634
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
## [1] 0.2161208
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
##   6.315802   8.363126 
##  (2.753143) (3.794443)
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
## [1]  0.55136032 -0.05812794  0.60704378 -0.38495415  2.10281755  0.37736439
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
##    [1]  0.1252792512 -0.6050682561  0.5707592092  0.4222036692 -0.0417584022
##    [6] -0.2077233637  0.7319694979  0.2266711647  0.5362952217  0.3220284171
##   [11]  1.0229156574  0.3610961167 -0.1457095551  0.0883964245  0.1803518101
##   [16]  0.0827228296 -0.6786614531  0.0302633075 -0.4365675772 -0.1224341418
##   [21]  0.1766727694  0.1405204848  0.1568660965 -0.2485356982 -0.5917868727
##   [26]  0.1024718009  0.2713288573 -0.0176743882  0.9221695369  0.5295822677
##   [31] -0.1758539486 -0.0124918310  0.5188546238  0.2664736696 -0.4590882083
##   [36] -0.4556262939  0.2103251537  0.0356281248  0.6242062387  0.1745586329
##   [41]  0.2285538768  0.0289399466  0.6781195132  0.8043340320  0.0304091651
##   [46]  0.6699515961  0.4553907797  0.4528600832  0.2962337363  0.1312301456
##   [51]  0.5752409251  0.5576627173 -0.3143777929 -0.3101658817 -0.5518489689
##   [56]  0.3993154825 -0.0236785642  0.0096770556  0.1590505828  0.8484822680
##   [61]  0.1103721755  0.4009961711 -0.1449548171  0.0033613420  0.7323044344
##   [66]  0.3398296618  0.1172153168 -0.4163560770 -0.1253571045  0.1897483829
##   [71]  0.2397282592 -0.1274132078 -0.6132568516 -0.5097669894  0.3524038503
##   [76]  0.1020824451  0.5985835096  0.0448579689 -0.1085944778  0.0672703539
##   [81] -0.6034999462 -0.1092232340  0.2932096153 -0.0083688116  0.0056300398
##   [86] -0.2269149278  0.0350368190 -0.0248151385  0.5858695970  0.1296324709
##   [91]  0.3034259846  1.1386511284  0.7169776709  0.1124803635  1.4787732487
##   [96] -0.6868106765  0.4243340344  0.1465395307 -0.4773346918 -0.0224489412
##  [101]  0.1457328758  0.9003158931 -0.0364015602 -0.1094224167 -0.0944076299
##  [106]  0.0350797784  0.5576808005  1.2387057831  0.5708413660  0.3359646094
##  [111]  1.0579295175 -0.3588899592  0.0716313295  0.2530341822  0.9293209927
##  [116]  0.1556110018  0.2256591929  0.1988527980 -0.2297109072 -0.2620758387
##  [121]  0.3941133461  0.4900866981 -0.4031078308  0.2507724636  0.3434151080
##  [126]  1.0136749634 -0.1457190291  0.1315871725  0.1784082416  0.4358598737
##  [131] -0.0192404149  0.3221744382  0.0775545531  0.3057529676  0.4024228636
##  [136]  0.6348391526 -0.1908874343 -0.4073775891  0.1498741977  0.2250082638
##  [141] -0.2408200404  0.3475560598 -0.0267447455  0.1360458121  0.0176305476
##  [146]  0.1873298561  0.1596658473  1.0764032569  0.0781093896  0.2398767675
##  [151]  0.1281753848 -0.7178650454  0.2695745972 -0.1769260779  0.2431268632
##  [156] -0.2966914870 -0.2128984590  0.3092572830  0.1752330793  0.8615279177
##  [161]  0.2601415705  0.1044310303  0.5304806455  0.4401903870  0.7736211262
##  [166] -0.6932022668  0.5315365353  0.3691448777 -0.0807269432  0.3500238624
##  [171]  0.3044413706 -0.3178881980 -0.2502721906  0.0291376006  0.9067381324
##  [176]  0.1501169014  0.0726108692 -0.3834148370 -0.2272685039 -0.4776693652
##  [181] -0.4098023777 -0.3568402282  0.5977390999  0.5335322964  0.2243113282
##  [186]  0.8081544963  0.5764018455 -0.5244260527 -0.1065304325  0.2579715013
##  [191] -0.1682347772 -0.3875627104 -0.3289581573  0.5492879718 -0.4389221284
##  [196]  0.4064532516 -0.0478038638  0.3744513005  0.4549753582 -0.3217892230
##  [201]  0.5988467185  0.8360229047  0.4225688012 -0.0156057147  0.4035955149
##  [206] -0.0206144002  0.5095913961 -0.7858263680  0.3428340622 -0.1654911292
##  [211] -0.2702593530 -0.1041622732  0.5319235631 -0.7525661061 -0.2802958064
##  [216]  0.3046992425  0.3214016661  0.4315628814  0.0979487596  0.1222805113
##  [221]  0.5668756927  0.7253945200 -0.1871141771  0.0470318137  0.1128592418
##  [226]  0.6465184907  0.4616026449  0.6121269456 -0.2714848602  0.3367722655
##  [231]  0.0102788689  0.2135996194  0.9908205973  0.1759394527  0.6357648361
##  [236]  1.1663185047  0.1502537812  0.2354172297 -0.1869840130 -0.4752404223
##  [241]  0.3004414642  0.8425345934  1.6457302752  0.3093424284  0.4288141788
##  [246] -0.4955065270 -0.0153924743  0.0500514835  0.2991651543  0.1285827161
##  [251] -0.0504573596 -1.0305407227  0.0445814810  0.0840473352 -0.1084634784
##  [256]  0.2557319293 -0.0009736092 -2.2022998464 -0.0098106175  0.4111333210
##  [261]  0.4377364324 -0.1185060607  0.7723523516  0.5433216908  0.6834303830
##  [266] -0.2721140226 -0.7360273450  0.2111952761  0.2658886821 -0.1322398598
##  [271] -0.0612219514 -0.1013524334  0.3579587034  1.5350374687 -0.9257985778
##  [276]  0.2891612887 -1.2428423746 -0.3241049901  0.1813641952  0.1388255014
##  [281]  1.1215586328  0.3322348367  0.5798602387  1.2585718371  0.4178651611
##  [286]  0.5737604671  0.4809443282  0.2227502694 -0.0987893143  0.4654034949
##  [291] -0.1029442439  0.0002784182  0.0970056382  0.7833954070  1.0898464573
##  [296]  0.2588907325  0.3376780437  0.0441600358 -0.0433001232  0.2048022062
##  [301]  0.1337842893  0.0616712381  0.2930490002  0.1510575633  0.1290648314
##  [306] -0.1134391652 -0.1047348977 -0.1669393502 -0.1097678781 -0.2565544505
##  [311]  0.2203519281  0.1410697292  0.1681623516  0.2017398765  0.2324330290
##  [316]  0.3174125573  0.3821325217 -0.0840999439 -0.0237362606  0.6433066925
##  [321]  0.6200583969 -0.0360709435 -0.8547208363  0.1318041840  0.1186167670
##  [326] -0.1558362034 -0.7316774572 -0.8861867330  0.7522827986  0.3350218441
##  [331]  0.2510964665  0.1970833618  0.0808614904 -0.0764658726  0.7884039261
##  [336]  0.7705362191  0.1935594705  0.3358458057  0.9304316730  0.6481477017
##  [341]  0.4186912770 -0.1726505469 -0.4937517312  0.1462487983 -0.1931696211
##  [346]  0.4224406687  0.6065765857  0.6965586569 -1.2433374740  0.5659879712
##  [351]  0.8608881018  0.6009184222 -0.2630752994  0.2917061454  0.2213462566
##  [356] -0.1091557873  0.2086928918  0.0804431469  0.4621545652  0.4608535279
##  [361]  0.6371443732  0.8154187155  0.5367495351 -0.5918586121 -0.7664625799
##  [366]  0.6018120258 -0.2031960424 -0.1031171394 -0.1738296226  0.4931547257
##  [371]  0.2165289500  0.7907781059  0.3941133461  0.5327437643  0.1716391293
##  [376]  0.2560268357  0.3297193343 -0.1861500762  0.3138887829  0.4692423514
##  [381]  2.3487293311  0.0997393159  0.9006113916  0.6736832056  1.3474436407
##  [386] -0.1142007278  1.6383301351 -0.2453046737  0.2094952544  0.0813456666
##  [391]  0.4730091537  0.3977133690  0.3465853287  0.4430930046  0.2929818311
##  [396]  0.1560330449 -0.3063323996  0.1980398083  0.2624253944  0.2702254001
##  [401]  0.1897219789  0.5411197425 -0.0236785642  0.4649996301 -0.2040780406
##  [406]  0.4428446097  0.5412455364 -0.0330078607  0.1854577899 -0.3429870265
##  [411] -0.4859841991  0.4143303728  0.2025173043  0.5250113953 -0.4043310911
##  [416] -0.1817977836 -0.0124918310  0.4755326671 -0.5609300636  1.8075133356
##  [421]  0.0344086076 -0.2845000486  0.1123970939  0.4898029243  0.5058665058
##  [426]  0.3584919237  0.4429625761  0.1565471213  0.1628477408 -0.1447282311
##  [431]  0.1852648784  0.1927341817 -0.4229081597  0.0663195883  0.1167281541
##  [436]  0.0208049717  0.0302633075 -0.0806846279  0.3980409105  0.2522173413
##  [441] -0.1226046071  0.5116658616  0.3404203317 -0.0873646920 -0.1744694941
##  [446]  0.6662023577 -0.1264626466 -0.1549045389 -0.3556056078 -0.4943929688
##  [451] -0.2139905122 -0.4487074802  0.6587812163  0.6338370894 -0.1041045101
##  [456]  0.2563220498 -0.1065909708  0.1238511905 -0.3284094348  0.2237662755
##  [461]  0.0540503657 -0.1200636775 -0.0611113330 -0.1908874343  0.1487101360
##  [466] -0.0637561079 -0.3859366253  0.9014314213 -0.5777792192  0.1652021204
##  [471]  0.2142555581  0.8763199769 -0.1188188807  0.4065667404  0.8326513536
##  [476]  0.3174125573  0.7342006971  0.2233831699 -0.2510490354  0.3489444463
##  [481]  0.5014234612  0.3417804025  0.0216841082  0.4698360317 -0.2109833584
##  [486] -0.2807339008  0.3046644958  0.1757171495 -0.0919281202 -0.0156441018
##  [491]  0.0222507087  0.7489217721  0.2935980154 -0.0742524524 -0.1924286805
##  [496] -0.2764909299  0.3604183211  0.2501080120  0.4448971470  0.6365239903
##  [501]  0.5470684341  0.6209482621  0.2443386515  0.1646411952  0.3869338188
##  [506]  0.4213652009  0.6959831392 -0.2681799995  0.5923823518  0.4220593003
##  [511] -0.2174770672  0.0940481652  0.2559364772 -0.5499069319 -0.4560607240
##  [516]  0.8558325504 -0.4084421403 -0.6371250313 -0.3528200101 -0.1231503968
##  [521]  0.4144459248 -0.3776953110  1.0915104550  0.8132444777 -0.0857219000
##  [526] -0.5037558302  0.3851900438  0.5197781521  0.2201619650  0.0797178531
##  [531] -0.0006012210  0.8582074673  0.4531165960  0.3670913596  1.2329262954
##  [536]  0.2927036231 -0.4256723719  0.0181843044 -0.4465445933  0.0404707277
##  [541]  1.1033914183  0.8356296853  0.1758591483  0.4738805910  0.1747722455
##  [546]  0.7228614131 -0.4714752655 -0.0438879574  0.3149255427  0.0820154623
##  [551]  0.4171310860  0.3544691456 -0.0961445777  0.3645966227  0.5096111721
##  [556]  0.5541120518 -0.3396723154  0.2205692772 -0.2495316190 -0.0826430752
##  [561] -0.0004407475 -0.1338976572 -0.1558244033 -0.7308681801  0.4871334799
##  [566] -0.3326775048  0.4195547672 -0.5511029831  0.0870122779 -0.6685087842
##  [571]  1.3368358405 -0.1550062993  0.2487049474 -0.4536960474  0.4444869953
##  [576] -0.6998731177 -0.2632173497  0.4395446381 -1.5330791739  0.5899931181
##  [581]  0.3937602544  0.0960232753  0.5277601179  0.0804183573  1.1068687091
##  [586] -0.5346440235  0.0043901997  0.1880988877 -0.6628595580  0.7493256288
##  [591]  0.2411841597 -0.5352010020 -0.0442875962  0.1485435330  0.1245745994
##  [596]  0.2815686952 -0.3405139731  0.8216622468  0.4639737538  0.5534711641
##  [601]  0.0979548758  0.6033517818 -0.4616959444 -0.4641913350 -0.2233432293
##  [606] -0.0951389914 -0.5450810515 -0.0547073986  1.3206872616 -0.0311505062
##  [611]  0.0107161203  0.5941296443  0.4815489112 -0.4250655484 -0.8679365671
##  [616]  1.3138891523 -0.0848551480  0.4423318211  0.1915160105  0.3060225182
##  [621] -0.1413839317  0.1081520362  0.6384014486  0.4224357754  0.2772533538
##  [626] -0.1799542426 -0.7108957607 -0.7293417084  0.4429625761  0.6326544514
##  [631]  0.0297941701  0.5561084511  0.8070585824  0.5079087435  0.4087078145
##  [636] -0.3547118790 -0.3643651434 -1.1173160798 -0.1075511754  0.0001517658
##  [641]  0.7080370572  0.6454710649 -0.2117538357 -0.5503535324 -0.1614555315
##  [646] -0.1439312542  0.1794387396 -0.1327500879 -0.0697074980  0.9935342614
##  [651] -0.0614718543  0.1857068454  0.8937463160  0.4445610037  0.1500134088
##  [656] -0.3473166208 -0.2721332718  0.1924871742  0.1836447026  0.3158217284
##  [661]  0.3124695161  0.3523754459 -0.4418107675  0.4498180524  0.2844023336
##  [666]  1.4212769692  0.8773357871 -0.1813535777  0.6518955310  0.7314348376
##  [671]  0.0900507493  0.8444619774 -0.1451486316 -0.3709450814 -0.0249994359
##  [676] -0.0732681144  0.2410262425 -0.1120521122  0.2467921801 -1.1622991305
##  [681] -0.7608839326  0.4652997177 -0.3547118790 -0.0460179138  1.0033081310
##  [686]  0.1782378137 -0.1764973755  0.8730418778 -0.3933446980  0.3711047825
##  [691]  0.1163889901 -0.3192902670  0.2423017482 -0.0691149912 -1.1190959235
##  [696]  0.0383061315  0.5780156726  0.1560309748  0.2852616148 -0.2810267130
##  [701]  1.4065611398 -0.1179442721  0.9645737947  1.1144839219  0.6447715953
##  [706]  0.2184138219 -0.2135392067  0.2639140275 -0.3709450814  0.3358458057
##  [711]  0.0261666346 -0.2268197238  0.4102520873 -0.2237224849  0.2074650779
##  [716]  0.6360643110  0.1366483169  0.6282150724  0.2394045340  0.3055753435
##  [721]  0.7400229593  0.0610435338 -0.2950687481 -0.1057842065 -0.2874634296
##  [726] -0.5656513467  0.4201177453  0.1781549104 -0.2142618486  0.4079070952
##  [731]  0.3432642997  0.5048336139 -0.0311121765  0.5549122250 -0.2907710140
##  [736] -0.1400891567 -1.0864139775 -0.3745110109  0.1003655423  0.5558950123
##  [741]  0.6036887069 -0.5801196887  0.4569132570  0.3413492881  0.1397835352
##  [746] -0.4256723719 -0.1868617595  0.5401783388  0.3877448977  0.3557233707
##  [751] -0.3259507166 -0.1823069180  0.2515790045  0.9978745530 -0.0029599440
##  [756]  0.1378584135 -0.5077064126 -0.4900641973  0.9111876699  0.9221985105
##  [761]  0.1091834103 -0.5484377944 -0.0466865063 -0.5550778024 -0.3321466943
##  [766] -0.4824335150  0.6401679368  0.1306571894  0.0533582740  0.1385348263
##  [771]  0.7644189081  0.4350257461  0.4540317228 -0.2113536742 -0.3517090666
##  [776] -0.5237933611  0.3859614383  0.1920064934  0.2113633121  0.2050898045
##  [781] -0.0101388653  0.4435186550 -0.2159749472  0.1554726717  1.7483288806
##  [786]  0.6771980768  0.6224202622 -0.2658192932  0.9605425664  0.5367495351
##  [791]  0.4055417500  0.5569651528 -0.0533654901  0.1909711706  1.0688291864
##  [796]  0.2717390022  0.3252581193 -0.3515731615 -0.0556166256  0.5067558926
##  [801]  0.8505916490 -0.2286610743  0.3474065736 -0.5703611374  0.7942897810
##  [806]  0.7485784945 -0.2037378229 -0.1061423065  0.0810384626 -0.2008857789
##  [811] -0.1953851620  0.1661122046  0.5267205840  0.4078823751 -0.5829781969
##  [816] -0.4717316938  0.1988894913  0.5328651995  0.1639628540  0.8544720645
##  [821]  0.3456038734  1.0412493806  0.4085063940  0.1827389414  0.1011368431
##  [826]  0.2719778429  0.2706072459  0.3122315680  0.1839726092  0.1605613274
##  [831] -0.7608839326  0.7816434844 -0.1339545338  0.2607571672  0.1264453169
##  [836]  0.1405307394 -0.5504876137  0.3185286904 -0.6467117321  0.2445011954
##  [841]  0.0348704348  0.4659731821  0.0136740264  0.1184888112 -0.8998640448
##  [846] -0.1372018224  0.7025910187 -0.3152750288  0.4077317893  0.9755860255
##  [851]  0.3993449195  0.4524431964  0.6399312347  0.3033860905  0.1299793703
##  [856] -0.8762254676  0.7417315041  0.7693324910 -0.0345692527  0.7364071266
##  [861]  0.3223415482  0.2266638785  0.1267312294  0.5196771364  0.7580671844
##  [866] -2.3462396998  0.4770805400  0.2641801840  1.6925687914  1.0135448112
##  [871] -0.0761105005  0.4430032514 -0.1113906785  0.4426307966  0.1797593786
##  [876] -0.2838464239  0.5730259442 -0.5044033417 -0.2693255979  0.0801217047
##  [881] -0.2383932320  0.8345009230  0.5314169750  0.2124840455  0.2036502929
##  [886]  0.0170331352  0.5406929028  1.0229156574 -0.3094944643 -0.2760430386
##  [891] -0.2494791993  0.3579171996 -0.4512337530 -0.2426548223 -0.4723446919
##  [896]  0.2899084538  0.5313141996  0.5643293544  0.6459022360  1.2853751413
##  [901]  0.1126448504  0.3807076494 -0.3767425673 -0.0706336539 -0.2542221548
##  [906]  1.0284234581  0.9525858198  0.8601809896 -0.2315597781  0.0004641584
##  [911] -0.0517114238  0.1287442938 -0.3528428848  0.0809271116  0.0300668226
##  [916] -0.5017523976  0.0176109832  0.5335322964  0.2118428119  0.2209404515
##  [921] -0.3404646851 -0.3690380540 -0.2441793806  0.6853659675  0.2057011232
##  [926]  0.9123059298  0.1588679246  0.1994275283 -0.1804953329  0.5922932036
##  [931] -0.5082572746 -0.3030188372  0.5221150392  0.4666061863  0.2696359834
##  [936] -0.2655476603  1.1870418986 -0.1988797821  0.1100213511 -1.0125596587
##  [941] -0.1100325200  0.0874825990 -0.1037453318  0.3471556488 -0.0486226748
##  [946] -0.8763827862 -1.0352027686 -0.4340211167  0.2287682177 -0.3465441495
##  [951] -0.7857061651  0.2956694914  0.3536811047 -0.2940707299  0.7765266181
##  [956]  1.1814629439  0.4270308334  0.3127258991 -0.3540710086 -0.2567589802
##  [961]  0.2354504088  1.1193121459 -0.3660958590  0.4688576491  0.4523129881
##  [966] -0.2757722103  0.4533516658  0.3661061720  0.6282965519  0.0205239149
##  [971]  0.2482902636  0.2211414064  0.5903978992  0.1889495900  0.3083234515
##  [976]  0.3496518206  1.2775069800 -0.8982289428 -0.0338859426  0.1037692037
##  [981]  0.0778021680 -0.6520136382  0.2038082467 -0.0268419464 -0.4042478617
##  [986]  0.7606326360  0.2229311299  0.2008081128 -0.0724668212  0.7666409709
##  [991]  0.2332531232  0.3455997394  0.1539223084 -0.1629435705 -0.6607141210
##  [996] -0.1380060229 -0.1345525026  0.9285580724  0.7855303102 -0.0437577027
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
```

```r
fit2
```

```
##       mean          sd    
##   0.75518725   0.28977874 
##  (0.09163608) (0.06479044)
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
## [1] -0.04592120  0.25160638 -0.06160416 -0.12488372 -0.55407859  0.62160655
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
## [1] -0.0063
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9344322
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
## t1*      4.5 0.008108108   0.9532164
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 2 3 6 8 9 
## 2 1 3 2 2
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
## [1] 0.0028
```

```r
se.boot
```

```
## [1] 0.9481198
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

