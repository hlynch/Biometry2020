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
## 2 4 6 7 8 9 
## 1 2 3 1 2 1
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
## [1] 0.0011
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
## [1] 2.712659
```

```r
UL.boot
```

```
## [1] 6.289541
```

Method #2: Simply take the quantiles of the bootstrap statistics


```r
quantile(xmeans,c(0.025,0.975))
```

```
##  2.5% 97.5% 
##   2.7   6.3
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
##    [1] 3.6 3.8 6.2 5.1 5.1 5.1 5.6 5.5 4.9 6.6 4.5 3.8 3.5 4.4 3.5 5.3 5.0 3.3
##   [19] 2.7 5.1 2.8 4.4 4.2 6.4 4.6 4.7 3.5 3.0 2.3 3.7 5.1 4.9 4.4 3.7 4.1 5.5
##   [37] 4.4 5.0 5.9 4.5 3.4 3.9 5.2 4.0 5.1 5.0 5.6 5.9 5.3 4.3 5.9 4.6 3.3 4.9
##   [55] 4.7 5.2 3.2 6.7 4.5 5.1 3.0 6.5 4.9 2.9 5.6 4.3 3.5 2.6 3.9 4.8 4.0 5.7
##   [73] 5.2 3.1 3.1 3.6 5.0 4.0 2.9 4.6 3.5 4.4 4.1 6.9 4.6 5.8 3.6 4.0 4.5 6.1
##   [91] 5.2 4.6 3.2 6.2 4.8 5.2 3.7 3.3 4.5 3.7 6.4 4.6 5.6 7.6 5.1 4.8 3.9 5.5
##  [109] 5.3 4.6 4.3 4.8 4.2 3.6 5.4 3.1 4.1 5.2 3.6 5.2 4.5 3.7 4.3 3.4 3.3 5.1
##  [127] 3.9 4.5 4.0 4.8 5.0 4.8 4.3 3.4 3.4 3.9 4.7 6.3 5.9 4.2 5.5 6.4 5.7 5.2
##  [145] 4.5 4.0 6.0 3.9 5.2 3.3 5.3 4.0 4.3 4.9 5.4 1.9 3.1 4.5 5.3 3.3 4.3 4.5
##  [163] 3.4 3.8 3.3 4.9 4.1 5.1 4.0 3.8 3.6 5.3 4.5 4.9 4.4 4.8 4.3 5.1 4.4 3.8
##  [181] 6.9 5.0 3.8 5.0 4.8 4.2 5.3 4.4 2.8 5.7 5.2 5.4 4.7 3.9 4.6 6.7 4.5 4.7
##  [199] 4.2 4.0 6.9 6.2 3.7 3.1 4.6 1.9 4.6 3.9 5.7 3.9 3.5 5.1 5.6 4.9 5.6 3.9
##  [217] 2.2 4.6 2.7 4.5 4.5 4.8 4.9 5.4 4.3 4.6 2.9 4.6 5.1 5.3 5.4 3.3 5.2 5.0
##  [235] 3.8 4.0 3.4 5.5 5.3 5.7 5.2 5.8 6.1 4.7 5.5 3.3 4.9 4.3 3.8 4.5 4.8 5.6
##  [253] 6.2 5.6 3.9 3.3 5.1 4.7 5.7 3.4 4.6 5.3 5.6 4.2 5.3 4.4 3.8 4.6 4.4 4.8
##  [271] 4.4 4.9 4.0 2.8 2.9 3.9 4.1 4.5 3.2 4.6 5.7 4.1 3.5 6.3 3.7 5.2 5.0 4.2
##  [289] 4.5 3.7 5.1 4.1 3.0 5.2 4.0 5.4 3.8 5.1 3.2 5.9 6.7 4.1 4.2 3.6 4.3 4.3
##  [307] 6.2 5.2 2.8 6.2 4.3 4.5 4.1 4.3 3.7 4.1 4.3 3.5 4.0 3.6 3.4 3.3 5.1 4.1
##  [325] 6.0 4.3 5.3 4.7 4.2 3.8 3.5 4.4 5.0 5.2 4.8 3.3 4.7 4.7 3.6 4.3 6.2 4.3
##  [343] 6.2 4.6 5.4 4.0 4.7 6.1 6.2 6.7 5.3 5.4 4.2 2.6 4.7 5.3 4.9 4.8 4.6 5.5
##  [361] 3.8 3.7 4.6 5.1 4.0 3.1 5.7 4.1 3.6 3.9 4.9 4.6 5.3 3.3 6.8 3.2 5.5 4.2
##  [379] 3.3 5.1 4.5 6.0 4.6 3.0 4.1 5.0 3.5 5.9 5.0 3.8 4.7 4.6 4.9 4.3 5.2 5.1
##  [397] 2.8 4.1 4.8 4.8 4.5 5.2 5.5 3.1 4.5 6.1 6.2 5.2 6.3 2.6 4.4 6.0 5.3 5.5
##  [415] 4.3 5.0 3.5 4.0 5.3 3.4 4.6 5.1 5.3 4.4 4.1 4.1 5.1 4.9 3.6 5.4 4.2 4.1
##  [433] 3.3 4.3 3.9 5.0 3.2 4.2 2.5 5.5 4.9 4.3 4.6 3.3 3.3 4.4 5.1 3.8 5.4 4.5
##  [451] 3.6 4.7 5.6 7.2 4.3 3.7 4.7 4.2 4.9 4.3 5.3 4.5 6.2 4.1 4.5 5.1 3.8 4.8
##  [469] 3.9 3.8 4.0 4.1 7.2 4.3 4.5 5.2 3.7 3.0 3.9 4.8 4.0 2.7 6.9 5.5 5.4 3.6
##  [487] 5.4 4.6 2.4 5.4 6.1 5.0 4.8 4.9 3.8 3.8 4.1 5.3 4.9 7.1 4.5 3.9 4.3 4.4
##  [505] 5.7 5.1 4.3 6.2 5.0 4.9 4.5 4.4 4.8 5.6 5.2 5.0 4.3 4.6 3.6 4.5 4.8 3.3
##  [523] 4.9 4.8 4.3 3.0 2.8 4.7 3.9 4.7 5.4 5.5 4.1 5.4 3.4 4.8 3.9 4.5 5.9 6.3
##  [541] 5.1 4.8 5.6 4.1 3.7 4.7 5.5 4.9 6.0 3.1 4.8 3.7 6.0 5.4 4.8 6.0 4.0 3.2
##  [559] 4.0 4.1 2.8 3.1 5.1 4.1 4.6 4.4 4.0 4.0 2.9 3.4 3.9 5.3 5.4 5.7 4.3 3.4
##  [577] 2.9 4.8 5.3 4.4 5.8 5.0 2.8 4.7 3.3 4.3 3.8 5.2 5.5 3.6 4.3 4.9 5.0 5.1
##  [595] 5.9 5.7 2.9 3.6 4.1 4.0 5.0 4.4 5.7 3.1 4.5 5.1 4.3 2.7 4.0 4.3 3.6 4.3
##  [613] 5.3 5.8 4.1 3.7 4.7 4.8 3.7 5.1 3.4 4.4 5.5 5.1 5.1 4.5 4.9 4.9 5.7 4.8
##  [631] 5.2 4.3 4.8 3.5 4.7 6.1 4.5 4.3 5.1 3.0 4.0 3.9 4.4 4.3 3.8 4.3 3.7 4.8
##  [649] 4.8 5.2 5.4 5.6 5.1 4.7 4.8 3.1 5.3 4.2 5.6 6.0 2.9 3.8 5.4 3.9 3.9 5.5
##  [667] 4.8 4.2 2.3 4.0 5.4 4.3 3.5 4.8 4.9 4.5 3.6 3.3 5.2 6.1 4.2 3.3 4.2 4.4
##  [685] 5.0 3.6 3.9 5.3 6.3 4.2 5.5 5.3 4.2 5.6 3.7 3.8 5.2 3.0 3.5 4.7 4.5 4.2
##  [703] 6.0 4.6 2.9 4.6 4.5 3.7 5.1 3.4 5.1 4.3 5.3 4.5 3.3 5.0 4.0 3.9 3.6 4.3
##  [721] 6.2 3.8 3.5 3.2 5.1 4.1 3.6 4.6 4.1 4.9 4.7 5.1 5.2 4.8 4.5 5.8 4.9 3.8
##  [739] 3.5 4.0 3.9 3.6 4.8 5.7 4.8 5.2 5.3 3.0 3.6 4.0 5.4 5.1 5.0 4.2 3.1 5.5
##  [757] 4.7 3.9 5.1 4.2 4.1 5.0 3.8 2.5 7.6 4.5 5.9 4.8 4.1 3.6 4.4 3.9 4.8 4.1
##  [775] 3.7 3.8 5.3 4.1 4.0 3.7 4.1 3.5 3.8 6.1 4.6 4.7 5.9 5.9 3.5 4.1 5.2 4.3
##  [793] 4.2 4.2 5.4 4.2 4.3 5.1 4.3 2.5 5.9 4.1 4.6 3.1 5.3 3.5 3.9 3.9 4.3 4.6
##  [811] 5.1 4.5 3.2 5.8 4.9 4.6 3.5 4.8 4.8 4.3 5.8 3.5 4.9 3.9 4.2 5.2 4.1 3.9
##  [829] 5.6 3.6 4.7 5.2 4.3 5.5 3.1 3.2 3.7 5.8 3.8 4.0 4.0 4.8 4.2 5.6 4.8 4.7
##  [847] 4.8 6.4 3.8 5.0 3.6 2.8 4.7 4.2 5.5 5.3 4.2 5.2 4.7 5.5 6.0 5.6 3.3 2.1
##  [865] 4.7 4.9 4.0 4.2 6.3 5.1 3.9 4.3 5.4 3.3 4.6 3.0 5.3 4.2 3.2 3.7 4.7 4.4
##  [883] 3.3 4.5 3.1 4.0 3.4 4.0 4.2 4.3 4.9 3.9 3.6 4.3 3.9 5.3 5.3 2.9 3.5 3.9
##  [901] 4.0 3.8 5.4 5.1 3.0 2.6 3.6 4.1 5.0 4.6 5.2 3.6 4.9 4.6 4.2 3.1 5.0 5.3
##  [919] 4.7 4.1 5.6 5.1 4.6 4.9 4.5 4.5 5.4 3.8 4.3 4.6 4.5 3.9 4.1 3.9 4.5 4.7
##  [937] 4.5 5.1 3.8 5.0 4.5 5.1 4.0 4.6 4.3 5.0 4.5 4.5 3.4 6.2 4.1 4.8 3.3 4.6
##  [955] 6.5 4.7 4.6 3.7 5.0 3.4 5.4 3.5 4.5 4.5 5.4 5.2 3.3 5.2 4.1 6.0 4.9 4.1
##  [973] 5.5 2.7 4.5 5.4 4.4 4.9 4.2 4.9 5.2 4.0 4.8 3.4 4.9 4.8 5.3 4.6 4.6 4.4
##  [991] 4.7 4.5 5.4 5.3 3.0 2.0 5.5 5.1 5.6 4.3
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
##    [1] 4.7 3.3 3.3 4.4 4.7 4.2 6.1 5.4 4.8 4.8 3.0 7.2 3.8 4.7 5.4 4.8 3.4 4.3
##   [19] 4.8 4.1 5.1 5.1 4.7 5.2 4.4 3.8 4.3 5.3 4.0 5.1 3.1 3.6 4.8 5.8 5.3 3.9
##   [37] 4.8 4.6 4.1 3.4 3.8 5.5 4.0 5.6 4.0 3.2 2.9 5.8 5.9 6.4 4.0 5.3 5.5 4.1
##   [55] 4.5 3.5 4.0 3.9 4.4 5.0 3.4 5.2 4.1 4.2 4.6 4.9 5.3 3.2 4.1 3.5 3.8 3.9
##   [73] 4.8 4.1 4.5 5.4 4.6 4.6 3.1 6.8 5.6 3.9 3.5 4.3 4.3 3.6 5.0 3.9 6.1 4.0
##   [91] 4.5 2.2 4.9 4.2 4.9 4.8 6.6 5.0 4.3 3.7 3.8 6.3 5.5 4.2 3.9 3.6 4.2 3.6
##  [109] 4.2 5.0 3.8 5.3 3.6 3.9 4.9 4.6 3.7 4.1 3.1 3.6 3.6 4.3 6.1 4.2 4.4 4.2
##  [127] 3.8 5.3 4.9 3.9 6.3 5.3 3.2 3.7 5.5 4.7 3.3 5.0 5.4 4.0 3.2 5.2 4.9 6.3
##  [145] 5.6 4.7 5.5 3.1 3.7 4.0 3.5 3.2 4.8 3.6 4.7 3.9 3.1 4.4 4.9 2.9 5.7 4.8
##  [163] 3.5 3.8 4.5 4.4 5.0 3.6 4.3 5.5 6.1 4.0 4.5 5.4 3.4 4.7 3.8 5.7 5.1 4.7
##  [181] 3.8 3.5 5.3 3.7 3.3 5.5 3.2 3.6 4.2 4.5 2.7 3.5 4.0 4.1 5.9 4.6 6.0 2.9
##  [199] 5.7 4.5 4.3 4.2 5.0 4.3 5.6 4.2 4.6 5.2 3.9 5.9 4.2 5.3 4.2 5.4 6.1 4.2
##  [217] 4.7 5.6 3.1 4.4 4.3 3.2 3.6 2.5 4.9 4.8 5.0 4.8 4.7 4.6 3.6 4.1 3.3 5.3
##  [235] 4.5 2.2 3.0 4.2 3.6 3.0 6.0 4.2 4.8 5.3 4.6 6.4 4.5 4.0 3.6 3.5 4.0 5.5
##  [253] 4.1 3.9 5.1 4.2 4.7 3.6 4.3 5.8 6.3 3.5 4.8 4.4 5.0 4.9 4.6 5.8 5.7 4.1
##  [271] 5.7 4.8 4.2 4.5 4.8 5.8 3.7 5.2 3.6 2.3 3.8 4.0 3.9 4.1 4.3 5.3 5.3 5.2
##  [289] 6.2 6.0 5.7 2.7 5.1 4.1 4.4 5.4 4.9 3.8 3.3 3.8 5.4 5.0 3.4 5.0 5.9 5.2
##  [307] 5.7 5.2 4.5 5.0 6.1 3.8 5.5 5.6 3.4 3.7 4.8 4.2 4.5 4.0 4.6 5.6 4.6 4.7
##  [325] 5.0 4.1 3.6 3.7 5.5 4.6 4.6 5.0 5.3 3.6 5.1 3.2 2.9 3.4 4.1 4.9 3.4 3.1
##  [343] 3.6 6.4 3.4 4.6 5.3 4.9 5.3 5.8 4.4 3.9 4.0 4.4 6.1 4.7 5.5 4.1 2.8 3.8
##  [361] 3.1 4.0 3.3 6.1 5.3 4.5 4.5 5.5 4.2 6.1 3.3 3.4 5.0 5.4 3.7 4.8 3.8 5.0
##  [379] 5.9 5.3 4.1 5.7 3.7 4.8 5.4 5.8 2.9 4.0 3.6 5.1 6.4 4.0 5.3 6.5 4.5 4.7
##  [397] 2.7 3.4 3.2 3.8 4.3 4.0 4.2 3.9 2.9 5.9 5.3 4.3 5.7 4.5 5.1 4.0 6.5 4.8
##  [415] 4.0 5.3 3.2 4.5 4.4 4.3 4.8 4.6 4.0 3.2 5.6 4.8 3.4 4.7 5.6 5.0 4.8 3.7
##  [433] 3.9 3.5 4.6 3.7 4.3 4.7 4.1 4.3 4.0 5.7 3.2 3.9 4.1 4.4 5.6 3.9 4.0 3.7
##  [451] 4.8 4.2 4.8 3.8 3.4 5.1 6.2 4.7 4.5 3.9 4.0 5.5 4.1 3.8 2.7 4.4 5.7 3.4
##  [469] 4.9 4.4 4.8 4.9 4.8 4.8 4.1 5.1 4.8 3.1 5.6 4.7 4.2 5.8 5.3 4.0 4.9 3.3
##  [487] 5.1 5.4 3.2 4.6 4.2 4.7 5.4 3.5 6.3 5.8 3.9 5.5 5.8 6.9 3.8 5.3 3.4 5.7
##  [505] 4.0 5.6 5.3 3.1 3.9 5.9 5.0 4.0 4.1 5.1 5.1 3.4 5.0 4.9 4.6 4.2 3.4 3.5
##  [523] 3.5 6.1 5.3 5.2 4.5 5.4 5.0 4.4 5.2 5.0 3.7 3.0 4.8 2.8 2.9 3.6 5.9 4.8
##  [541] 5.9 5.6 5.3 3.8 4.4 4.3 4.2 3.3 4.5 5.9 6.3 5.0 5.4 4.3 4.3 2.5 3.7 3.7
##  [559] 6.1 4.2 3.8 3.3 5.7 3.8 5.0 3.8 4.0 5.3 5.4 4.1 4.8 4.0 4.3 5.8 6.2 3.9
##  [577] 3.2 3.8 5.4 4.2 4.9 4.2 1.7 5.9 4.1 4.4 4.0 4.6 5.1 3.3 5.6 3.8 4.4 4.9
##  [595] 5.2 3.1 5.5 4.5 2.4 4.7 4.3 3.5 4.3 5.4 4.7 5.4 4.1 3.3 4.7 5.6 4.7 6.1
##  [613] 4.4 4.7 4.1 5.9 3.9 6.3 6.0 4.2 3.0 5.1 3.4 4.1 3.9 4.7 4.1 5.2 6.1 4.8
##  [631] 3.9 5.0 3.1 5.1 3.8 4.7 5.7 5.3 4.4 3.7 4.1 4.9 4.2 4.3 3.8 3.2 4.6 4.9
##  [649] 4.6 4.7 4.9 4.9 2.4 4.5 4.3 4.8 4.6 3.6 3.3 4.7 5.2 4.8 3.9 3.6 3.0 4.6
##  [667] 4.9 4.1 4.0 3.2 5.2 5.6 5.0 4.5 4.9 3.9 5.3 3.8 5.5 5.6 3.3 4.8 4.6 3.9
##  [685] 5.7 3.9 4.2 3.3 5.5 4.3 4.7 5.8 4.4 5.4 3.9 3.6 3.0 5.2 3.0 5.3 6.0 5.0
##  [703] 4.9 5.3 4.1 4.9 3.1 5.1 3.8 3.8 4.1 4.6 5.4 5.3 6.0 5.3 3.5 4.8 4.6 5.6
##  [721] 3.2 5.2 3.5 5.0 3.8 6.5 3.3 4.2 3.0 6.3 3.7 2.8 4.8 3.7 4.1 5.5 6.1 5.4
##  [739] 4.4 3.8 5.5 2.6 4.3 4.6 4.1 3.3 3.3 4.2 3.7 4.2 4.6 3.7 2.9 4.4 4.9 6.0
##  [757] 4.7 4.1 4.4 4.8 6.0 4.9 5.0 4.4 5.1 3.6 6.4 1.8 4.2 3.4 5.1 3.2 3.0 5.3
##  [775] 4.7 5.3 5.6 4.8 5.2 6.4 3.4 4.7 5.8 3.9 3.2 5.3 4.3 4.6 4.8 3.4 3.5 4.5
##  [793] 4.3 3.7 2.9 4.9 5.5 3.6 5.0 3.3 4.2 3.3 4.5 5.1 4.1 4.6 3.3 5.2 3.4 4.3
##  [811] 5.0 3.9 4.6 3.9 4.9 3.8 4.7 3.8 5.8 5.6 3.2 3.3 3.0 5.8 4.3 5.1 3.1 3.9
##  [829] 3.9 4.3 5.8 4.8 7.5 6.6 4.2 3.8 4.5 2.8 3.8 4.2 4.0 4.3 4.8 5.1 5.0 4.1
##  [847] 5.1 4.7 3.5 4.0 5.2 4.5 5.5 4.8 5.1 4.8 3.6 4.2 4.4 3.0 4.4 6.2 7.5 4.5
##  [865] 4.8 4.5 4.6 5.2 3.5 4.1 6.1 4.2 4.8 4.5 4.4 2.9 4.4 5.0 5.1 5.2 4.7 3.9
##  [883] 5.1 4.8 4.5 5.1 4.7 4.6 4.8 4.4 3.8 3.9 3.7 6.0 3.8 4.2 3.8 4.6 4.3 4.9
##  [901] 4.9 3.4 4.7 3.5 5.8 4.4 2.8 4.6 3.5 4.4 4.0 3.5 3.8 4.3 4.0 4.4 4.0 5.2
##  [919] 4.8 4.4 4.5 3.6 6.0 5.7 4.0 3.8 4.9 4.6 4.2 5.6 4.3 5.0 4.9 4.9 3.4 4.2
##  [937] 3.3 4.1 4.3 5.3 4.8 5.6 6.4 4.9 5.3 4.5 6.4 4.5 2.9 3.8 3.5 4.8 3.7 5.1
##  [955] 4.9 3.4 4.4 4.3 5.4 4.3 5.1 4.6 5.3 4.0 3.5 5.8 4.3 4.1 5.1 3.2 4.3 4.5
##  [973] 5.7 5.6 3.7 4.5 3.4 6.0 3.8 2.3 3.9 1.5 6.7 3.8 7.4 3.9 5.2 5.0 5.7 3.5
##  [991] 3.9 4.8 5.0 4.6 5.5 5.3 5.2 4.8 6.8 6.0
## 
## $func.thetastar
## [1] -0.0031
## 
## $jack.boot.val
##  [1]  0.45266106  0.36235632  0.33242424  0.11823204  0.08132530 -0.05367847
##  [7] -0.19376855 -0.24040698 -0.40572289 -0.48804348
## 
## $jack.boot.se
## [1] 0.933479
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
##    [1] 3.8 3.7 4.3 3.8 5.3 3.7 3.6 5.5 5.1 4.4 4.6 4.9 3.2 2.9 4.2 4.2 6.4 2.9
##   [19] 3.2 3.0 4.9 3.7 5.2 4.8 4.8 2.5 5.1 3.6 4.6 3.1 4.6 4.9 5.1 3.4 4.6 3.8
##   [37] 3.4 3.8 4.1 5.0 2.9 4.0 5.5 6.9 4.4 5.8 2.9 5.4 4.4 5.5 5.6 4.8 4.3 4.3
##   [55] 3.8 3.3 4.6 7.0 3.6 6.2 5.2 3.8 4.6 3.8 5.1 6.1 5.5 3.6 4.0 5.3 4.1 5.2
##   [73] 4.6 4.4 2.9 4.0 4.7 4.3 5.0 5.7 5.1 4.2 4.2 5.1 5.4 4.9 5.6 4.2 3.5 7.0
##   [91] 3.6 5.0 4.8 5.2 4.8 4.5 5.2 3.5 4.2 3.3 4.2 4.4 3.6 4.3 6.3 3.5 5.1 3.6
##  [109] 4.2 4.5 4.7 4.4 2.5 4.0 5.2 6.6 5.7 4.5 3.3 4.0 3.7 3.6 3.9 5.4 3.1 4.6
##  [127] 5.3 4.7 3.0 4.6 5.2 3.5 5.8 4.1 4.2 4.7 4.7 4.7 4.5 5.1 4.4 4.1 4.5 6.2
##  [145] 4.1 5.2 4.9 4.0 1.8 5.6 2.9 4.6 4.7 3.9 5.5 5.7 4.9 4.3 4.6 2.6 3.7 4.7
##  [163] 4.5 5.0 3.9 5.9 2.7 4.0 4.5 4.9 3.9 4.8 5.4 6.0 4.1 3.2 4.0 4.2 4.2 4.2
##  [181] 6.3 4.8 5.0 4.8 2.9 5.2 2.5 4.6 4.2 4.9 5.2 4.4 5.4 3.9 5.9 4.9 5.9 3.6
##  [199] 5.4 2.1 4.9 5.9 3.6 3.7 3.1 3.5 3.9 4.1 5.3 4.7 3.4 6.5 4.8 4.4 4.4 5.1
##  [217] 5.8 3.9 4.6 4.1 3.0 3.7 4.9 4.2 4.2 3.2 5.5 5.6 3.4 4.3 3.2 4.5 4.3 4.8
##  [235] 4.3 3.6 3.8 5.5 5.6 4.9 5.5 6.4 2.7 4.2 4.7 4.5 3.0 3.9 3.9 2.7 3.4 4.8
##  [253] 3.3 6.9 5.8 4.8 4.0 2.3 5.4 4.5 4.4 4.0 4.5 5.0 5.0 4.4 5.4 2.8 4.6 4.7
##  [271] 5.0 5.3 4.0 4.8 5.3 3.8 3.4 4.3 4.3 4.7 4.4 2.7 4.5 5.2 4.2 4.7 3.2 3.7
##  [289] 4.2 4.4 4.8 3.9 4.6 6.0 5.2 4.2 3.6 4.7 4.7 5.7 4.6 5.7 5.1 2.4 3.7 3.9
##  [307] 5.1 4.4 6.1 5.3 3.5 4.1 4.5 3.5 4.3 3.7 4.7 3.8 5.0 3.4 5.7 5.6 4.2 3.9
##  [325] 5.3 5.8 4.9 5.3 4.8 5.0 5.3 5.1 4.9 3.5 3.1 3.6 6.3 4.5 4.8 4.8 6.2 4.1
##  [343] 2.5 5.4 4.8 4.8 3.3 5.0 5.1 6.0 5.0 5.4 5.5 4.2 3.6 5.0 4.3 2.4 4.4 4.0
##  [361] 2.0 4.8 5.7 3.4 6.1 5.0 4.1 6.0 4.5 4.6 5.2 3.5 5.5 5.2 3.5 4.1 2.8 3.2
##  [379] 5.5 4.5 5.6 4.2 5.0 4.4 4.4 5.1 3.7 4.6 5.5 5.4 4.4 3.0 4.5 6.1 4.8 4.9
##  [397] 3.6 5.3 5.0 3.2 5.0 3.0 6.0 3.5 4.4 5.3 6.0 4.9 5.4 4.1 4.7 4.5 4.6 5.2
##  [415] 2.3 4.8 4.8 4.9 2.9 4.2 5.5 3.9 4.4 3.1 5.0 4.3 4.2 4.0 3.3 5.7 5.4 4.0
##  [433] 5.0 5.0 4.3 4.0 4.7 3.4 2.8 4.3 5.1 5.4 4.0 3.6 4.0 5.7 4.9 5.1 5.2 2.9
##  [451] 4.8 3.5 4.7 4.1 4.5 5.7 6.3 5.6 5.6 5.4 5.0 4.5 4.7 5.4 3.5 4.2 4.2 4.9
##  [469] 3.9 3.7 5.3 3.7 6.9 4.4 3.9 4.0 4.0 4.4 3.3 5.5 3.4 4.6 2.8 5.0 5.0 5.5
##  [487] 3.0 3.8 4.0 6.8 4.6 4.7 4.6 5.0 2.9 5.7 5.6 3.0 4.8 4.6 4.9 4.9 3.7 4.8
##  [505] 4.6 5.4 4.2 3.6 4.2 5.0 4.2 6.0 4.5 5.6 3.7 4.3 5.6 5.1 3.1 3.4 5.3 3.0
##  [523] 4.0 4.7 5.4 5.2 5.0 5.1 5.9 5.0 4.5 3.9 4.9 4.8 4.7 5.1 3.0 4.0 4.4 3.4
##  [541] 3.4 5.4 4.7 4.8 3.5 4.2 4.2 3.7 4.0 3.7 5.0 3.5 5.3 4.0 3.6 6.1 4.0 3.8
##  [559] 4.7 3.5 5.0 4.3 4.8 5.7 5.0 2.1 5.2 5.7 3.7 5.2 5.5 4.5 5.8 3.8 3.7 5.1
##  [577] 5.1 4.8 6.4 5.3 4.5 5.2 4.1 3.5 4.4 5.1 4.0 4.9 4.8 5.4 3.5 5.4 6.0 5.5
##  [595] 4.9 4.8 4.3 4.3 5.1 5.0 3.4 5.9 5.2 5.3 3.5 4.3 3.2 4.8 4.5 4.1 4.9 3.3
##  [613] 4.9 3.6 4.6 2.5 4.2 4.5 2.8 4.6 3.6 4.9 3.4 4.1 4.9 4.2 4.7 5.2 5.1 6.0
##  [631] 4.3 3.5 3.5 5.2 4.8 4.8 3.4 5.2 4.5 4.9 2.7 4.6 2.8 4.1 5.3 4.9 4.5 5.0
##  [649] 5.4 4.3 4.1 4.7 4.0 5.3 6.5 4.8 4.7 5.0 5.6 5.1 5.1 4.9 4.3 4.0 3.7 6.0
##  [667] 3.5 4.3 4.8 5.0 5.5 2.9 5.5 4.3 4.5 3.6 4.8 4.5 3.5 6.3 4.8 4.8 5.2 5.8
##  [685] 4.2 6.1 3.0 4.9 3.7 4.8 4.3 4.3 6.3 4.5 6.1 5.2 4.3 2.9 5.8 7.3 5.1 5.5
##  [703] 4.1 4.2 4.2 6.3 4.7 4.6 3.8 5.1 3.7 3.5 5.3 3.6 4.8 5.2 4.5 4.3 4.6 5.3
##  [721] 3.1 5.7 5.7 3.4 5.6 4.1 4.7 5.9 5.4 6.2 4.9 4.0 5.7 3.6 3.6 4.0 4.6 3.8
##  [739] 6.1 4.8 3.6 4.4 4.4 4.1 6.5 5.6 3.9 5.5 5.7 5.8 4.6 2.8 3.8 2.5 4.7 5.2
##  [757] 5.1 4.1 4.8 5.5 4.4 5.5 3.4 4.4 5.7 6.2 3.2 4.6 4.7 2.4 5.1 5.0 5.2 5.7
##  [775] 4.3 6.2 3.8 4.3 4.7 4.9 3.6 4.0 3.1 5.3 5.8 3.6 5.2 4.7 2.5 5.0 5.5 3.9
##  [793] 2.6 3.7 4.8 4.5 5.6 5.2 3.9 4.9 4.2 3.3 4.4 5.2 6.1 5.7 3.9 4.6 3.7 3.4
##  [811] 5.6 3.8 4.0 5.2 4.1 4.5 4.4 3.6 2.5 5.6 6.5 4.0 3.6 5.7 4.1 5.0 5.7 4.5
##  [829] 3.3 3.4 5.4 4.3 4.6 5.1 3.7 5.3 4.7 5.5 5.1 4.5 3.9 5.6 4.5 3.6 4.2 3.6
##  [847] 4.6 4.8 6.7 4.7 3.3 4.6 5.6 4.7 4.3 4.3 4.5 3.8 5.2 6.7 3.8 5.1 5.3 4.4
##  [865] 4.7 4.1 5.3 3.7 5.2 4.3 4.2 4.4 2.7 5.5 4.3 3.6 3.8 4.5 3.8 5.3 4.8 5.7
##  [883] 5.8 4.7 5.2 3.9 3.3 5.2 4.1 3.8 4.5 4.9 5.5 3.8 4.2 4.6 4.8 5.0 4.9 2.8
##  [901] 3.9 4.3 4.5 4.9 6.3 4.0 5.5 4.7 4.6 5.5 4.1 3.4 4.3 5.7 3.7 6.1 4.7 5.5
##  [919] 4.2 4.2 3.3 2.6 5.4 4.3 5.2 4.8 4.9 4.2 3.7 5.1 4.2 3.6 5.3 5.5 6.0 4.6
##  [937] 3.5 4.2 4.3 5.0 5.3 3.2 3.8 5.9 5.0 3.1 4.2 5.5 4.4 4.0 3.8 4.7 3.5 3.2
##  [955] 4.9 4.1 3.9 5.9 5.4 5.9 4.7 2.5 5.1 6.3 4.8 4.3 6.5 4.3 4.8 3.4 3.9 5.5
##  [973] 5.2 3.5 6.4 5.8 6.3 5.1 3.8 5.0 5.4 4.4 6.6 4.1 3.0 3.6 4.3 4.7 5.8 3.0
##  [991] 5.0 4.7 5.7 4.0 6.4 5.4 2.9 4.7 6.4 4.2
## 
## $func.thetastar
## 72% 
## 5.1 
## 
## $jack.boot.val
##  [1] 5.500 5.400 5.300 5.300 5.200 5.016 4.900 4.800 4.700 4.500
## 
## $jack.boot.se
## [1] 0.9402549
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
## [1] 1.233616
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
##      shape       rate   
##   2.1144918   3.8177633 
##  (0.8812392) (1.7946857)
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
## [1] 1.2118737 1.2994186 1.2053020 1.4735390 0.1508720 0.8204117
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
##    [1]  1.6522644880  0.1660145129  1.3802901151 -0.3801000466  0.7876642164
##    [6]  1.9553297601 -0.3195775480  1.4559984740  1.0296487938  1.6223743985
##   [11]  0.9021354567  1.0648961477 -0.3813353643  0.8619979805  0.8684948476
##   [16] -0.5843150345  1.7316698217  1.0612746239  0.7356582753 -0.8057765225
##   [21]  1.4935750394 -0.1364411771  1.3734892193  1.6272434181 -0.1068476059
##   [26]  0.7647242309  0.9477231441  0.2271762230  0.9609127989  0.8632308570
##   [31]  1.2892804667 -0.5807154598  0.9537927467 -0.1802961877  0.8171183176
##   [36]  0.5306183493  0.8296647466  1.1807370600  1.2517494753  1.5228220578
##   [41]  1.3560948026  0.1477593600 -0.0609929250  1.2533138736  0.7881194293
##   [46]  1.1213196724  1.0722399408  0.9147920502  1.1259733605  1.0802215398
##   [51]  0.4544161816  0.6420421904  0.5215812719  0.2815866989  1.3049757966
##   [56] -0.3087638420  1.2743763471  0.4627136472  0.9088188593  1.4584550505
##   [61]  0.8041563059  1.2340183170  0.8345665199 -0.1836289635 -0.0895798616
##   [66]  1.3984746436  0.1897664086  0.6999960479  1.0876970902  1.5778359254
##   [71]  0.1250976537 -0.1689333528 -0.4142012511  0.8137086001 -0.2896972764
##   [76]  0.1739663011  0.6222550483 -0.4678066076  0.9035045033  1.0119868160
##   [81]  0.8825986089 -0.2526403523 -0.0420208165  0.5272348741  1.0737749293
##   [86]  0.0685768227  1.1579442726  1.1577369287  0.6223109945 -0.3664784148
##   [91]  1.7954257122  0.7498188721 -0.3130429354  0.8974250742  0.1657728422
##   [96]  0.9908272598 -0.1769655493  0.5847647844  0.1386763925  0.8612828483
##  [101]  1.5057006630 -0.0589824046  0.5994983244  0.9020359478 -1.0270561231
##  [106]  0.4062566006  1.1054788499  0.9037038144  0.5780435545  1.9316177316
##  [111]  1.0446745102 -2.3512206156  1.0117359459 -0.5084252824  0.6451312227
##  [116]  1.6654188481  0.8447184118  1.1407769478  1.7619337885  1.3496073016
##  [121]  1.5932697135 -0.2742122351  0.2155044538  0.6472842071  0.0894816748
##  [126]  1.4671275077  0.1871444727  1.7051910284 -0.4219341091  0.7638395012
##  [131]  0.8156076344  1.7400996319  1.6492774707 -0.1253072709  0.2767215180
##  [136]  1.2434915012  1.5841614045 -0.4033662431  1.2778289007 -0.5013407409
##  [141]  0.3460907962  0.4786364688  1.8011853803  0.8786697927 -0.5453210296
##  [146]  1.4104195584 -1.0107712621  0.8019300922 -0.8003348129  1.4504149700
##  [151]  1.3803054497  1.5355493889  0.8520093247  1.3873040884  1.4449820626
##  [156] -0.3602516724  1.1197755386  1.3608946978  0.6982510654  0.5655231711
##  [161]  0.3538805908  0.6615597446 -2.3018706586  0.7984901203  0.4112237632
##  [166]  0.6762311544 -0.6419913873  0.7494865838  1.5973614146  0.7615822170
##  [171] -0.7506973621  1.3299872008  0.5600048920 -0.6671325350  0.9700809753
##  [176]  1.9378647479  1.1665100330 -0.6279997630 -0.2697654626  0.7422868704
##  [181] -0.9026098996  1.6728668020  0.6119862130 -0.4299293931  0.1776835496
##  [186]  0.7615117783  0.6423680977  0.4765177450 -0.2856413952  0.5724769432
##  [191]  0.5291757254  0.6335641686  0.9068264509  0.4037386430  2.0401903017
##  [196]  0.9557326799 -1.0303566539 -0.1124109032  0.7877426431  1.4638497056
##  [201]  1.0269267009  0.7253925903  0.6801440520  0.4256958704 -0.4004220413
##  [206]  0.5761452314  0.5921755933  0.9809648248  0.4784839114  1.0620311726
##  [211]  0.8498224976 -0.3279178646 -0.4408470052  0.4348960696  0.1072935122
##  [216]  0.1195450211 -0.0720202575  1.8591651790  0.7335828379  1.1191518370
##  [221] -0.2134127294 -0.5272094560  1.3303746239  0.7638461615  0.6225658747
##  [226] -0.0812377906  1.6089319094  1.1729744122  0.8095489309  0.2312919506
##  [231]  1.2070261768  1.1167412722  0.0827560090 -1.1229917250  1.3147895543
##  [236] -0.1800070476  1.2609737128 -1.4932817156  1.6257510375  0.4931644146
##  [241]  0.2163285671  1.8107025681  0.5114865283  0.8449990374  1.3406359858
##  [246]  1.5612246067 -0.9169040961  1.1339950498  1.2329981528  1.3822110703
##  [251] -1.1249864104 -0.6602904184  1.1996245459  1.7619982116  2.0676721154
##  [256]  0.9582437936  1.3211178692  1.0201728186 -1.1278072607  1.2583540678
##  [261] -0.4223921037  0.0008551674 -0.3984948708  0.0572413859  1.2136394726
##  [266] -0.2866574169  1.1513051878 -0.4261817482  0.5903966224  1.1862765790
##  [271]  0.9193231230  1.1563812739  1.3814638516  1.0679710897 -1.1573523408
##  [276]  1.2574166839  0.8739466406  1.1234173676  1.5272147382  1.0610506210
##  [281] -0.4156642692  2.0930715276  0.1075821108 -0.1584095428 -0.3520333331
##  [286] -0.4998374192 -0.3688174417  1.9057168377  1.0510706901  1.4408787780
##  [291] -0.1163057224  0.7241527588  1.1945637030  0.0371507675  1.0059533858
##  [296]  1.0874573325  2.1835009859 -1.0392436203  1.3430808378  0.0123884532
##  [301]  1.2115635870  2.0398135573  0.7753464643  0.7747278507 -0.0489004032
##  [306] -0.4368074805  1.6878786990  1.5786453196 -0.3715529017  1.7311139895
##  [311] -0.4452348354 -0.7314868652  1.2144312822  0.0705598843  1.4627601793
##  [316] -0.2081141308 -0.3625707773  1.1005578828  0.8593557319  0.8222961152
##  [321]  0.2222310709  0.8195523909  1.7309534337  0.4776998558  0.9375984883
##  [326]  1.4667419162  1.6392349600  0.0905951705  0.2261752043  0.3034328147
##  [331]  1.8312961074 -0.3214357696  1.4969396422  0.8678999798  1.0752376706
##  [336]  1.2689590417  0.8847149985  1.0955723192  0.8975121872  1.1256970728
##  [341]  0.0061235947  1.3791316893  0.6239079948  1.6170467465  1.4225832105
##  [346]  0.6155035332 -0.0821761564  1.3000402287 -0.6238378330  1.3400226223
##  [351]  1.1646320377  1.1590291672  0.9023520353  1.3507654348  1.1593677469
##  [356]  0.8836849374  1.1700123080 -0.5507786645 -0.6315337872  1.9698760512
##  [361] -0.7198785182  0.2984131402  1.2143525287  1.2122902787  1.2320636346
##  [366] -0.2847346671  0.5158224774  0.3653411159  0.9167347041  1.4611988642
##  [371]  1.4640130141  0.6334352884  0.5635105880  1.3873951738 -0.2819883796
##  [376]  0.9554111031  0.0436769736  0.2767765863 -0.2890299857  0.2615593916
##  [381] -0.4004220413 -0.5569831288  1.1004495974  0.9181515029  0.2982291361
##  [386]  0.9519360572  0.5981385634  1.0754673914  0.9829932692  1.5935640643
##  [391] -0.2111015040 -0.7193150103  1.9990518146  1.3200428678  1.5801417821
##  [396]  0.4710712143 -0.7749429636 -1.3330425016  0.9612559804 -0.4726108361
##  [401]  0.9871966488  1.1808909402  0.1385371209  0.9456679020  1.0301301488
##  [406]  1.4867263401  0.7628821983  1.4892913026  0.2725845584  1.4156909447
##  [411]  1.3343507865  1.0859832652  0.2241785862  0.7972518151  1.5634826154
##  [416]  0.1122384699  1.1823715286  0.0010722357  1.6320360840  1.0963436275
##  [421]  1.2473610898 -0.5634001022  1.1476552316 -1.1752657220  0.8090161325
##  [426]  1.4380823647 -0.1500572140  1.0197651629  1.2774435685 -0.7181657821
##  [431]  1.5595618868 -0.4313821880  1.6340398042  1.4861405364  0.8501265816
##  [436]  0.1401856909  1.8823695272  1.0606714276 -0.3351686711 -0.0874286086
##  [441]  1.0880125839  1.5118465196  0.7093191632 -0.1016001928  1.1042794557
##  [446]  1.2662613906  2.1241907611  0.9150505133  1.1757633947  0.4170756219
##  [451]  0.6957601686  0.6485689472  0.9347717708  0.8726699048  0.8089238913
##  [456]  2.2585132355 -0.5333518186  1.2014923995  0.7110011257  1.7594191421
##  [461]  1.4372028709  1.0254484094  0.2021142920  0.7921102567  0.2426707056
##  [466]  1.2455805066  0.9689244761  0.3442266664  0.3165396260  1.2545836881
##  [471]  1.2906611114 -0.5283765245 -0.6022938284 -0.6742015605 -0.9668596673
##  [476]  0.8365132724 -0.4237536636  1.3034714582  1.7720613407  0.6541952833
##  [481]  0.5804883138  1.7494799296  1.0393179890  1.4045110815  1.1689966048
##  [486]  1.0422392534  0.9970251178  1.3957905449 -0.0801530200  1.6201791714
##  [491]  0.8598524582 -1.0622128311  1.1159305534  0.0142467918  0.4370523361
##  [496]  0.1498747810  0.6858358308 -0.0116273074 -0.3642615301  1.3335117104
##  [501]  1.2422603341  1.6857739177  0.1444441223  0.9185378958 -0.0059824035
##  [506] -0.2805164453  0.9141781169  1.3009078200  0.3850769757  0.0105666371
##  [511]  1.2006645867  1.3000085304 -0.3206833296  1.3814001829  1.0727435490
##  [516]  2.1813901344  0.0017415993  1.6472414187 -0.6325366036 -0.4494325206
##  [521] -0.1718993275 -0.8903257191  0.1167034555  1.0057554941  1.2726786491
##  [526]  1.2482522592  1.2782809295  1.5844016178  0.9745239005 -0.3286902153
##  [531]  0.7768450949  0.0804116152  0.8109692192  0.9141162292  1.4448117558
##  [536]  1.7341585792 -0.2132829913 -0.2143499673  1.3447343578  1.4062852681
##  [541]  1.3076051380  0.0304068468  1.2250330206  1.2320636346 -0.1206010577
##  [546]  1.4585044351  2.0130374141  1.1350672029  1.2452248007  1.1455858328
##  [551]  1.9076020445  1.1574955868  0.9159636203  0.1214115647  0.9904151201
##  [556]  0.7950353785  1.1442140624  0.9889273374 -0.2813014218  0.2108636420
##  [561]  0.0943165072  1.6213493743  0.4218052730  0.7597449355  1.2684502895
##  [566] -0.2418210683 -0.4107965116  1.0479198236  1.9342552615  0.0126209288
##  [571]  0.1254568976 -0.0599349962  1.0222350717  0.0461877059  0.0370784135
##  [576] -0.4521932770  1.0497888216  0.5018253421  0.9974804109  1.7333635292
##  [581]  1.0988759254  1.0031363855 -0.8059804131 -0.5623221854  0.4162995529
##  [586]  0.7067778610  1.0644538006  1.4874534353 -0.2057285820  1.9730737742
##  [591] -0.4156376900  1.7253094541 -0.3751605421 -0.0116622388  0.3936929933
##  [596] -0.2115603692  1.3558723380  1.2312609822 -1.1345675518  0.2982291361
##  [601]  1.6681025696  1.2384283318 -0.1243538196  1.5593188576  1.0500960214
##  [606]  0.5184894007 -0.2252696113  1.7490411373  1.3111827286  2.1098907725
##  [611] -0.2871637386  0.5028289427  1.2378032174  0.8356143527  1.4795680220
##  [616] -1.1180562669  0.6026105283 -0.4864012727  0.3731941489  0.7044884377
##  [621]  1.3346511275  1.2700997585 -0.4010120744 -0.0193699352 -0.4078460373
##  [626]  0.9081031094  1.0119435660  0.9659589845  1.7198827488  0.4006729478
##  [631] -0.5471676892 -0.2017399659  1.3411485461  1.0369017166  0.9266755515
##  [636]  1.0769030819 -0.5099518771  1.4418297925  1.5510962050  0.8026256738
##  [641] -0.8123317171  0.6626028358 -0.1867489297  1.1121345578 -0.0447827085
##  [646] -1.0549449269  0.1530659223 -0.7338796541 -0.5075162432  0.3495728876
##  [651] -0.1475360026  0.7751784351  1.2972964316  0.7884264410  0.8422276216
##  [656]  0.2827111522  1.4250594845 -0.0327098612  0.6813931607  0.9091835531
##  [661]  1.0591967124  0.0171209748 -2.3607912463  1.1065212459 -0.2652714150
##  [666]  1.5464397612  0.9686135001  1.2658583679 -0.0347446393  0.1933334614
##  [671]  0.8594233182  1.6995409933  0.9979840813  0.6776422870 -0.4569683601
##  [676] -0.1267494054  1.2383031215 -0.5115201713  0.9314219771  0.8952429119
##  [681]  1.1447494424 -0.0517364299  0.9542955201 -0.2328397535  1.2939129118
##  [686] -0.9158495519 -0.4535159215  1.4776500815  0.9414206184 -0.0736773526
##  [691]  0.6233334375  0.2507900017  1.3798964073 -0.8042400069  1.1621614331
##  [696]  1.4259075589 -0.0527494571  1.4879569852 -0.4627850407  0.9353799483
##  [701]  1.3076886428  1.7489621289  1.6770481448  1.5807728671  1.2788587709
##  [706]  1.3632690683 -0.2005040621  1.2615134189  1.6654519182 -0.2673897888
##  [711] -1.0939676595  1.1613266795 -0.1741314673  1.2173561170  1.2231112602
##  [716]  0.5303279548  0.0301806557  0.9974677409  1.0482336806  0.3168610123
##  [721]  1.6822469304  0.9241927186 -0.7417867240 -0.0128921369  1.0244085963
##  [726]  0.8683972403  0.7136831411  0.6866914580  1.3350734917  1.1104802486
##  [731]  0.9180215298  0.7809178075  1.5385079928  1.4548277775  1.5261017515
##  [736]  0.9962186622  0.4536660969  0.6373142524  1.6674703096 -0.1435501898
##  [741]  0.6956075130  0.1371352425  1.0186782645  0.2933774225 -0.4146189660
##  [746]  1.0146712238  1.2886745814  1.4781754440 -0.3879055543  1.7319780485
##  [751] -0.1286911544  0.9996873472  0.4755461195 -0.8138706490  1.0862999514
##  [756]  0.6581050113  1.8571157805  1.2601838685  0.9219736309  1.4144961155
##  [761]  0.3529164422  0.8829424418  1.4671469504  1.3698589849  1.6376554348
##  [766]  1.8664223337 -0.1011300952  0.7933451275  0.5658107059  0.8251857938
##  [771]  1.7802467305  0.8413874531  1.5982303971  0.9578965900  1.3048320182
##  [776]  0.1642807700  0.7705524032  0.7897386846 -0.7535576091 -0.7586059838
##  [781] -0.6868185246  1.0700405777  1.3100692929 -0.1887320558  0.6389029353
##  [786] -0.7122936178  0.0849243148  1.3707956913  0.7559331595 -1.2164858089
##  [791]  0.7670964359  1.2029047074  0.9171380367 -0.5315342770 -0.0433782677
##  [796]  1.2156765966 -0.7782174592  1.7445455642 -0.4997777058  0.7841228078
##  [801]  1.6063608132  1.2454895352  0.0007021181  1.6205289909  0.8514111958
##  [806]  1.6394574707  0.6876675263  1.6878807060 -0.2900654977  0.0918886242
##  [811]  0.5461347386  1.2496133802 -0.0552302699  0.3079812105 -1.2428302190
##  [816]  0.9761483706  1.2668527132  1.7134568016  0.7884264410  1.0963538031
##  [821]  1.6770979948  1.5633625213  0.7249492710  0.8790333986  0.0184973107
##  [826]  1.4845841707 -0.3688174417 -0.2390314303 -0.4651892563  1.1045771827
##  [831]  1.3176004925  0.4049671458  0.8619979805  0.9988832270  0.9114044090
##  [836]  1.6476161503  1.5915256762  1.4361543924  0.8650229352  0.0123884532
##  [841] -0.2675695448 -0.1584648543  1.1145622767  1.3618137112  0.1355262162
##  [846] -1.5012101486  1.0686372365 -0.3767847782  1.0269267009  0.9187505425
##  [851]  1.0596881336 -0.9011214941  1.7360703162  0.4651703128  1.5011525312
##  [856] -0.4595494025  0.7207777292  1.6218576246  1.2539107657  1.7979502754
##  [861]  0.5299665183  0.9322726034 -0.0805315274  1.2271928699  1.3445389214
##  [866]  1.7492683664 -0.0810544703  0.4491616558 -0.0546327228  0.8772393993
##  [871] -0.0940841009  1.7231660987 -0.8502123106  0.5655440991  0.3427227075
##  [876]  0.5824106827  1.4664954155  0.6758199494  1.1559132151  1.1167412722
##  [881] -0.3777378520  1.2636476040  1.3492893964  0.9009707359 -0.3546149616
##  [886]  0.7818223449  1.2444672594 -0.0020346136  0.0808418014  2.0415538074
##  [891]  0.4186365363  0.7932673899  1.3634225549  0.7002719218 -1.0918042566
##  [896] -0.7899141508  0.6117711844 -0.8768294582  1.3834877368 -0.1227262918
##  [901]  0.0417806568  1.1926435225 -1.0107712621 -0.1996453791 -0.5343612530
##  [906] -1.1121437032  1.0398131032  0.8352495153  1.0259138184  0.2536529013
##  [911]  0.8345047629 -0.3684381795 -0.1686293222  0.6834856666  2.0108096197
##  [916]  0.8714621321 -0.0608197352  1.7198659163  1.5666552970 -0.5885434824
##  [921]  1.1851699753 -0.0344828135  0.8334866596  1.4520464806  1.2144206491
##  [926]  1.7211156029  1.0712833429  1.3753603637  1.1603399877  1.7960876111
##  [931]  1.7756104509  1.3587666914  0.8432652430  0.9041724932  0.3617130059
##  [936]  1.0103134208  1.0207683607 -0.2122357558  0.1758231418  0.9715623809
##  [941]  1.2484079290  0.9180215298  1.1769910078  0.5393976606  1.5315802668
##  [946]  1.1404888621  0.7985529163  1.4479636341 -0.1800070476  1.0884138157
##  [951]  0.9160497075  0.5212516480 -0.1893195250  0.2805520282  1.8685452939
##  [956]  1.6913908908  0.5957606845 -0.4502897657  1.2141443016 -0.2902500194
##  [961]  0.8406844439 -0.4241961444  1.0140129615  0.8296647466  1.9838086442
##  [966] -0.2902500194  1.7989227107  1.0336848454  0.4175514227  0.1960171101
##  [971]  1.0563900500  0.9732870631  1.6563546553  0.9929011579  1.0949009115
##  [976]  0.2619753436  1.3356249457  1.7072318139  1.1845587432 -0.4223649795
##  [981]  0.8846982953  0.6997330954 -0.5914970273 -0.5189080901  1.4845841707
##  [986] -0.4911681361 -0.5163013087  1.3071561932 -1.1167923931  0.8590514348
##  [991] -0.4290326008  0.7526002639  1.0510706901  0.2508572564  0.7568185648
##  [996]  0.4753130231  1.1511541059  0.9487177579 -0.6823651365  0.9115301974
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
##   0.55385677   0.37567176 
##  (0.11879784) (0.08399945)
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
## [1]  0.19992249  0.61116183 -0.05259683  0.48108786 -0.03772860  0.22739899
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
## [1] 0.0011
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9194096
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
## t1*      4.5 -0.04794795   0.8894017
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 6 8 
## 2 2 1 1 3 1
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
## [1] -0.0278
```

```r
se.boot
```

```
## [1] 0.8791473
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

