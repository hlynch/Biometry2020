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
## 0 2 3 4 5 6 9 
## 2 1 2 1 2 1 1
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
## [1] -0.0093
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
## [1] 2.675714
```

```r
UL.boot
```

```
## [1] 6.305686
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
##    [1] 3.6 5.7 5.0 6.1 3.6 4.5 3.2 4.0 2.8 5.6 4.2 5.2 4.7 5.3 4.1 4.0 5.1 4.2
##   [19] 5.4 5.5 3.9 4.6 4.8 3.7 6.2 6.3 3.2 5.2 4.9 4.0 4.2 5.6 4.8 4.9 3.9 3.5
##   [37] 5.7 5.3 4.7 3.7 3.1 3.8 4.6 5.5 4.5 4.1 4.6 4.8 4.1 4.9 4.0 2.7 5.1 2.5
##   [55] 4.3 3.7 5.0 3.7 4.3 4.3 3.6 4.5 4.7 4.4 4.7 4.3 5.9 3.1 5.3 6.5 3.9 3.7
##   [73] 5.3 5.1 5.1 3.8 5.4 4.6 4.1 5.3 5.1 6.0 6.3 3.9 4.0 3.9 3.3 4.6 5.9 3.1
##   [91] 5.3 5.2 5.8 3.6 4.4 3.9 4.5 3.6 4.8 5.6 3.4 5.5 4.1 3.8 3.5 5.8 2.6 1.9
##  [109] 4.6 4.8 5.5 6.3 4.7 4.0 3.9 3.1 7.7 5.0 4.0 3.1 3.6 4.4 3.9 3.9 4.9 4.4
##  [127] 4.2 3.8 3.5 4.5 4.4 4.2 5.0 3.4 5.6 4.8 3.7 3.9 4.5 3.3 4.1 3.0 3.6 4.5
##  [145] 4.3 3.9 4.5 4.7 5.9 4.1 5.7 4.6 4.5 4.0 3.8 4.3 3.5 5.1 4.9 3.7 5.1 4.1
##  [163] 5.2 5.0 6.5 3.9 5.2 3.9 5.4 4.9 4.8 5.2 6.1 5.1 5.2 6.0 5.1 4.8 5.2 4.4
##  [181] 3.0 5.5 3.2 4.2 5.7 6.2 4.6 6.0 2.8 4.8 3.2 6.3 5.2 5.0 3.4 4.9 4.3 3.3
##  [199] 5.1 3.9 3.7 4.6 4.6 3.6 3.6 2.2 4.6 4.2 4.1 4.2 6.4 5.2 4.8 5.0 5.5 4.2
##  [217] 4.1 4.5 5.8 4.3 4.5 5.8 5.0 4.5 5.4 5.2 5.7 4.5 4.6 5.7 4.7 4.3 5.6 2.7
##  [235] 5.7 4.3 5.0 3.1 5.2 4.7 5.8 4.9 4.4 5.9 5.2 3.4 5.0 5.5 6.1 4.6 4.0 4.7
##  [253] 4.4 3.6 4.9 5.2 3.6 6.3 5.0 6.0 5.1 5.1 5.5 4.7 3.2 5.0 2.5 5.1 3.6 5.5
##  [271] 5.1 2.9 3.6 5.1 4.0 3.8 4.0 4.1 5.0 4.8 5.1 5.6 4.3 5.2 3.4 5.4 5.6 4.1
##  [289] 3.2 5.4 3.1 2.6 4.3 4.1 4.5 4.6 4.0 5.5 4.0 3.8 4.1 4.0 4.5 5.3 5.1 5.4
##  [307] 4.2 4.8 4.9 4.2 3.7 5.1 5.1 4.2 2.8 4.6 5.4 4.8 4.9 5.0 3.8 3.7 3.1 3.8
##  [325] 5.6 2.7 4.1 4.7 4.1 4.1 3.5 4.6 4.4 6.0 5.0 5.1 5.9 4.7 5.1 5.3 3.8 4.3
##  [343] 2.8 4.3 4.4 4.8 4.9 4.5 6.7 5.2 3.4 3.2 3.6 4.2 3.3 6.3 3.7 5.1 3.6 4.2
##  [361] 5.1 5.0 3.8 4.2 4.0 5.0 4.0 6.0 4.2 5.5 5.1 5.6 3.9 4.5 4.2 4.5 2.5 6.0
##  [379] 3.4 4.7 6.1 4.8 4.1 5.5 3.5 5.0 4.5 4.6 5.0 5.3 3.5 3.2 4.1 4.3 4.4 4.7
##  [397] 6.1 4.8 4.9 4.4 2.8 4.6 3.8 5.8 4.7 3.1 3.6 5.8 4.8 4.5 3.2 4.5 4.2 3.1
##  [415] 4.1 4.1 4.8 4.3 5.0 4.0 5.0 4.3 4.3 3.7 3.6 4.5 4.0 3.8 5.2 4.3 2.9 4.3
##  [433] 4.5 4.7 2.9 4.8 3.9 4.2 4.2 3.0 5.0 4.8 3.7 4.4 5.8 4.4 3.8 2.7 5.8 3.6
##  [451] 4.3 2.9 6.4 3.5 4.8 4.8 5.6 6.1 3.3 3.9 4.7 3.5 5.5 5.5 3.9 4.2 4.8 5.3
##  [469] 5.6 4.2 5.7 4.1 4.9 4.1 4.0 3.5 4.7 3.1 5.1 5.0 3.6 3.8 4.2 5.4 2.9 5.7
##  [487] 4.9 3.2 5.1 4.5 2.9 6.0 4.6 3.1 5.2 3.9 6.1 3.4 6.0 6.3 4.8 4.8 3.9 4.7
##  [505] 4.6 3.1 4.7 3.8 2.7 5.2 4.0 4.9 3.9 5.5 4.6 5.6 4.2 4.9 3.8 3.5 5.1 4.8
##  [523] 4.4 3.4 5.2 5.8 4.6 5.2 5.6 3.5 4.4 6.1 5.3 5.9 4.2 5.6 6.0 4.0 4.8 4.4
##  [541] 4.7 4.0 5.3 3.1 2.2 5.7 4.5 5.3 3.6 4.2 5.0 4.7 3.8 4.2 4.3 3.4 3.6 4.4
##  [559] 4.3 4.4 3.3 5.1 5.2 5.1 5.1 3.2 2.9 4.2 3.6 4.2 4.0 5.0 5.4 4.7 2.2 5.1
##  [577] 4.7 4.9 5.4 4.6 4.4 4.0 4.4 4.9 3.4 4.6 4.9 3.8 3.8 4.1 3.5 5.6 3.8 4.6
##  [595] 4.5 6.6 4.2 4.2 3.5 5.2 3.8 4.0 5.5 4.4 3.9 3.5 4.1 5.2 3.3 4.1 4.0 5.3
##  [613] 3.7 6.1 4.4 5.0 4.4 4.8 5.3 5.4 5.2 3.4 3.6 5.1 2.4 6.3 5.8 5.9 5.3 4.9
##  [631] 5.1 3.4 3.9 4.5 5.0 5.0 6.3 4.8 5.5 5.1 4.0 4.8 4.6 3.0 4.7 3.9 4.6 4.2
##  [649] 2.7 4.6 3.8 4.5 4.8 4.7 5.2 5.2 3.5 3.7 3.9 5.0 4.4 4.7 4.5 2.0 4.9 3.8
##  [667] 3.5 5.8 5.7 4.3 3.6 4.8 5.4 3.9 5.7 3.9 4.0 4.7 4.5 3.3 5.0 4.3 4.4 5.4
##  [685] 4.9 3.0 4.8 3.9 3.7 3.8 4.0 3.3 3.6 6.0 5.1 5.3 5.7 3.2 4.3 3.9 4.6 3.0
##  [703] 4.4 5.8 4.5 6.1 4.3 5.6 5.1 4.3 5.0 4.7 5.4 5.6 4.3 4.7 6.3 5.8 3.6 3.3
##  [721] 4.1 5.6 2.4 4.8 4.9 2.8 5.3 4.9 4.1 4.5 2.5 4.3 5.5 4.0 3.3 3.8 3.8 4.5
##  [739] 3.9 4.7 5.8 3.6 5.2 4.9 4.4 3.9 5.7 2.6 4.6 3.7 4.2 5.5 6.2 3.6 4.2 4.0
##  [757] 5.6 3.7 5.1 4.1 4.0 4.7 3.8 5.9 4.4 3.4 4.9 4.1 4.4 3.6 6.2 3.6 4.3 4.1
##  [775] 2.9 3.6 3.1 5.6 6.6 5.7 4.8 5.6 3.9 5.3 3.9 5.0 3.9 4.8 5.5 4.9 3.0 3.6
##  [793] 6.3 3.3 5.4 5.1 6.2 4.7 4.5 3.5 4.9 2.5 3.6 5.1 5.9 3.4 5.0 4.1 4.2 3.8
##  [811] 3.8 3.8 3.5 3.6 4.3 4.6 4.9 5.5 3.2 5.8 4.6 4.2 4.9 5.0 4.2 4.0 5.2 3.1
##  [829] 3.9 4.5 5.0 4.4 2.9 4.8 4.5 4.1 5.5 4.0 4.7 3.9 5.2 5.8 4.1 4.6 4.6 3.6
##  [847] 3.6 4.2 4.8 5.7 4.6 3.6 4.3 3.9 4.6 4.3 4.7 3.6 3.8 5.9 5.2 4.7 5.8 5.0
##  [865] 4.9 4.6 5.1 3.9 2.5 4.5 6.2 5.0 5.6 2.8 3.0 3.6 4.4 5.1 5.5 4.6 5.9 4.5
##  [883] 5.0 3.5 6.9 4.3 4.4 3.9 2.8 3.0 5.2 4.1 5.1 5.4 6.6 5.5 5.3 3.0 4.2 3.8
##  [901] 3.8 4.7 3.3 5.4 1.9 4.6 4.7 3.6 3.0 4.9 6.0 3.3 4.5 4.4 6.6 2.9 3.8 3.9
##  [919] 7.0 5.8 3.3 3.0 5.8 4.1 5.1 4.4 4.1 4.5 5.3 4.1 4.7 5.3 6.0 4.4 4.2 5.1
##  [937] 4.6 5.4 3.9 4.8 3.7 2.9 7.5 5.5 4.3 3.8 4.0 4.4 4.0 4.9 4.2 4.9 5.3 5.2
##  [955] 5.2 5.2 4.2 4.5 4.4 6.3 4.5 5.2 4.6 4.2 4.6 5.0 4.3 4.8 4.8 2.3 5.7 4.0
##  [973] 5.7 4.9 4.6 6.4 4.1 3.2 4.3 3.7 3.5 4.5 4.1 5.6 5.9 4.3 5.0 5.5 4.9 3.8
##  [991] 5.9 4.8 5.1 3.8 5.0 4.3 4.7 4.5 5.1 6.8
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
##    [1] 5.1 6.1 4.7 5.8 3.1 3.6 5.6 4.1 4.2 5.2 5.8 5.1 4.6 2.4 4.8 4.7 3.2 5.1
##   [19] 3.5 3.2 3.9 4.3 3.4 4.0 4.4 4.1 4.3 5.3 5.6 4.7 5.4 5.6 5.6 4.3 3.5 4.9
##   [37] 5.2 4.3 4.2 4.5 5.6 3.6 4.2 5.4 3.3 5.0 4.0 3.0 4.2 4.1 4.6 3.9 4.7 4.5
##   [55] 4.0 4.1 4.4 5.6 6.2 4.5 5.1 5.6 4.5 3.2 5.1 4.6 5.0 5.1 4.6 4.7 3.4 5.1
##   [73] 3.5 4.8 3.4 4.3 5.0 4.8 3.6 2.9 5.3 3.8 3.2 4.6 4.5 4.0 4.1 5.3 3.3 5.1
##   [91] 3.9 5.5 5.3 5.1 3.2 5.2 3.4 5.4 5.3 3.7 4.3 3.1 4.3 3.0 4.4 4.9 4.3 5.1
##  [109] 5.0 5.1 5.7 5.8 4.5 2.7 4.7 3.5 4.0 6.0 5.0 4.0 5.4 4.2 3.5 5.1 5.5 5.9
##  [127] 3.5 3.6 5.3 3.5 3.3 5.2 4.4 3.0 4.8 4.4 4.0 3.9 4.6 5.1 3.8 3.6 5.1 5.3
##  [145] 5.5 3.9 5.2 5.6 2.7 3.7 6.1 3.0 4.8 4.5 5.7 5.2 5.3 5.0 3.7 5.3 4.9 3.8
##  [163] 5.5 5.0 5.7 3.1 5.3 4.3 6.2 4.3 3.9 4.0 4.8 5.1 3.3 4.8 2.8 5.2 4.1 3.8
##  [181] 3.3 5.9 3.0 3.4 4.5 3.2 4.3 2.7 3.4 4.8 4.8 4.3 2.7 5.6 4.1 3.0 5.8 2.3
##  [199] 4.9 4.2 5.8 2.9 5.5 5.7 4.3 6.2 5.4 5.9 4.6 3.7 6.6 5.3 4.7 5.1 5.9 4.8
##  [217] 5.7 4.0 5.2 4.3 6.2 6.2 3.1 4.4 5.4 4.2 3.0 5.0 4.9 4.9 5.1 4.8 4.3 4.6
##  [235] 5.1 2.5 3.7 3.4 4.4 5.3 4.7 4.7 4.8 5.0 5.8 4.7 4.4 3.8 5.2 4.9 4.9 4.5
##  [253] 4.3 4.7 4.2 2.9 5.0 4.7 5.5 4.0 5.3 4.2 4.7 4.5 5.2 5.3 4.3 4.5 4.2 5.8
##  [271] 4.6 5.5 3.4 3.0 4.8 4.9 3.7 3.5 4.5 3.9 5.0 3.2 3.9 4.4 2.6 4.1 3.1 5.2
##  [289] 4.0 3.2 6.0 4.1 5.5 2.6 3.6 4.3 5.5 5.7 5.3 3.9 3.5 5.3 5.3 3.9 5.0 6.4
##  [307] 5.7 4.8 4.7 6.3 3.7 4.0 4.5 4.4 3.5 5.4 5.8 5.4 3.2 4.7 4.5 5.8 3.0 3.5
##  [325] 4.4 5.7 4.7 4.3 4.3 5.3 4.5 5.7 5.5 4.0 4.9 4.6 4.4 3.4 5.4 5.7 4.1 3.9
##  [343] 3.6 4.9 3.7 4.1 4.6 4.1 3.9 5.9 3.0 4.3 5.8 5.7 3.3 2.9 2.4 5.0 4.1 4.7
##  [361] 4.9 5.0 6.1 5.8 4.7 5.3 4.0 3.9 3.1 4.5 3.8 4.9 4.0 4.1 3.8 4.2 3.1 4.8
##  [379] 5.9 3.6 4.2 4.3 6.1 5.1 5.5 3.8 4.9 3.2 4.0 4.4 3.2 4.5 5.9 2.3 3.3 4.4
##  [397] 5.3 3.4 4.1 4.3 5.3 4.8 3.7 3.4 4.1 5.4 3.9 4.5 6.5 4.6 4.2 4.4 3.7 3.8
##  [415] 6.0 5.8 5.8 4.0 6.1 4.1 4.3 5.0 3.9 5.1 3.7 5.4 4.3 4.3 5.6 5.7 3.7 4.1
##  [433] 4.2 5.6 4.9 4.5 6.8 3.1 5.8 4.2 4.3 3.5 3.8 4.8 3.9 5.0 3.3 5.0 3.6 3.8
##  [451] 5.0 5.0 5.0 6.3 6.3 4.7 4.0 4.1 5.6 2.9 3.9 4.9 4.3 4.0 5.2 3.6 3.6 4.6
##  [469] 3.5 4.1 4.5 4.8 5.1 2.6 3.1 3.6 4.5 6.8 4.9 4.4 2.8 4.7 3.6 4.9 3.8 5.3
##  [487] 3.4 4.7 4.9 4.3 4.1 5.8 4.9 5.9 2.7 3.7 3.7 5.3 3.7 2.2 4.7 3.8 3.9 5.2
##  [505] 4.6 4.2 5.7 4.0 4.0 5.9 4.8 2.6 5.9 4.9 3.7 6.4 5.2 5.7 5.1 5.4 4.0 5.9
##  [523] 5.2 4.9 5.0 2.2 3.7 4.6 3.1 4.7 4.4 3.1 4.6 3.3 6.1 4.8 5.0 3.7 5.8 3.6
##  [541] 5.4 2.3 4.2 5.7 4.3 4.3 3.9 4.3 4.4 4.3 3.5 4.8 5.4 2.8 4.2 4.2 5.7 6.3
##  [559] 2.5 5.3 3.6 4.2 4.8 4.9 4.1 2.7 5.6 3.8 4.0 4.5 5.1 4.3 5.9 5.7 4.5 4.3
##  [577] 4.3 3.8 2.1 3.5 3.1 3.5 4.7 4.7 5.2 3.4 6.5 4.2 4.2 4.9 3.4 4.7 5.0 5.4
##  [595] 5.1 4.3 3.4 4.0 3.7 4.8 4.5 2.6 2.8 4.9 5.9 5.2 5.6 3.2 4.1 3.8 4.2 4.5
##  [613] 5.8 6.1 5.2 4.4 3.9 4.1 4.8 3.5 4.4 5.5 3.9 4.3 4.5 4.7 5.4 3.2 6.1 4.4
##  [631] 3.1 3.2 5.4 5.1 4.3 4.2 3.7 3.5 5.5 4.6 4.6 4.2 4.1 2.4 4.5 4.2 2.5 5.5
##  [649] 3.9 5.1 4.8 4.6 6.4 3.7 5.9 6.2 4.5 5.1 5.2 3.4 3.3 6.5 5.0 4.8 4.4 4.9
##  [667] 4.7 4.2 4.2 4.9 4.5 4.5 5.0 5.4 6.0 4.6 5.2 4.8 2.4 3.4 5.0 6.1 4.2 4.2
##  [685] 4.3 3.8 3.5 5.1 4.7 4.2 5.0 6.1 4.6 5.2 4.7 6.1 4.5 3.8 3.9 4.5 5.1 4.6
##  [703] 5.3 4.5 4.7 4.9 4.1 5.5 4.1 4.7 4.6 4.6 5.4 4.6 4.9 5.0 2.7 6.1 4.0 3.0
##  [721] 3.9 4.7 3.3 3.4 5.3 5.4 5.3 4.1 5.7 4.5 4.4 5.3 5.7 3.7 3.7 5.1 4.0 3.2
##  [739] 4.5 4.3 2.5 3.7 6.1 6.0 3.6 4.8 4.3 5.6 3.5 2.3 3.7 4.6 3.7 4.9 4.8 4.5
##  [757] 4.8 3.8 5.2 3.7 4.1 5.6 4.6 4.8 4.3 5.4 2.8 3.6 4.8 5.4 4.6 3.8 3.6 6.3
##  [775] 5.4 2.5 4.4 5.8 5.3 5.6 3.7 4.3 4.9 5.7 4.4 6.2 4.9 4.2 3.4 2.9 5.4 4.0
##  [793] 4.4 4.0 5.0 4.3 5.3 3.7 6.0 4.8 4.5 4.8 4.8 4.7 3.7 3.8 4.6 5.7 4.4 4.8
##  [811] 4.2 3.8 3.6 7.1 4.0 3.2 4.4 3.8 4.3 5.8 5.1 4.5 5.6 4.3 6.5 6.0 4.6 4.2
##  [829] 4.1 4.4 5.0 5.4 6.3 4.3 3.6 4.5 4.6 5.1 3.5 6.0 5.4 5.0 3.4 5.3 5.2 4.6
##  [847] 5.1 3.4 4.4 4.2 3.8 5.2 3.6 5.1 3.6 2.6 4.1 5.5 5.4 4.7 4.1 3.8 5.2 4.1
##  [865] 3.5 6.2 4.2 3.9 4.4 3.2 5.5 3.8 5.6 3.4 3.3 5.7 5.1 4.2 4.4 3.4 3.6 4.4
##  [883] 5.9 3.7 4.3 2.8 4.5 4.5 3.7 3.3 4.6 3.5 4.2 4.9 5.1 6.2 3.6 5.8 4.9 2.5
##  [901] 4.3 5.4 5.5 5.1 4.8 4.6 5.2 4.6 3.6 3.2 4.6 3.2 3.7 4.5 3.9 3.9 3.1 4.4
##  [919] 5.0 3.9 3.9 3.9 3.8 3.6 5.9 4.3 5.3 5.4 5.6 4.4 2.8 3.6 3.7 5.2 4.2 2.6
##  [937] 4.2 3.5 4.0 3.7 5.0 3.8 4.5 4.6 4.6 3.4 3.6 5.6 3.5 4.6 3.8 3.8 5.0 4.3
##  [955] 6.5 3.5 5.7 2.8 5.3 4.4 3.7 4.4 4.4 5.3 4.8 5.3 4.6 5.0 5.6 6.3 4.0 4.6
##  [973] 4.1 4.2 3.4 3.6 4.2 4.5 5.4 4.5 4.9 5.0 4.3 3.3 3.2 5.0 7.0 3.1 4.1 2.4
##  [991] 5.9 3.7 4.4 4.0 5.5 4.7 4.8 4.2 3.4 4.5
## 
## $func.thetastar
## [1] -0.014
## 
## $jack.boot.val
##  [1]  0.52551320  0.30614525  0.35172414  0.16617647  0.10888252 -0.06186441
##  [7] -0.21146132 -0.34568690 -0.42479564 -0.52492918
## 
## $jack.boot.se
## [1] 1.019955
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
##    [1] 4.4 6.0 4.7 6.1 4.4 4.6 3.6 4.5 4.0 5.2 4.9 4.4 4.5 5.3 3.9 4.0 3.2 3.9
##   [19] 4.2 4.3 4.5 2.6 4.8 4.2 5.1 4.6 5.5 3.9 4.7 4.6 4.4 5.5 2.9 5.3 4.8 4.5
##   [37] 4.6 5.1 3.0 5.7 3.8 4.0 5.3 4.3 4.3 5.0 4.5 4.9 3.7 4.5 3.2 5.3 3.9 5.2
##   [55] 4.9 5.8 3.5 4.5 3.4 3.6 4.2 5.3 4.9 5.1 5.1 3.9 5.2 3.9 4.1 4.4 4.9 4.9
##   [73] 6.0 4.7 3.6 4.8 4.5 4.2 4.0 4.2 5.5 2.8 4.9 5.4 5.1 3.4 5.8 5.2 3.6 5.0
##   [91] 4.4 4.6 4.5 5.6 5.7 5.6 5.5 5.6 4.9 4.5 5.2 4.9 4.9 4.2 3.9 5.2 6.4 3.4
##  [109] 4.4 2.5 4.5 6.4 3.5 5.0 5.3 3.6 3.4 4.4 5.3 4.1 5.3 4.7 5.8 5.4 5.5 4.3
##  [127] 4.0 4.1 4.3 4.9 5.1 4.5 4.8 4.8 3.2 4.8 4.9 4.0 4.0 4.3 5.6 3.4 3.8 2.9
##  [145] 4.4 3.2 3.7 4.0 4.1 4.4 4.0 5.4 4.2 3.9 2.8 4.2 5.0 3.9 2.8 4.9 5.3 3.8
##  [163] 4.2 5.4 4.2 3.1 5.1 2.8 4.7 2.8 4.7 6.0 5.1 6.0 3.9 4.4 5.5 6.2 5.5 5.0
##  [181] 4.9 5.8 4.8 4.6 5.1 4.1 3.0 5.7 3.2 4.3 4.4 3.1 5.9 4.5 5.3 4.3 3.8 5.3
##  [199] 5.2 3.0 4.9 3.5 4.5 3.8 2.4 4.0 3.9 4.3 3.0 3.3 2.6 4.6 3.4 4.2 6.0 6.0
##  [217] 4.8 5.3 4.9 4.6 5.7 2.9 4.3 2.6 5.0 4.3 5.2 4.7 2.6 3.6 5.4 5.3 4.6 4.6
##  [235] 4.5 4.3 4.3 4.4 4.8 5.6 5.3 6.3 4.4 4.2 4.5 5.2 5.8 3.6 5.3 6.6 3.4 4.9
##  [253] 3.7 4.4 3.5 6.7 3.6 4.2 3.9 5.7 5.1 3.0 4.3 4.6 3.4 5.6 5.2 2.3 3.9 3.3
##  [271] 4.8 3.2 5.0 5.8 4.1 3.9 4.4 2.5 4.7 5.1 6.1 3.1 3.7 4.4 4.4 5.1 3.7 3.2
##  [289] 3.8 5.1 4.1 5.3 5.7 4.8 4.5 6.2 3.8 4.0 4.2 4.6 4.5 5.2 4.2 4.4 4.0 5.6
##  [307] 6.2 5.6 4.8 5.4 5.6 5.2 4.3 3.2 2.8 4.8 4.3 4.7 3.3 3.2 4.4 4.8 3.5 3.4
##  [325] 4.5 3.2 5.1 4.8 5.6 5.0 4.9 5.9 4.7 5.7 5.1 4.5 6.2 4.5 5.6 3.0 4.3 5.7
##  [343] 2.8 4.5 4.2 4.8 4.7 4.2 4.8 3.1 5.5 3.8 3.8 3.4 3.0 4.6 5.0 2.4 4.6 4.1
##  [361] 3.6 4.8 5.1 4.2 4.2 4.6 6.2 4.8 5.0 5.1 3.6 4.7 5.0 4.9 4.4 5.5 3.3 2.8
##  [379] 3.6 5.6 4.2 4.5 4.2 5.0 4.6 3.2 4.2 5.9 4.3 5.4 4.5 4.9 7.2 5.0 2.5 5.2
##  [397] 4.4 3.6 4.6 4.5 5.0 6.0 4.2 5.8 3.8 3.4 4.3 4.6 3.4 2.7 4.4 5.0 4.9 4.1
##  [415] 4.5 4.1 4.3 5.4 3.7 4.7 4.9 3.3 6.1 3.5 3.8 4.9 4.3 4.0 5.1 4.6 4.0 4.1
##  [433] 4.5 1.9 5.7 3.7 4.5 4.4 6.3 4.1 4.6 4.3 5.8 3.3 4.0 4.2 4.9 3.6 3.9 3.3
##  [451] 3.2 4.7 4.6 4.4 4.4 4.7 4.3 2.5 4.5 4.2 5.5 5.0 6.2 5.7 3.8 4.9 3.8 5.4
##  [469] 2.6 6.2 4.2 5.9 5.9 3.1 3.8 4.4 4.8 4.1 5.2 4.0 4.8 2.8 4.4 4.4 5.2 4.5
##  [487] 4.3 4.3 4.7 4.4 4.4 3.5 6.0 5.2 5.2 4.8 5.1 3.6 4.7 5.1 5.0 5.3 3.2 5.9
##  [505] 5.5 4.7 6.0 3.7 4.8 4.5 4.8 4.1 4.9 4.4 4.7 4.7 3.4 3.8 3.3 4.6 5.8 4.1
##  [523] 4.8 4.7 4.5 5.8 4.7 5.2 3.4 5.0 5.4 3.8 5.5 3.8 4.9 4.8 4.1 6.0 4.8 3.4
##  [541] 3.3 4.2 4.5 3.6 3.8 4.2 2.9 5.1 5.4 4.2 5.4 4.0 4.0 2.4 4.0 4.9 5.4 4.4
##  [559] 5.8 5.6 4.7 3.0 3.7 3.6 2.9 5.5 2.9 3.8 4.9 5.1 3.9 4.4 3.9 2.8 5.2 4.4
##  [577] 4.8 5.8 4.8 3.9 4.3 4.1 3.6 3.5 4.7 5.5 5.3 4.2 2.8 4.2 4.3 5.3 4.6 5.8
##  [595] 4.2 3.7 2.7 3.8 3.3 4.7 3.3 6.4 6.9 4.5 3.7 4.2 4.1 5.1 4.0 5.4 5.8 2.8
##  [613] 6.5 6.0 4.4 5.0 5.4 3.6 4.8 5.6 3.9 5.1 3.6 4.3 5.2 3.9 5.9 4.4 4.9 4.9
##  [631] 3.6 3.6 3.4 3.7 4.1 4.2 4.9 5.6 4.5 6.0 4.1 3.9 5.1 5.2 4.9 3.9 4.8 4.8
##  [649] 5.1 4.8 3.8 5.5 5.0 4.5 4.7 6.0 4.2 4.7 4.1 4.9 4.5 3.8 5.3 4.7 6.3 4.5
##  [667] 3.5 4.9 5.1 3.0 4.9 5.6 6.1 3.8 4.6 4.5 4.4 6.1 2.3 3.1 4.4 5.7 5.4 4.0
##  [685] 5.4 3.1 3.3 2.9 3.8 3.6 4.5 4.3 5.2 4.5 4.4 5.2 6.1 3.8 5.7 5.2 3.6 5.1
##  [703] 3.2 4.5 2.5 6.4 5.1 5.0 2.8 4.0 4.6 5.2 5.5 4.4 5.2 4.1 4.7 3.3 5.5 3.5
##  [721] 1.9 6.1 4.2 4.3 3.7 5.0 4.8 4.7 5.1 4.5 4.1 4.0 4.0 5.2 4.6 4.3 4.6 5.2
##  [739] 3.4 4.5 4.8 5.7 4.4 4.2 5.6 4.5 4.8 5.1 4.4 2.3 5.8 4.1 3.3 4.6 3.6 5.7
##  [757] 6.2 5.7 4.0 4.5 3.3 3.6 4.2 4.3 4.7 4.6 5.5 4.3 4.5 6.0 5.4 4.2 4.5 4.3
##  [775] 4.3 6.5 3.7 2.3 6.6 5.4 3.6 5.5 4.8 6.1 5.0 3.5 5.1 3.7 4.9 3.6 3.9 4.3
##  [793] 3.0 3.6 5.3 4.5 6.0 3.9 3.9 5.8 6.4 3.5 4.8 5.5 4.1 5.1 5.1 5.1 4.9 5.1
##  [811] 4.4 4.0 3.7 5.5 4.9 3.5 3.7 3.1 3.6 3.4 4.3 2.4 5.4 4.7 5.6 4.0 4.5 4.5
##  [829] 4.5 4.1 4.7 5.4 5.2 5.5 4.6 6.1 4.7 5.0 4.4 2.9 4.8 5.0 3.9 6.4 4.1 5.3
##  [847] 5.2 4.6 5.4 5.4 2.9 5.8 4.6 5.8 5.5 5.3 3.3 5.7 3.1 4.4 5.9 4.9 3.3 4.5
##  [865] 5.4 4.5 5.2 3.2 4.4 7.0 3.4 5.8 4.5 4.3 2.6 4.7 4.9 3.7 4.7 5.1 2.8 4.4
##  [883] 4.5 3.8 3.4 3.5 5.2 3.3 5.4 5.4 4.8 3.8 5.1 4.2 4.8 4.1 3.3 4.8 5.8 4.8
##  [901] 3.8 5.5 4.7 3.6 3.2 6.0 4.3 5.3 3.9 4.3 3.3 2.4 3.5 2.1 3.0 6.3 4.5 3.6
##  [919] 3.8 4.3 3.2 5.5 4.0 3.5 3.9 4.1 4.5 6.9 4.3 3.9 4.7 3.7 5.0 5.2 5.6 4.3
##  [937] 4.2 3.7 4.4 5.0 4.6 5.0 2.0 4.8 5.1 3.6 3.1 3.5 5.2 5.7 5.7 3.2 4.0 3.6
##  [955] 3.6 5.4 4.3 4.5 5.8 4.5 4.3 5.3 4.8 5.7 3.8 5.1 5.1 3.1 4.5 5.3 4.8 5.2
##  [973] 4.4 4.7 4.8 4.6 3.5 4.6 4.5 3.9 5.6 4.1 4.0 2.4 5.0 4.2 4.5 4.7 3.4 5.1
##  [991] 4.5 3.8 4.2 4.7 4.4 4.7 5.4 2.2 3.8 3.5
## 
## $func.thetastar
## 72% 
##   5 
## 
## $jack.boot.val
##  [1] 5.4 5.4 5.4 5.2 5.1 4.9 4.9 4.8 4.7 4.5
## 
## $jack.boot.se
## [1] 0.9104395
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
## [1] -0.1699313
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
##   4.736668   7.950992 
##  (2.048065) (3.626810)
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
## [1]  0.9539188  0.9767475  1.5318866  0.2929803  0.3473668 -0.2335794
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
##    [1]  1.201763235 -0.861005575  0.104332617 -0.228194167  0.113456031
##    [6]  1.157943223 -0.430915229  0.811841715 -0.543412820 -0.492914701
##   [11]  0.210874974  0.483818880 -0.270185486 -0.286562237  0.448554766
##   [16] -1.342834640  1.305670232 -0.152956062  0.237273893 -0.097435686
##   [21] -0.040547825  0.336303734 -0.779237745 -0.244907615  0.370626924
##   [26] -0.290091199 -0.406579551  0.150945969  0.016235859 -0.387593947
##   [31] -1.179854226 -0.217922869 -0.352626103 -0.483988420  0.080753874
##   [36] -0.894758144 -0.504187919 -0.095249350 -0.185533657  0.342479562
##   [41]  0.092368087  0.161873837  0.395193561  0.153207146  0.085112766
##   [46] -1.114759170  0.185700697  0.064276584 -0.483198236  0.426380859
##   [51] -0.646099413 -0.086770157  0.260664046 -0.604248390 -0.643065358
##   [56] -0.366559622  0.266157529 -0.090115894 -0.370044626 -0.597178986
##   [61] -0.696789421 -0.117944114  0.002124252 -0.207680951 -0.339964422
##   [66] -0.502803816  0.274053072 -0.539020817 -0.396545323  0.367717092
##   [71] -0.041945275 -0.094404298 -0.018922937  0.475020570 -0.262810901
##   [76]  0.979460551  0.425269629  0.287256323 -0.642045117 -0.730055102
##   [81]  0.096006165  0.611328456 -0.960637369 -0.331427960 -0.093282413
##   [86]  0.250050472  0.242768068  0.156735579  1.065977366 -0.178942239
##   [91] -0.331129628 -1.776632345 -0.628321620 -0.632012267 -0.041593869
##   [96] -0.642107934 -0.039893325  0.231684320 -0.234696876  0.119115677
##  [101] -0.363503055  0.823390726 -0.419557158  0.072773390 -0.326228836
##  [106]  0.364370250 -0.193625052  0.229366051 -1.133362133  0.164554267
##  [111] -0.355977110 -0.669073837  0.755655398 -0.341239379 -0.375445446
##  [116] -0.239758238  0.487497124 -0.186681292 -1.960032555 -0.361602328
##  [121] -0.840801399 -1.101044033  0.075340482 -0.215081874 -1.356634112
##  [126]  0.207868565  0.218008096 -0.854527636  0.109813886  0.033104616
##  [131] -0.484239645 -0.688061741 -0.538743103 -0.099276519  0.379725431
##  [136] -0.087207688 -0.706253027 -0.539231325  1.084361976  0.016326755
##  [141] -0.197968347 -0.723343529 -0.716733208  0.112707866 -1.041591859
##  [146] -0.299799042  0.264478981 -0.050741918 -0.556345608 -0.097934808
##  [151] -0.703262375 -0.361084838  0.413239545  0.483253966 -0.458539792
##  [156] -0.167038643  0.347802229 -0.064433200 -0.241097619 -0.742116826
##  [161]  0.266789783  0.030317814  0.194133484 -0.115979232 -0.149226328
##  [166] -0.255321434 -0.562262737  0.113391204  0.118381741  0.407309152
##  [171] -0.983621346  0.542346588  0.018026953  0.289029892 -0.573206212
##  [176] -0.691064707 -1.320346182 -0.559553365 -0.450588439 -0.219644930
##  [181]  0.660284742 -0.763054659  0.071417704 -0.047823002 -0.745041213
##  [186]  0.320410721  0.381224656 -0.659457937 -1.067630341  0.090024774
##  [191] -0.627954606  0.593202512 -1.427177256 -1.898356711  0.214116135
##  [196] -0.862704135  0.183452169 -0.363488275 -0.785718292 -0.394719998
##  [201] -1.221782189 -1.182559056 -0.729467819 -0.391206760  0.569505220
##  [206] -0.021503654 -0.049862057 -0.354699955 -0.611001938 -0.201357205
##  [211] -0.118304184 -0.690341349 -0.252972196 -0.826335646 -0.506152398
##  [216] -0.079749052 -0.381939592 -0.252690913  0.004664174 -0.609134814
##  [221]  0.539941419 -0.512512796 -0.103437633 -0.117621181 -0.732273065
##  [226] -0.545514234 -0.331893562  0.321245043  0.093551480 -0.113348794
##  [231] -0.821277340  0.291803309 -0.001495972 -0.520664838 -0.467855855
##  [236] -0.477825084  0.047125370 -0.436448877 -0.293347943  0.085744415
##  [241] -0.091967470 -0.267928320  1.099691095 -0.263331134 -0.148974274
##  [246] -0.645832593 -1.247450645 -0.946317934 -0.905981275  0.022748275
##  [251]  0.313183493 -0.436811848  0.586936290 -0.598686568  1.143074672
##  [256] -0.323831273 -0.084617856  0.468356106 -0.268461013 -0.107236305
##  [261]  0.040894976 -0.647950237  0.018975960  0.262488794  0.394121749
##  [266]  0.177451363 -0.493362935 -0.648311849  0.341522022  0.070010422
##  [271]  0.406801620 -0.261126537  0.424075041  0.534687466  0.034567373
##  [276] -0.162639304  0.688449186  0.126460032 -0.500973807 -0.394284439
##  [281] -0.634451079 -0.166098416 -0.665012073 -1.132234460 -0.770729314
##  [286]  0.346502469 -0.708837160 -0.268130283 -0.001699364 -0.259415439
##  [291] -0.637140913 -0.659999457 -0.150965675  0.542793772 -0.983368400
##  [296] -0.391685293 -0.404413211 -0.190354627 -0.129864638  0.489095551
##  [301]  0.128224332  0.122089101 -0.649037480  0.483301677  0.385077406
##  [306] -0.899416522 -0.042074809  0.596578249  0.053602257 -0.192049273
##  [311] -0.411868274  0.433142812  0.177395298 -0.269393095 -0.236972804
##  [316] -0.564365289  0.095683131  0.721200708  0.707177561  0.077203375
##  [321] -0.113607379 -0.282576836 -0.374287904 -0.870302860 -1.127900600
##  [326] -0.013115937 -0.701860250 -0.226817820  0.450186936  0.567205739
##  [331] -0.212723238 -0.224890391 -0.507091089 -0.522075408  0.166294539
##  [336] -1.570618463  0.412744021  0.007507941  0.164358206 -0.281963377
##  [341]  0.888660464  0.173762582 -0.872323566 -1.401278515 -0.049692235
##  [346] -1.846135694  0.524469964  0.102862087 -0.285652288 -0.280098767
##  [351] -0.619011756 -0.158589682 -0.670489192 -0.167001768 -0.300477793
##  [356] -0.867258138 -0.505555777 -0.348429014 -0.930927120 -0.097754152
##  [361] -0.089005787 -0.620984548  0.364435470 -0.294095690 -0.240731750
##  [366]  0.746074575 -0.613338431  0.270217815  0.356207961  0.387982094
##  [371] -0.266236174 -0.409234493 -0.414107938 -0.569991874 -0.665223543
##  [376] -1.171890614 -0.179987678 -0.176419150 -0.382267002 -0.405584490
##  [381]  0.281256194  0.845421313  0.824353588  0.449050368 -0.353990182
##  [386]  0.646049642 -0.570806420  0.091423681  0.046964614 -0.035620467
##  [391] -0.430825333 -0.726457021 -0.787771987  0.362086681 -0.152402300
##  [396]  0.258163610 -0.588234727  0.709408390  0.478103767  0.276807570
##  [401] -0.533687980 -0.257956899  0.013784328  0.943631862 -0.179327403
##  [406] -0.474834638  0.168467039 -0.708272600 -0.156090986 -0.665012073
##  [411] -0.018314447 -0.322300908  0.216380458 -0.848539184 -0.574599035
##  [416] -0.391032381 -0.593749923  0.104272109 -0.042678044 -0.553721139
##  [421]  0.026960007 -0.613027327 -0.137439641 -0.006889218 -0.814039104
##  [426] -0.176390903 -0.023737226  0.258163610 -0.295734792 -0.683403372
##  [431]  0.533350423 -0.495778529 -0.490805741 -0.590955197  1.063553216
##  [436] -0.206061990 -0.891019769  0.411536931 -0.398398858  0.221637835
##  [441] -0.155520487 -0.816403278 -1.152839415 -0.291043343 -0.355977110
##  [446] -1.382626007 -0.369829372 -0.235668290 -0.226838352 -0.289496655
##  [451] -0.779764007  0.667288239  0.456523835  0.064276584 -0.298922867
##  [456] -0.193320860 -0.132848962 -1.090281704  0.904478807  0.172973870
##  [461] -0.111591366 -0.121079946 -0.930668203 -0.474796708  0.204416964
##  [466]  0.052838246 -0.457131682 -0.793635769 -0.827140182 -0.651991462
##  [471]  0.468496832 -0.181128929 -0.833122591 -0.020765610 -0.212723238
##  [476] -0.296965258 -0.109972100 -0.060093210 -0.182662989  0.627048274
##  [481]  0.327077320 -0.321605192 -0.048351762 -0.823523054  0.026451811
##  [486] -0.269786109 -0.788368284 -1.878855219 -0.370567020 -0.207722748
##  [491] -1.208680521  0.198439470 -0.189671990 -0.403622344  0.598181578
##  [496] -0.056725737 -0.050611219  0.283902949 -0.116880023 -0.488797404
##  [501] -0.070006269 -0.470503950 -1.369014803  0.399935706  0.131578239
##  [506] -0.049012479 -0.037422701  0.355709356 -0.475949575  0.189085386
##  [511] -0.308209554  0.582616555  0.623295102 -0.622069343  0.276188313
##  [516] -0.311356734  0.271717454 -0.387307093 -0.021415251 -0.091960892
##  [521] -0.160198481 -0.734338246 -0.259532074  0.495261554 -0.645894457
##  [526] -0.076713032  0.304193686 -0.961995622 -0.513102349 -0.072641403
##  [531]  0.323674849 -0.593959520  1.034965590 -0.206004750  0.647520162
##  [536] -0.278594317 -0.647300686 -0.069596254 -0.462764585 -0.896331196
##  [541]  0.257536542 -0.003123178  0.238703516 -0.104736238 -0.077318961
##  [546] -0.757518973 -0.294762780 -0.271136447  0.273370261 -0.134931548
##  [551]  0.747976319  0.008876829 -1.019501007  0.869503632 -0.934128531
##  [556] -1.467428320  0.105939163 -0.033351433 -0.093759919 -1.353517238
##  [561] -1.312212935 -0.463683366 -0.351324297  0.044779301  0.458383443
##  [566] -0.527023426 -0.432106684 -0.099584298 -0.665336378 -0.374289487
##  [571] -0.344806609 -1.555135589 -0.849901713  0.150696053  0.141955525
##  [576] -0.929099997 -0.388597543  0.255271167 -0.208463576 -0.320429836
##  [581] -0.338856595 -0.235602576 -0.719030296 -1.099507350  0.523719164
##  [586] -0.993904520  0.982675140 -0.276889376 -0.928448310  1.197218912
##  [591] -0.362409499 -1.126890606 -0.717322567 -0.833644707 -0.932175965
##  [596] -0.111627631  0.039714099  0.639279193 -0.021750216 -0.065782035
##  [601] -0.096763118 -0.149226328 -0.791719437  0.204704531 -0.011224961
##  [606]  0.179087974  0.466042498  0.100087592  0.545289866 -0.095032648
##  [611]  0.604936862 -0.501424247 -1.286020134 -1.254965264 -0.253915520
##  [616] -0.729204010  0.440219185  0.447768240 -0.095850779 -0.082266783
##  [621] -0.324992704 -0.506599873 -0.760485148 -0.499916502 -0.250817087
##  [626]  0.003933547 -0.209521588  0.006574690 -1.384993883 -0.580696189
##  [631] -0.455932322  0.019371520 -1.242273984 -0.376012625 -0.999290251
##  [636] -0.120673288  0.975231503 -0.984855770 -0.099328399  0.187624374
##  [641] -0.153346458 -0.133530107 -0.283771578 -0.248935036 -1.435632369
##  [646]  0.415055501 -0.160358154 -0.736313230 -0.033890648 -0.657697391
##  [651]  0.501790817 -0.411350237 -0.249922513  0.445853343  0.033410272
##  [656]  1.214991202  0.062900749 -0.279712375 -0.251736215 -0.558281549
##  [661] -0.118388192 -0.185352864 -0.778642985  0.305974937 -0.483453764
##  [666] -0.345223146 -0.746923366 -0.068174378 -0.132767846 -0.456764026
##  [671]  0.739573747 -0.076586637 -0.189283063 -0.934128531 -0.689376326
##  [676] -0.340180469 -0.073022420 -1.270372336 -0.501925333  0.978763731
##  [681] -0.648444595 -0.809743610 -0.535110238  0.448952452 -0.197581768
##  [686]  0.372686488  0.820811608 -0.571341808 -0.912694836 -0.502812915
##  [691]  0.309562390 -0.202249958  0.266326856 -0.451947574 -0.100435210
##  [696] -0.208433371 -0.904908975  0.363968267 -0.855569650  0.251225022
##  [701] -0.720912644 -0.255805242 -1.065437228  0.512209730  0.203014553
##  [706] -0.085828978  0.839005428  0.004824186  0.448313365 -0.314063466
##  [711] -0.133978812 -0.194161792  0.767641665 -0.748107508  0.983848244
##  [716] -0.703956890 -0.530575748 -0.146854447 -0.626357024 -0.903828011
##  [721]  0.343606598 -0.301701831  1.251864902 -0.581574997 -0.062433872
##  [726] -0.049862057 -0.379155915  0.390798531 -0.962206752  0.316583156
##  [731] -0.445287784  0.476887159 -0.943459699 -0.219253153  0.307668199
##  [736] -0.617990275 -0.579780483  1.191901818  0.184072317 -0.506475719
##  [741] -0.210657240 -0.297461228 -0.394012802 -0.773626331  0.427703896
##  [746] -1.264286222 -0.982637545 -1.144698747 -0.407553593 -0.238206586
##  [751] -0.224947796 -0.564256561  0.103220861  0.241258225 -0.117456776
##  [756] -0.950541786 -0.295598259 -0.355334265  0.249654675  0.283242521
##  [761]  0.380218496  0.270296692 -1.518091751  0.616259529 -0.711993692
##  [766] -0.029680005 -0.618027738  0.165734812 -0.503333214 -0.887014719
##  [771] -0.028506347 -0.076098791  0.612869260  0.151564833 -0.558823638
##  [776] -0.643413149 -0.574372214 -0.552462574 -0.244702850 -0.502860961
##  [781] -0.998296935 -0.153852984 -0.222811108 -0.295598259 -0.651261822
##  [786] -0.135726421 -0.482253288  0.402867501 -0.075713008  0.398283740
##  [791] -0.149596400  0.275034945  0.186023275 -0.721354496 -0.830297599
##  [796] -0.542813778 -0.251245612 -0.077448770 -0.213891240 -0.177570387
##  [801] -0.098585055 -1.427423594 -0.647021839  0.074596256 -0.967016799
##  [806] -0.888455462 -0.194795504  0.149163759  0.607826920 -0.102743361
##  [811] -0.630980793 -0.379196047 -0.565091306 -0.700189341  0.046273137
##  [816]  0.080549544 -0.319962424  0.069591665  0.579708408  0.796680694
##  [821] -0.725553610 -1.501055446 -0.406707828 -0.153731934  0.155836853
##  [826] -0.778642985  0.221469762 -0.003057723 -0.579558767  0.735404968
##  [831]  0.368786117  0.128081732 -0.024928368 -0.123355293 -0.093890397
##  [836] -0.377932332 -0.017381186 -0.719200641 -0.194161792  0.386433208
##  [841] -1.082501986  0.032053043 -0.238496482 -1.069273802 -0.064394229
##  [846]  0.306555236 -0.114500482 -0.284277747  0.012881325 -0.015001216
##  [851] -0.706430389 -0.442824286 -0.654304838 -0.592250472 -0.627864004
##  [856]  0.602334277 -0.379839857  0.353733418  0.002463575 -0.271990530
##  [861] -0.412614439 -0.496858889 -0.006971843 -0.494388509  0.219774153
##  [866]  0.433788374 -0.971892960 -0.232641061  0.606802684 -0.518301291
##  [871] -0.098085907 -1.800141265 -0.457810472 -0.193452389  0.205832264
##  [876] -0.493119489  0.369457371 -0.523387927 -0.611436940 -0.491467580
##  [881] -0.721012852 -0.644549964  0.269054068 -0.231147850 -0.275875060
##  [886] -0.215344947  0.433684126 -1.989089388 -0.806865097  0.700308847
##  [891] -0.638273360 -0.988313997  0.765724688  0.259057519 -0.266576354
##  [896] -0.338490395 -0.092601230 -0.248329088 -0.087893946  0.252447497
##  [901] -0.786943641 -0.017347401  0.200823994 -0.191411166 -0.035373056
##  [906]  0.186403683  0.419447135 -0.941762785 -0.727791202 -0.083386823
##  [911]  0.247329691 -0.352265329 -0.010080127  0.309331063  0.287929176
##  [916] -0.598804044  0.388386922  0.055879837 -0.551215364 -0.227159813
##  [921]  0.585669138 -0.814946822 -0.133199509 -0.365604048 -0.310884784
##  [926] -0.428925661  0.242329169  0.154547543 -0.523446162 -0.528349453
##  [931]  0.590359512 -0.129293335 -0.266279964 -0.084863154 -0.497895546
##  [936] -0.208414796 -0.179795381 -0.488649724 -0.326025892 -0.647595949
##  [941]  0.401677268  0.297846941 -0.445163655 -0.023398125  0.307499068
##  [946] -1.460328675 -0.405277899 -1.028813200 -0.603198197  0.208474304
##  [951] -0.165311153 -0.441659486 -1.793567627  0.196690068 -0.095434334
##  [956]  0.427639770  0.199007952 -0.031890517 -0.536379665  0.331295303
##  [961] -1.610527154 -0.027395645  0.467139350 -0.986550218  0.810420860
##  [966] -0.140027682 -1.458720280  1.053467059  0.607296372 -0.659671576
##  [971]  0.434959765  0.069292093  0.151113271 -0.539244895  0.276475533
##  [976] -0.208964595 -0.349149068  0.416296754 -0.776305224 -0.959082038
##  [981]  0.369700819 -1.040053492  0.802627630 -0.118777954 -0.015343966
##  [986] -1.647282326  0.102075634 -0.027621948  0.043390979 -0.049010240
##  [991]  0.295256995 -0.814600760 -1.239365181 -0.823272800  0.251248052
##  [996] -0.427751146  0.323938443 -0.486914443 -0.931650587  0.567518523
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
```

```r
fit2
```

```
##       mean          sd    
##   0.59573455   0.24877038 
##  (0.07866810) (0.05562236)
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
## [1] -0.4657030 -0.8825495 -0.3765631 -0.8048611 -0.7248213  0.7383707
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
## [1] 0.0186
```

```r
sd(results2$thetastar)  #the standard deviation of the theta stars is the SE of the statistic (in this case, the mean)
```

```
## [1] 0.9038534
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
## t1*      4.5 0.02252252   0.8963091
```

One of the main advantages to the 'boot' package over the 'bootstrap' package is the nicer formatting of the output.

Going back to our original code, lets see how we could reproduce all of these numbers:


```r
table(sample(x,size=length(x),replace=T))
```

```
## 
## 0 2 3 4 5 8 9 
## 1 1 1 2 1 2 2
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
## [1] -0.0316
```

```r
se.boot
```

```
## [1] 0.9047873
```

Why do our numbers not agree exactly with those of the boot package? This is because our estimates of bias and standard error are just estimates, and they carry with them their own uncertainties. That is one of the reasons we might bother doing jackknife-after-bootstrap.

The 'boot' package has a LOT of functionality. If we have time, we will come back to some of these more complex functions later in the semester as we cover topics like regression and glm.

