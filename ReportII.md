# Report II

Luhuan Wu, Xiaohui Li

June 4, 2017.   



## 1. Introduction

In this week, there are four things we did:

1. We revise the way we generate the in-degree and bi-degree sequences by generating the unifrom distribution in [0, 1].
2. We calculate the conditions for the parameters alpha, beta for the algorithm to converge.
3. We come up with two ways to evalute the correlation between page rank and betweenness centrality.  
   - plots of pagerank and betweennes centrality to get direct observation
   - Statistical measurement: rank correlation.
4. We set step-wise parameters to observe in which interval, page rank would overlap, or have positive correlation with, betweenness centrality, and in which interval it does not.



## 2. Procedures

### 1. Generating bi-sequences

$$
P(D^+ = d^+, D^- = d^-) = P(Poisson(W^+) = d^+, Poisson(W^-)=d^-)\\\text{where }(W^+,W^-) \text{ are jointly regularly varying.}\\P(W^+ > x) = (\frac{x}{b})^{-\alpha},\qquad x >b\\P(W^- > x) = (\frac{x}{c})^{-\beta}, \qquad x >c \\\qquad (\alpha, \beta > 1)\\W^+ = a(W^-)^d\\
z=F(x) =1-P(W^->x)=1-(\frac{x}{c})^{-\beta}\\
\Rightarrow x=F^{-1}(z)=c(1-z)^{-\frac{1}{\beta}}
$$

So setting
$$
U \sim \text{Uniform}(0, 1), W^{-1}=F^{-1}(U)=c(1-U)^{-\frac{1}{\beta}}
$$
We get the desired random variable.

### 2. Convergence conditions

In _ReportI_, we have conducted the relationships between the parameters by using the distribution and equal expectation of
$$
(W^+, W^-)
$$
as given below:
$$
d =  \frac{\beta}{\alpha}\\
\text{If } \alpha \neq \beta,  \qquad \qquad b = (\frac{\alpha (\beta - 1)}{(\alpha-1)\beta})^{\frac{\beta}{\alpha - \beta} }a^{\frac{\alpha}{\alpha - \beta}}\\
\qquad \qquad \qquad  c =  (\frac{\alpha (\beta - 1)}{(\alpha-1)\beta})^{\frac{\alpha}{\alpha - \beta} }a^{\frac{\alpha}{\alpha - \beta}} \\
\text{If } \alpha = \beta, \text{ then }a = 1 \text{ and we could arbitralily choose } b=c.	\qquad \qquad \qquad (3)
$$
In Algorithm 2.2, if we want 
$$
\frac{\Delta_n}{n} \rightarrow 0 \qquad \qquad \text{as }n \rightarrow \infty
$$
we must have 
$$
E[(W^+ - W-)^{1 + \epsilon}] = 0   \qquad \text{ for any given fixed } \epsilon >0
$$
We choose one sufficient condition:
$$
\epsilon =1 \text{ that is }Var(W^+ - W^-) = E[a(W^-)^d -W^-]^2  < \infty \\
\Rightarrow \beta >2, \beta >2d, \beta >d+1
$$

### 3. Evaluation of the correlation between page rank and betweeness centrality

#### 3.0 networkx implementation for page rank and betweenness centrality

[page rank_document](https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.algorithms.link_analysis.pagerank_alg.pagerank.html)

[betweenness centrality_document](https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.algorithms.centrality.betweenness_centrality.html?highlight=betweenness%20centrality#networkx.algorithms.centrality.betweenness_centrality):

> Betweenness centrality of a node ![v](https://networkx.github.io/documentation/networkx-1.10/_images/math/53c4a26799ff6cdae6f13d6b1fcf961660d83169.png) is the sum of the fraction of all-pairs shortest paths that pass through ![v](https://networkx.github.io/documentation/networkx-1.10/_images/math/53c4a26799ff6cdae6f13d6b1fcf961660d83169.png):
> $$
> c_B(v) = \sum_{s,t \in V} \frac{\sigma(s,t|v)}{\sigma(s,t)}
> $$
> where ![V](https://networkx.github.io/documentation/networkx-1.10/_images/math/ea879891318d7d3a04f0d3a7b0ad892a9b63c8e8.png) is the set of nodes, ![\sigma(s, t)](https://networkx.github.io/documentation/networkx-1.10/_images/math/f724eb676c09147ce1219ad0322502ba9136ce5a.png) is the number of shortest ![(s, t)](https://networkx.github.io/documentation/networkx-1.10/_images/math/a89084d7e1d70e18f88801fc46310063a5410644.png)-paths, and ![\sigma(s, t|v)](https://networkx.github.io/documentation/networkx-1.10/_images/math/f8be1c463f76d2bf1d7c06ea64d9ae7877f657ce.png) is the number of those paths passing through some node ![v](https://networkx.github.io/documentation/networkx-1.10/_images/math/53c4a26799ff6cdae6f13d6b1fcf961660d83169.png) other than ![s, t](https://networkx.github.io/documentation/networkx-1.10/_images/math/11ea1225465563f6ec492f643e46328d31d8c3f0.png). If ![s = t](https://networkx.github.io/documentation/networkx-1.10/_images/math/8ca592c5a04ef883ff9b192f44aa5edda9da03de.png), ![\sigma(s, t) = 1](https://networkx.github.io/documentation/networkx-1.10/_images/math/cc560c61381cf8c03e64b279f8bde36e3df96352.png), and if ![v \in {s, t}](https://networkx.github.io/documentation/networkx-1.10/_images/math/256407d40e389f52991c0d7b415e3487ee9dc6cd.png), ![\sigma(s, t|v) = 0](https://networkx.github.io/documentation/networkx-1.10/_images/math/c7d3eb38e416b28a38660834405a9c99d97ad661.png)

#### 3.1 Plots

We define two ways to plot the nodes to help us observe the correlation.

1. We plot the page rank of nodes in decreasing order, and then plot those nodes' betweenness centrality in order.
2. We plot the betweenness centrality of nodes in decreasing order, and then plot those nodes' page rank in order.

#### 3.2 Statistical measures 

##### 3.2.1 Rank correlation — Spearsman's test

>In [statistics](https://en.wikipedia.org/wiki/Statistics), a **rank correlation** is any of several statistics that measure an **ordinal association**—the relationship between [rankings](https://en.wikipedia.org/wiki/Ranking) of different [ordinal](https://en.wikipedia.org/wiki/Ordinal_data) variables or different rankings of the same variable, where a "ranking" is the assignment of the labels "first", "second", "third", etc. to different observations of a particular variable. A **rank correlation coefficient** measures the degree of similarity between two rankings, and can be used to assess the [significance](https://en.wikipedia.org/wiki/Statistical_significance) of the relation between them.
>
> —- <cite>[Wikipedia_Rank correlation][1]</cite>

We use [Spearsman's test][2] to measure the rank correlation between page rank and betweenness centrality. It basically measures  to the distance between the nodes' page-rank ranking and betweenness-centrality ranking.

Spearsman's test of two ranked data sequence evaluated by different ranking algorithms return two values: _correlation_ and _p-value_, while the corresponding null hypothesis is that two ranked data are statistically independent, which means two ranking algorithm has no correlation. 

The smaller the p-value, it's more likely to reject null hypothesis. 

A correlation of +1 or −1 occurs when each of the variables is a perfect monotone function of the other.

So we feed the page rank of nodes and betweenness centrality of nodes to Spearsman's test and observe the correlation and p-value to judge the correlation between two page rank and betweenness centrality.

##### 3.2.2 Measuring overlapping nodes

Given the graph size n, we calculate the overlapping number of nodes, m, of first k page-rank-ranked nodes and first betweenness-centrality-ranked nodes, and thereby get the _**percentage of overlapping nodes**_ p = m/k. That is
$$
p = \frac{\text{# nodes in both top k page-rank-ranked and top k b-c-ranked nodes}}{k}
$$

### 3. Set up parameters intervals to generate the models

According to the relationship of the parameters, we actually only need to set up 3 parameters to get the total 6. So we set the parameters intervals as follows:

```python
a_range = np.arange(0.6, 1.6, 0.2) # 5 choices: 0.6, 0.8, 1.0, 1.2, 1.4
beta_range = np.arange(2.5, 3.5, 0.2)  # 4 choices: 2.5, 2.7, 2.9, 3.1
d_range = np.arange(0.8, 1.3, 0.2)  # 3 choices: 0.8, 1.0, 1.2
```



## 3. Results and Analysis

### 3.1 Statistical test  

#### 1) K-S test

the **Kolmogorov–Smirnov test** (**K–S test** or **KS test**) is a nonparametric test of the equality of one-dimensional probability distributions that can be used  to compare two samples (two-sample K–S test).

#### 2) Wilcoxon signed-rank test

The **Wilcoxon signed-rank test** is a non-parametric hypothesis test used when comparing two related samples, matched samples, or repeated measurements on a single sample to assess whether their population mean ranks differ. 

### 3.1 Comparison of degree sequence

 We use the above statistical test to check whether the generated simple graph have asumptotically the same distribution as previous bi-degree sequence distribution.

#### a. In-degree



#### b. Out-degree



We could find that both tests' p-value in both in-degree and out-degree are not small enough to reject the null hypothesis, this gives sense to receive the alternative that the two samples (graph bi-degree sequence and original bi-degree sequence) are generated from the same distribution.



### 3.2 Correlation between page rank and betweenness centrality.

#### a. Full information in [result.txt][result]

#### b. Test for spearman correlation

**Spearman's rank correlation coefficient *ρ*** is a nonparametric measure of rank correlation statistical dependence between the ranking of two variables. It assesses how well the relationship between two variables can be described using a **monotonic function**.

Under the null hypothesis of statistical independence (*ρ* = 0), we could caculate the p-value to check whether *ρ* is statistically significant. If the p-value is small enough, we could reject the null and receive the dependence between the two ranking.

We use different pairs of alpha and d to generate simple directed graphs and test their pagerank and betweenness centrality statistical dependence.

## 4. Conclusion

#### a. Appropriate algorithm

Given WLLN's requirement for finite second moment constraint the parameters, we could generate the sample bi-degree sequence following the distribution at a fast enough speed. The simple directed configuration model generated from erased algorithm has the asympotically same degree distribution as the original bi-degree sequence degree distribution.



#### b. Betweeness centrality and page rank

Betweeness centrality and page rank ranking is statistically dependent.

- When the mean degree of the nodes get higher, the overlapping percentage tends to be higher.

  ​

-  When average degrees of nodes are bigger than, say around 2, 



[1]: https://en.wikipedia.org/wiki/Rank_correlation
[2]: https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
[result]: https://github.com/leahwu/Go_graphs/blob/master/results1.txt

####  


​			
​		
​	





