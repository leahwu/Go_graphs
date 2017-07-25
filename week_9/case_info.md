### Case I

alpha = beta = 2.1, E=3, d=0,0.25,0.5,0.75,1.0

GiantComponent:  3720
GiantComponent:  3802
GiantComponent:  3928
GiantComponent:  3903
GiantComponent:  3883

sample corr:

0.0097,	0.2368,	0.3835,	0.627,	0.9806

computed degree corr(theoretical):

0.0,		0.2899,	0.6402,	0.8687,	0.9345

对应的图是1，2，3，4，5

以及corr1 （degree corr随d的变化）

观察： 基本上total degree的影响大于betweeness centrality（差不多这两个）， 大于out degree大于in degree 大于page rank

但是特例：d=0的时候outdegree 的影响小于in degree

当d增大的时候， 五条线越来越接近，意味着这五个指标相关性越来越强

### Case 2

alpha=2.1, beta=2.3, E=3, d=0,0.25,0.5,0.75,1.0

GiantComponent:  3869
GiantComponent:  3886
GiantComponent:  3928
GiantComponent:  3943
GiantComponent:  3953

sample corr:

0.0228,	0.1512,	0.5584,	0.6563,	0.7122

computed degree corr(theoretical):

0.0,		0.2326,	0.5286,	0.7181,	0.7717

对应的图是1，2，3，4，5

以及corr1 （degree corr随d的变化）

观察结果如上

