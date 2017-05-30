#Report I.  

Xiaohui Li, Luhuan Wu  
May 29, 2017


##1.Introduction:

The task is to generate a simple directed configuration model given in-degree and out-degree distribution. 


##2.Algorithms:  
The algorithms below except 2.2.3 are from [Directed Random Graph with Given Degree Distributions][dcm_algorithm]. 

2.1 Algorithm to generate bi-degree sequence that satisfies:
$$
\sum D^+ = \sum D^-
$$
2.2 The algorithms to transfrom the multigraph to the simple graph:

​	2.2.1 Repeated algorithm: keep generating random dcm untill obtaining a simple graph.

​	2.2.2  Erased algorithm: remove self-loops and parallel edges of the dcm 

​	2.2.3 Revised repeated algorithm: we modified the repeated algorithm so that it stops adding pairing the edges that result in self-loops or parallel edges. But we are not sure whether this approach is reasonable.

##3.Procedures:
### 3.1 Anlaysis and calculations

First, we need to generate the bi-sequences that have the following powerlaw distribution:
$$
P(D^+ = d^+, D^- = d^-) = P(Poisson(W^+) = d^+, Poisson(W^-)=d^-)\\\text{where }(W^+,W^-) \text{ are jointly regularly varying.}\\P(W^+ > x) = (\frac{x}{b})^{-\alpha},\qquad x >b\\P(W^- > x) = (\frac{x}{c})^{-\beta}, \qquad x >c \\\qquad (\alpha, \beta > 1)\\W^+ = a(W^-)^d
$$
Second, we can calculate the relationship between the parameters from the joint and margin distribution:
$$
\begin{align}
P(W^+ >x) &= P(a(W^-)^d>x) \\
&=P(W^->(\frac{x}{a})^{\frac{1}{d}}) \\
&=x^{-\frac{\beta}{d}}(ac)^{\beta}\\
&= x^{-\alpha}b^{\alpha}
\end{align} \\
\Rightarrow d = \frac{\beta}{\alpha},\qquad c = (\frac{b}{a})^{\frac{\alpha}{\beta}} \qquad (1)
$$
Furthermore, according to the conditions of Algorithm 2.1 to converge, the expectation of F and G should be equal. Conditional on W+ and W-, it suffices to let the expecation of W+ and W- to be equal. Thus, we have
$$
p(W^+) = \frac{d}{dt} (1-(\frac{x}{b})^{-\alpha}) = \frac{a}{b}(\frac{x}{b})^{-\alpha - 1} \\
E(W^+) = \int_b^\infty xp(W^+) dx = \frac{\alpha b^\alpha}{\alpha-1}b^{1-\alpha} \\
E(W^+) = E(W^-) \Rightarrow c = \frac{\alpha}{\beta}(\frac{\beta - 1}{\alpha - 1})b \qquad (2)
$$
By (1) & (2), we can derive
$$
d =  \frac{\beta}{\alpha}\\
\text{If } \alpha \neq \beta,  \qquad \qquad b = (\frac{\alpha (\beta - 1)}{(\alpha-1)\beta})^{\frac{\beta}{\alpha - \beta} }a^{\frac{\alpha}{\alpha - \beta}}\\
\qquad \qquad \qquad  c =  (\frac{\alpha (\beta - 1)}{(\alpha-1)\beta})^{\frac{\alpha}{\alpha - \beta} }a^{\frac{\alpha}{\alpha - \beta}} \\
\text{If } \alpha = \beta, \text{ then }a = 1 \text{ and we could arbitralily choose } b=c.
$$

So we only need to choose a, alpha and beta to initialize the model.

### 3.2 Progamming

We use python as our programming language and our programs develop upon the support of package [_networkx_](networkx) and its method [_directed\_configuration\_model_](dcm) specifically. We create 7 python files ([visit code repo][Go_graphs]) :

- PowerLawDistribution.py   

  To generate the specific bi-sequences that have  the specific powerlaw distribution.


- ValidDegree.py
  Input any bi-degree sequence, use algorithm 2.1 to make the sum of in-degree sequence equal to the sum of out-degree sequence.

- SDGErased.py. 

  Input any bi-degree sequence, use algorithms 2.2.2, the erased algorithm, which removes the self-roops and parallel edges to get simple directed graph. This module also provides methods to plot the degree distribution and the graph, etc.  

- SDGRepeated.py

  Input any bi-degree sequence, use algorithm 2.2.1, the repeated algorithm, which repeats generating the DCM untill obtaining a simple directed graph. This module also provides methods to plot.

  In actual implementation, we keep generating directed comfiguration model by calling method `directed_configuraion_model_revised` in DCMRevised.py until  `flag == True`.

- DCMRevised.py

  We modified the source code of [_directed\_configuration\_model_](dcm) by chaning only the part of random pariing the 'stubs', which stops as long as there exist self-loops or parallel edges. The idea is to stop generating the graph once it is added a self-loop or a parallel edge.

  The random pairing in the original source code is:

  ```
  while in_stublist and out_stublist:
          source = out_stublist.pop()
          target = in_stublist.pop()
          G.add_edge(source,target)
  ```

  We moditied this part to be:

  ```
   while in_stublist and out_stublist:
          source = out_stublist.pop()
          target = in_stublist.pop()
          if source == target or G.edges().__contains__((source, target)):
              flag = False
              print(source, target)
              break
          G.add_edge(source,target)
  ```

- DCMRevised2.py

  Again, we modified the random pairing part in the source code as mentioned above. The idea is to seek for new pairing once the graph is added a self-loop of a parallel edge.

  ```
      while in_stublist and out_stublist:
          source = out_stublist.pop()
          target = in_stublist.pop()
          if source == target or G.edges().__contains__((source, target)):
              print((source, target))
              continue
          G.add_edge(source,target)
  ```

- test1.py & test2.py

  Feed the parameters a, b, alpha, beta to the file and generate the desired results including the directed graph model and relevent plots.  


## 4. Results and Analysis

### 4.1 Repeated Algorithm

​	After running programs for 1 hour, we still fail to obtain the desired simple graph.

### 4.2 Erased Algorithm

### 4.3 Revised Repeated Algorithm





[dcm_algorithm]: https://arxiv.org/pdf/1207.2475.pdf
[networkx]: https://networkx.github.io/documentation/networkx-1.10/overview.html
[dcm]: https://networkx.github.io/documentation/networkx-1.10/_modules/networkx/generators/degree_seq.html#directed_configuration_model
[Go_graphs]: https://github.com/leahwu/Go_graphs