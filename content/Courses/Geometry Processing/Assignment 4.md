---
title: Smoothing by Laplacian
draft: false
tags:
---
****
Checkout this [link](https://github.com/pavanteja295/geometry-processing-csc2520-solutions/tree/main/geometry-processing-smoothing) for implementation of [assignment 4 ](https://github.com/alecjacobson/geometry-processing-smoothing).  I felt the smoothing formulation in the assignment is a bit confusing fiddling between  the functional form  (energy term) as well as the weak formulation. In this notes I tried to stick to a single view and derived the final update rules.

## Function smoothing on surfaces

1.  Informally, given a scalar function $f:\Omega \to \mathbb{R}$ over some surface $\Omega$ that is *noisy* (that is the function $f$ *fluctuates* quite a bit around the neighbourhood of any point $x \in \Omega$) we would like to find another function $u^{\prime}$  that is close approximation to $u$  but less *noisy* (that is the the function $u^{\prime}$ at some $x\in \Omega$ deviates less from  neighbhourhood average) . 
2. In image processing out of many criteria for smoothness a common operators to smooth a 2D image is a convolution with Gaussian filter. In this assignment we discuss an alternative operator to smoothing a function using Laplace operator.


## Laplacian operator 
1. Laplacian operator for $f \in C^2$ is a second order differential operator of form: 

$$
\Delta f  = \frac{\partial^2 f }{\partial x^2} + \frac{\partial^2 f }{\partial y^2} + \frac{\partial^2 f }{\partial z^2}
$$
1. Below I state a simple proposition that sheds some light on why Laplacian is useful for the smoothing problem. Informally,  Laplacian captures deviations of the function value from local neighbourhood averages


----

![[./prep1.png]]
![[content/Courses/Geometry Processing/prep1.png]]

![[public/static/prep1.png]]
![[../../../public/static/prep1.png]]

***Proof***
1. Consider the Taylor series approximation of $f$  at $x$, then for $\bar{x} \in S_{R}(x)$ we have :

$$
f(\bar{x}) = f(x) + (x - \bar{x})f^{\prime}(x) + \frac{(x - \bar{x})^2f^{\prime\prime}(x)}{2} + o(||x-x||^3)
$$
2. Taking the integral over $S_{R}(x)$ the middle term vanishes (points on either side of the diameter cancel each other) . Then we have 
$$
\begin{align}
\int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) &=  \int_{\bar{{x}} \in S_{R}(x)} \frac{(x - \bar{x})^2f^{\prime\prime}(x)}{2} +  o(||x-x||^3) \\
\int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) -  f(x) &= + \int_{\bar{{x}} \in S_{R}(x)} \frac{(x - \bar{x})^2f^{\prime\prime}(x)}{2} +  o(||x-x||^3) \\ \\
\int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) -  f(x) &= + \int_{\bar{{x}} \in S_{R}(x)} \frac{R^2f^{\prime\prime}(x)}{2} +  o(||x-x||^3) \\ \\
\int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) -  f(x) &= A_{n-1} R^{n-1} \frac{R^2f^{\prime\prime}(x)}{2} +  o(R^3) \\ \\
\frac{2}{A_{n-1} R^{n + {1}}} \int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) -  f(x)  + o(R^3)&= f^{\prime\prime}(x)  \\
\end{align}
$$
3. Applying the limit $R \to 0$ yields the result in the proposition. 
$$
\boxed{\Delta f = f^{\prime\prime}(x) = \lim_{R \to 0  } \;  \frac{2}{A_{n-1} R^{n + {1}}} \int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) -  f(x)  + \underbrace{  o(R^3)}_{\to 0} }
$$
4. Hence, this proves Laplacian capture the  deviations of a function  value from the local average of the shell.  The result needs to be extended to the case of multi variate functions but the spirit should remain the same.
5. The function $f$ in addition to being twice differentiable should also be Lesbegue integrable. 

---

3. Note that Laplacian operator is defined on functions over the Euclidean spaces. For Laplacian operator for functions on surfaces check out  Laplace-Beltrami operator [link](https://en.wikipedia.org/wiki/Laplace%E2%80%93Beltrami_operator)  . For this assignment we don't even use the Laplacian operator we only use the gradients of the defined step function on the mesh by exploiting the weak formulation.

### Laplacian as smoothing operator
1. Essentially in order to smooth a function  $f$ we can update the function to $f^{\prime}$ as :

$$
f^{\prime} = f+ \lambda \Delta f
$$
where $\lambda$ controls the rate of smoothening.  If we write the above update as a differential equation to smooth the function multiple times gives us a *weighted* version of the infamous [heat equation](https://en.wikipedia.org/wiki/Heat_equation)

$$
\boxed{\frac{\partial f}{\partial t} = \lambda \Delta f}
$$
^05ce41


Below we discuss two numerical methods to solve the above partial differential equation. 
 
# Finite Difference method (FDM) : 

1. I have been thinking along the lines of gradient descent to solve the above PDE. Given the update direction one could naively determine  $f_{t+1}$ by first determining Laplacian of the function $\Delta f_{t}$  at every point and update $f_{t}$ as :

$$
f_{t+1} = f_{t} + \lambda \Delta f_{t}
$$

2. Two step method for finite difference: 
	1. Discretise the Laplacian  on the given triangle mesh- See [Laplace Beltrami operator](https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Mesh_Laplacians)
	2. Determine $f_{t+1}$ by the above update rule (either using explicit or implicit Euler scheme)
	
3. However FDM seems to not be robust for complex geometries.


> [!todo|green] FDM vs FEM
> - [ ] Why is FEM efficient than FDM , mathematically? - Boundary conditions
> - [ ] Also come up with certain instances where FDM fails and FEM thrives - Complex geometries
> - [ ] This is an interesting and seems to be a fundamental question to the origins of FEM - See [this](https://scicomp.stackexchange.com/questions/290/what-are-criteria-to-choose-between-finite-differences-and-finite-elements) and [this ](https://en.wikipedia.org/wiki/Finite_element_method#Comparison_to_the_finite_difference_method)



#  Finite element method (FEM) : 

1. FEM seems to be like the most ubiquitous method used to solve PDEs for complex geometries. 
2. Here I only outline a crude formulation of FEM while leaving the crucial questions on its stability for a later point.

## Formulation 

1. Consider a continuous function $f(x , t)$ , $x\in \Omega, t\in \mathbb{R}$ that we would like to numerically solve for given some boundary conditions and partial differential equation constraints.  
2. Numerical solution of the PDE with its boundary constraints is determined by discretisation of the function $f$ wrt to the space variables $x$ and time variables $t$ 

3. Two key ideas behind finite element method : 
	1. Weak Formulation of PDEs
	2. Discretisation over space variables
	
	Consider the heat equation.   See [[#^05ce41]] .
	
## Weak solutions

1. Weak formulation solves PDEs by *weakening* the point-wise constraints defined over the entire $\Omega$  using integrals.  This leads to *weak* solutions where the solution satisfies the constraint with respect to some chosen test functions $v$.  
2. 
	1. **Strong solution** : If $f^*$ is a strong solution for the original heat equation then it requires the solution to hold the point-wise equality: 

$$
\forall x \in \Omega  ; \quad \frac{\partial f^*}{\partial t}(x) \;\; = \;\;  \lambda \Delta f(x)   \;
$$
	2. **Weak solution** : If $f^*_{w}$ is a weak solution then it requires that the $f_{w}^*$ must hold the following weak implied constraints with some test functions $v_{i} \in V$  where $V$ is a function space of choice.
$$
\forall v_{i} \in V \quad  \int_{x \in \Omega} \quad \frac{\partial f^*}{\partial t}(x) v_{i} = \int_{x \in \Omega} \lambda \Delta f(x)v_{i} 
$$

^weaksolutionheat

   The above [[#^weaksolutionheat]] can be reduced using Green's identity and integrating by parts for each test function $v_{i}$ 
   
$$
\begin{align}
 \quad  \int_{x \in \Omega} \quad \frac{\partial f^*}{\partial t}(x) v_{i} &= -\int_{x \in \Omega} \lambda \nabla  f(x) \nabla v_{i}  +  \oint_{\partial \Omega} \langle{\nabla f} , \;  n \rangle\; v_{i}  \\ \\
&= -\int_{x \in \Omega} \lambda \nabla  f(x) \nabla v_{i} + \cancel{\oint_{\partial \Omega} \langle{\nabla f} , \;  n \rangle\; v_{i} }_{=0} \\ \\
&= -\int_{x \in \Omega} \lambda \nabla  f(x) \nabla v_{i} 
\end{align}
$$
^weakformheat

The boundary integral is zero for closed meshes and if the mesh is not closed we choose test functions $v_{i}$ that are zero on the boundary edges of the mesh. This removes the surface integral term.


3.  Why is the weak formulation required  in the first place to solve for the PDE : 
	1. By careful selection of the function space of test functions restricted to some compact support functions the point-wise constraints over the entire space $\Omega$ are relaxed to integral forms and constraint is required to be satisfied distributionally over the support of the chosen test function.
	2. Theoretically,  by results like Lax-Milgram Theorem we know that the weak formulation has unique and existing solutions. 
	3. If there are additional boundary conditions given in the problem (although this is not the case for the given heat equation) by some wise choice of test functions the boundary conditions could be relaxed or totally ignored. 
	4. Under the integrals, using the properties like integral by parts, Green identity , costly second order operators can be reduced to simple first order operators. 


> [!todo|green] Weak formulation
> - [ ]  Formally, define the weak formulation. 
> - [ ]  With concrete examples showcase the advantage of using a weak formultion for finding the solutions.
> - [ ] Proof for Lax-Milgram theorem. 


### Discretisation 

1.  Given the weak form of heat equation (See [[#^weaksolutionheat]]) . The above weak formulation is defined over $f, v \in V$  , the functional space, that is infinite dimensional. 
2. One way to return the dimensionality of the problem is to restrict the function space $V$ to be a finite dimensional space. 
3. A common choice of $V$ is to be piece-wise *linear* functions with a linear function defined for every triangle element.
4. Usually, considering the mesh vertices  as $n$ control points $V$ could be defined as piece-wise linear functions using the control points as $v \in V$

$$
v(x) = \sum^n_{j=1} u_{j} \phi_{j}(x) 
$$
  with $\phi_{j}$ defined as:
$$
\phi_{j}(x) = \begin{cases} \\
1 \quad \text{if} \quad x \in \text{Mesh vertices} \\
\frac{A_{x}}{A} \quad \text{if} \quad x \in \Delta_{j}  \\ \\
0 \quad \text{otherwise} 
\end{cases}
$$
The above function space $V$ is parameterised by the value of the control points $\{u_{1}, u_{2} \ldots  u_{n}\}$.  Hence, with the specific choice of the basis functions the function space $V$ if finite dimensional.   Plugging in this assumption on the function space $V$ into [[#^weakformheat]]  and choosing $v_{1} = \phi_{1} , v_{2} = \phi_{2} , \ldots v_{n} = \phi_{n} ,$  gives :


$$
\begin{align}
-\int_{x \in \Omega} \lambda \nabla  f \nabla v_{i}  &= -\int_{x \in \Omega} \lambda \nabla  \bigg(  \sum^n_{j=1}f_{j} \phi_{j }\bigg)\nabla\phi_{i} \\ 
&= -\int_{x \in \Omega} \lambda   \bigg(  \sum^n_{j=1}f_{j} \nabla \phi_{j }\bigg)\nabla\phi_{i} \\ 
&= -\int_{x \in \Omega} \lambda   \bigg(  \sum^n_{j=1}f_{j} \nabla \phi_{j } \nabla\phi_{i}\bigg) \\
&= - \lambda \sum^n_{j=1}f_{j}  \underbrace{\bigg(   \int_{x \in \Omega}\nabla \phi_{j } \nabla\phi_{i}\bigg)}_{L_{ij}} \\ 
&= - \lambda \sum^n_{j=1}f_{j}L_{ij}

\end{align}
$$
 $L$  is defined to be *stiffness* matrix that approximates the integral. Note that the term $\nabla \phi_{j } \nabla\phi_{i}$   is compact and non-zero over the triangles containing the edge between vertex $i$ $vertex$ $j$.  The matrix L is independent on the function chosen in the function space $F$. 


Similarly, the left hand term can be reduced to : 

$$
\begin{align}
\int_{x \in \Omega} \quad \frac{\partial f}{\partial t} v_{i} &= \int_{x \in \Omega} \quad \frac{\partial }{\partial t}\bigg( \sum^n_{j=1} f_{j} \phi_{j} \bigg ) \phi_{i}  \\
&= \int_{x \in \Omega} \quad \sum^n_{j=1} \frac{\partial f_{j} }{\partial t}  \bigg( \phi_{j}  \phi_{i} \bigg) \\
&=  \quad \sum^n_{j=1} \frac{\partial f_{j} }{\partial t}  \underbrace{\bigg( \int_{x \in \Omega}\phi_{j}  \phi_{i} \bigg)}_{M_{ij}}\\ \\
&=  \quad \sum^n_{j=1} \frac{\partial f_{j} }{\partial t}  {M_{ij}}\\
\end{align}
$$

 $M$ is defined to be *mass* matrix.  Both $M$ and $L$ could be precomputed given the mesh geometry and is independent on the function $f \in V$.  This converts the infinite dimensional weak formulation to finite dimensional formulation as 

$$
\begin{align}
\forall v_{i} \in \{ v_{1} , v_{2} , \ldots , v_{n} \}\qquad\qquad & \sum^n_{j=1} \frac{\partial f_{j} }{\partial t}  {M_{ij}} =  \lambda \sum^n_{j=1}f_{j}L_{ij} \\  \\
&M \frac{\partial f }{\partial t}  = -\lambda L f \\ \\
&\boxed{\frac{\partial f }{\partial t}  = -\lambda M^{-1}L f }\\ \\
\end{align}
$$
^finalupdaterule




## Numerical solution of FEM

1. In summary , the PDE is reduced to a simple update rule using the stiffness matrix and mass matrix - $L, M$ that are computed from the mesh geometry (Checkout  [Assignment 4](https://github.com/alecjacobson/geometry-processing-smoothing) for details on computing the matrices). 
2. Explicit smoothing iteration  (forward euler method) using [[#^finalupdaterule]]  can be written as: 

$$
\begin{align}
f_{t+1} &= f_{t} +  \frac{\partial f_{t} }{\partial t}  \\
&= f_{t } - \lambda  M^{-1}L f_{t}  \\
\end{align}
$$


$$
\boxed{f_{t+1}= (Id - \lambda M^{-1} L)f_{t}} \\
$$

3. For more stable smoothing process, implicit smoothing iteration has been used (backward euler) in the [assignment](https://github.com/alecjacobson/geometry-processing-smoothing) (there is a typo in the sign of the update rule)

$$
\begin{align}
f_{t+1} - f_{t} = -\lambda M^{-1}Lf_{t+1} \\
\end{align}
$$
$$
\boxed{f_{t} = (Id   + \lambda M^{-1}L) f_{t+1}}
$$

> [!todo|green] Implicit vs Explicit methods
> - [ ] Difference between implicit and explicit
> - [ ]  When is one ideal than other
> - [ ] Concrete example for showing the difference between implicit vs explicit method 







> [!proposition]+ Laplacian capture deviations from local Average : Let $f:\mathbb{R} \to \mathbb{R}$ be twice differentiable function at $x$, $f \in C^2$ ,  then 
> $$
> \begin{align}
\Delta f(x) &= 
> \lim_{R\to 0 } \; \frac{2}{A_{n-1} R^{n + {1}}}\int_{\bar{{x}} \in S_{R}(x)} f(\bar{x}) - f(x)  \\
 &= \bar{f}_{S_{R}}(x) - \lim_{R \to 0 } \frac{2}{R^2} f(x)\\
> \end{align}
> $$
> where $S_{R}(x) = \{ \bar{x } \in R : || \bar{{x}} - x|| =R \}$ denotes the outer shell of a $R$-ball around $x$ , $A_{n-1}$ denotes the  hypervolume of the unit sphere. 
