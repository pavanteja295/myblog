---
title: Iterative Closest point
draft: false
tags:
---




*Note*: I skipped the math as the [assignment](https://github.com/alecjacobson/geometry-processing-registration) already well explained the internals. My implementation is uploaded [here](https://github.com/pavanteja295/geometry-processing-csc2520-solutions/tree/main/geometry-processing-registration)
## Problem of registration 

Given two sets of points   $X, Y \subset  \mathbb{R}^3$ namely source and target respectively the problem of registration aims to find a transformation function $T: \mathbb{R}^3 \to \mathbb{R}^3$ such that  the two sets of points $T(X) = \{ T(x) |  x  \in X \}$ and $Y$ are *indistinguishable* as per some *defined*  distance measures $d$.

$$ {\arg\min}_{T \in \mathcal{T}} \quad  d(T(X), Y)  $$

In computer vision, registration plays a crucial role in combining scans from multiple depth cameras that capture depth information of an object from unknown but different viewpoints. This fusion enables to form a unified 3D representation of the object, which is essential for applications such as 3D modeling, robotics, and augmented reality.  In this application $X, Y$ correspond to the point sets from the two view points and $\mathcal{T}$ is the space of **rigid transformations** where $T(\mathbf{x}) = R\mathbf{x} + t$

In this assignment we are interested to solve the problem of registration under rigid transformations.  


##  Distance measures 

The *distance* measure $d$ of the two point sets $Z, Y \subseteq \mathbb{R}^3$ in a Euclidean space can be defined using the Euclidean norm (  $l_{2}$ ) as one of the following:

### Hausdorff distance


$$d_{H}(Z, Y)  =  \sup_{z \in Z} \; \inf_{y \in Y} \quad l_{2}(z, y ) $$

It measures the worst correspondence between $Z$ and $Y$.  $d_{H}$  is a good indicator function as  $d_{H}(Z, Y) = 0$ , then the all the points in $X$  found a (not necessarily a unique) correspondence in $Y$ such that $z=y$.   However, the main drawback under practical settings is  its sensitivity to outliers.  A more optimisation friendly is to consider average distance measure 


### Average distance 

$$d(Z, Y) = \int_{z \in Z} \; \inf_{y \in Y} \quad l_{2}(\mathbf{z}, y )$$

Assuming $Y$ to be finite set of closed triangles (as per the assignment)  (a point cloud) or **compact set** ( finite set of triangles with closed boundaries as given in the assignment) for function $l_{2}$  (continuous over set $Y$) there should exists a minima over $Y$ . Let this minima be denoted by $P_{Y}(\mathbf{z}) = \arg\min_{y \in \mathbf{Y}}  \; l_2(\mathbf{{z} , \mathbf{y}})$  and call this projection of the vector $\mathbf{z}$. $$d(Z, Y) = \int_{z \in Z} l_2(\mathbf{z}, P_Y(\mathbf{z}))$$

### Discrete Optimisation problem

Further assuming the set $X$ is finite with the above defined distance measure we can write the problem of registration as the following minimisation problem restricting the space of transformations to rigid i.e  $\mathcal{T} = \{T: \mathbf{R}^3 \to \mathbf{R}^3 | \; T(\mathbf{x}) = R\mathbf{x} + \mathbf{t}\}$ where $R$ is a 3D rotational matrix and $\mathbf{t}$ is a translation vector

$$
\begin{align}
R^*, t ^* &= \underset{R \in SO_{3}, \; t \in \mathbb{R}^3}{\arg\min}  \;  \sum_{\mathbf{x} \in X} l_{2}\bigg(T(\mathbf{x}), \; P_{Y}(T(\mathbf{x}))\bigg) \\ \\
&= \underset{R \in SO_{3}, \; t \in \mathbb{R}^3}{\arg\min}  \; \sum_{\mathbf{x} \in X} \bigg|\bigg| \mathbf{z} - P_{Y}(\mathbf{z}) \bigg|\bigg|^2
\end{align}
$$


where $\mathbf{z} = T(\mathbf{x}) = R\mathbf{x} + t$

### Challenges and proposal

Two key problems with the above optimisation problem:
1. $P_{Y}(T(\mathbf{x}))$ does not have a closed form function as it requires iterating over the entire set $\mathbf{Y}$ to find the closest point projection. 
2. Non-linearity or not quadratic of the optimisation problem wrt to  $R, t$  due to  $P_Y(T(\mathbf{x}))$ 




To address the above problems ICP proposes:
1. ICP proposes Iterative optimisation by moving closer to the projection estimates from previous step  since  there is no closed form expression for the $P_{Y}(T(\mathbf{x}))$ 

The iterative form of ICP can be rewritten for step $t$ as :

$$
\begin{align}
R_{t+1}^*, t_{t+1} ^* &= \underset{R \in SO_{3}, \; t \in \mathbb{R}^3}{\arg\min}  \; \sum_{\mathbf{x} \in X} \bigg|\bigg| R\mathbf{z}_{t} + t - P_{Y}(\mathbf{R\mathbf{z}_{t} + t}) \bigg|\bigg|^2
\end{align}
$$

2. To make the optimisation problem  linear  $P_{Y}(T(\mathbf{{x}}))$  is chosen to be :
	1. **Point - point matching** - $P_{Y}(\mathbf{z}_{t+1})$ is constant across two different iterate steps.
	2. **Point - plane matching** -  $P_{Y}(\mathbf{z}_{t+1})$ is derived from $P_{Y}(\mathbf{z}_{t})$ assuming $Y$ is a planar surface.
 

### Point - point matching 

Point-point matching is generated and could be used for the case when $Y$ is also a point cloud (with no normal information). In case of point-point matching the assumption is that the projection of $\mathbf{z}_{t+1}$ and $\mathbf{z}_{t}$ is same. 

$$
\begin{align}
R_{t+1}^*, t_{t+1} ^* &= \underset{R \in SO_{3}, \; t \in \mathbb{R}^3}{\arg\min}  \; \sum_{\mathbf{x} \in X} \bigg|\bigg| R\mathbf{z}_{t} + t - P_{Y}(\mathbf{z_{t}}) \bigg|\bigg|^2
\end{align}
$$


### Point-plane matching
Point-plane matching is generated and could be used for the case when $Y$ also contains additional information like surface normals $\hat{n}$ . 

In this case we can extrapolate $P_{Y}(Z_{t +1})$ by assuming that the surface is locally planar and $Z_{t}$ and $Z_{t+1}$ are in close proximity and their projections lie on the same plane :

$$
P_{Y} (Z_{t+1})  =  Z_{t+1} +  \bigg \langle R\mathbf{x}_{t} + \mathbf{t} - P_{Y}(Z_{t}), \; \hat{n}_{P_{Y}(Z_{t})}  \bigg  \rangle \; \hat{n}_{P_{Y}(Z_{t})}
$$

The optimization problem in this case is: 

$$
\begin{align}
R_{t+1}^*, t_{t+1} ^* &= \underset{R \in SO_{3}, \; t \in \mathbb{R}^3}{\arg\min}  \; \sum_{\mathbf{x} \in X} \bigg|\bigg| \bigg \langle Z_{t+1} - P_{Y}(Z_{t}), \; \hat{n}_{P_{Y}(Z_{t})}  \bigg  \rangle \; \hat{n}_{P_{Y}(Z_{t})}\bigg|\bigg|^2
\end{align}
$$



### Ok but why does ICP algorithm work ?

At first glance I found ICP to be some wizardry. Although the optimisation problem in point-point ICP is the infamous least squares with fixed target (the usual case in supervised machine learning) at each optimisation step ICP finds new set of closest points $P_{Y}(z_{t})$ and optimize $R, t$ for these new projections.  To clarify this ICP lets put down the algorithm for ICP from 

 
$$
\begin{align*}\\ 
&\textbf{Algorithm: } \text{Iterative Closest Point}\\ 
&\textbf{Input: } \text{X, Y}\\ 
&\textbf{Output: } \text{R, t}\\
\\

& R, \mathbf{t} \Leftarrow \text{ random()} \\
& Z = X\\\\
&\textbf{While } \text{not converged}: \\
 & \quad \quad P_{Y}(Z_{t-1}) \leftarrow \text{ Step -1 : Determine the closest points to $Z_{t-1}$ on } Y\\
& \quad \quad \mathbf{R}_{t},  \mathbf{t}_{t} \leftarrow \text{ Step - 2  : Minimize } \bigg|\bigg| R\mathbf{Z}_{t-1} + t - P_{Y}(\mathbf{Z_{t-1}}) \bigg|\bigg|^2\\
& \quad \quad Z_{t} = R_{t}Z_{t-1} + t_{t}\\

\end{align*}
$$




1. Step - 1 for iteration  *t* : determines the closest points $P_Y(Z_{t-1})$ on $Y$ for all the points $Z_{t-1}$ such that distance between the two pairs is minimum.
2. **Step -2** for iteration  *t*:  Finds $R, t$ to minimize such that the distance between $P_{Y}(Z_{t-1})$ and $Z_{t}$ is less than that between $P_{Y}(Z_{t-1})$ and $Z_{t-1}$, 
3. Again, **Step - 1** for iteration  *t +1* :  determines the closest points $P_Y(Z_{t})$ on $Y$ for all the points $Z_{t}$ . Hence, the distance between $P_{Y}(Z_{t})$ and $Z_{t}$ is less than that between $P_{Y}(Z_{t-1})$ and $Z_{t}$, i.e 


$$ l_{2}(P_{Y}(Z_{{t-1}}), Z_{t-1}) \quad \leq \quad  l_{2}(P_{Y}(Z_{{t-1}}), Z_{t}) \quad \leq \quad l_{2}(P_{Y}(Z_{{t}}), Z_{t}) \quad \ldots$$


Hence, in the limit the distance between $Z_{t}$ and $P_{Y}(Z_t)$ **always** converge, hopefully to zero. 

### Does ICP always work for registration? 

1. **Misalignment due to sparsity and noise**: 
	1. If the point cloud sample $X$ is sparse then the final correspondence could be to some erroneous points in $Y$. 
	2. Similarly, if the point cloud sample $X$ is noisy then it could also lead to misalignment. (Could be interesting to do some experiments and analysis )