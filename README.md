# RMTCovEst
Implement the improved covariance estimation for large class of metrics from paper
'Random Matrix Improved Covariance Estimation for a Large Class of Metrics, ICML 2019'

# 2 versions of the code: 
Version 1 Implementation of the algorithm 2 for linear shrinkage initialization--> Faster implementation

Version 2 Implementation of the algorithm 1 for general initialization
# RMTCovEst
Implement the improved covariance estimation for large class of metrics

# 2 versions of the code: 
Version 1 Implementation of the algorithm 2 for linear shrinkage initialization--> Faster implementation

Version 2 Implementation of the algorithm 1 for genral initialization

# Archive content
1- The function implementing our method is called \texttt{RMTest.m} which takes as arguments the data matrix \texttt{X}, the \texttt{gradient\_check} and \texttt{plot\_cost} options, the initial point \texttt{C\_0} and the distance under investigation. The function returns the estimated covariance matrix \texttt{C\_est} and the cost function related to this estimated covariance.

2-  The main script comparing all algorithms for synthetic data is \texttt{CompareEst.m}.

3- The main script comparing the estimation methods for the LDA/QDA application is \texttt{ML\_applications.m}.

4-  folder \texttt{Manopt}: a Matlab toolbox for optimization on manifolds \cite{manopt}. Some of the functions in this folder were adapted to better suit our present problem.

5-  folder \texttt{Othermethods}: containing alternative estimation methods among which the QuEST methods QuEST1 and QuEST2 from \cite{LW15,LW18}, the Rao-Blackwell Ledoit-Wolf estimation methods \cite{chen2010shrinkage} and the Oracle Approximating Shrinkage estimation method \cite{chen2010shrinkage}.

6-  folder \texttt{Utilities}: contains supplementary codes for the implementation of LDA,QDA.

# Code CompareEst.m
The different options proposed to execute the script \texttt{CompareEst.m} comparing the different estimation algorithms are as follows:

1- The range of $n$ and the value of $p$

2- The covariance matrices metric ``\texttt{distance}'' (among \texttt{Fisher, Battacharrya, KL, log, log1st, t}) or for the precision matrices (among \texttt{Inverse\_Fisher,
    ,Inverse\_Battacharrya,Inverse\_KL,} \texttt{Inverse\_log1st,
    Inverse\_log,Inverse\_t})
    
# Archive content
1- The function implementing our method is called \texttt{RMTest.m} which takes as arguments the data matrix \texttt{X}, the \texttt{gradient\_check} and \texttt{plot\_cost} options, the initial point \texttt{C\_0} and the distance under investigation. The function returns the estimated covariance matrix \texttt{C\_est} and the cost function related to this estimated covariance.

2-  The main script comparing all algorithms for synthetic data is \texttt{CompareEst.m}.

3- The main script comparing the estimation methods for the LDA/QDA application is \texttt{ML\_applications.m}.

4-  folder \texttt{Manopt}: a Matlab toolbox for optimization on manifolds \cite{manopt}. Some of the functions in this folder were adapted to better suit our present problem.

5-  folder \texttt{Othermethods}: containing alternative estimation methods among which the QuEST methods QuEST1 and QuEST2 from \cite{LW15,LW18}, the Rao-Blackwell Ledoit-Wolf estimation methods \cite{chen2010shrinkage} and the Oracle Approximating Shrinkage estimation method \cite{chen2010shrinkage}.

6-  folder \texttt{Utilities}: contains supplementary codes for the implementation of LDA,QDA.

# Code CompareEst.m
The different options proposed to execute the script \texttt{CompareEst.m} comparing the different estimation algorithms are as follows:

1- The range of $n$ and the value of $p$

2- The covariance matrices metric ``\texttt{distance}'' (among \texttt{Fisher, Battacharrya, KL, log, log1st, t}) or for the precision matrices (among \texttt{Inverse\_Fisher,
    ,Inverse\_Battacharrya,Inverse\_KL,} \texttt{Inverse\_log1st,
    Inverse\_log,Inverse\_t})
    
3- The target matrix ``\texttt{Covariance}'' (among \texttt{dirac,Wishart, toeplitz}) and their parameters ``\texttt{param}'' if needed (for the Wishart and Toeplitz cases) 

4- The initialization point for the gradient descent algorithm denoted ``\texttt{initialization}'' (linear shrinkage from \cite{LED04} denoted ``\texttt{shrinkage}'', shrinkage from \cite{chen2010shrinkage} denoted ``\texttt{alternative shrinkage}'', QuEST denoted as ``\texttt{ledoit-wolf}'') or the identity denoted ``\texttt{manual}''

5-  Other binary option can be chosen (among 0/1): (\texttt{gradient\_check} to check if the gradient is correct, \texttt{plot\_cost} to see the cost/real distance during iterations.
Code \texttt{ML\_applications.m}
The different options proposed to execute the script \texttt{ML\_applications.m} are as follows:

1-  the data on which the LDA/QDA are applied denoted ``\texttt{dataset}'' for which the options are \texttt{synthetic} for synthetic data and \texttt{eeg} for eeg dataset.

2- The machine learning algorithms denoted ``\texttt{application}'' for which the options are \texttt{LDA} and \texttt{QDA}.

3-  For synthetic data, the examples of covariance under investigation. For both the covariance of the first and second class, \texttt{Wishart} and \texttt{toeplitz} are the two options.


# Reproducing the results of the article

The following sections detail the parameter setting to reproduce the figures of the main article.

## Figure 1
Script $\rightarrow$ \texttt{CompareEst.m}\\
covariance $\rightarrow$ toeplitz\\
parameter $\rightarrow$ 0.4\\
p $\rightarrow$ 512, n $\rightarrow$ 500\\
distance $\rightarrow$ Fisher\\
Initialisation $\rightarrow$ manual\\
plot\_cost $\rightarrow$ 1


## Figure 2

Script $\rightarrow$ \texttt{CompareEst.m}

covariance $\rightarrow$ Wishart/toeplitz(0.1)/toeplitz(0.9)/dirac

n $\rightarrow$ linspace(210,500,10), p $\rightarrow$ 200

distance $\rightarrow$ Fisher

Initialisation $\rightarrow$ shrinkage

plot\_cost $\rightarrow$ 0


## Figure 3
Script $\rightarrow$ \texttt{CompareEst.m}

covariance $\rightarrow$ Wishart/toeplitz(0.1)/toeplitz(0.9)/dirac

n $\rightarrow$ linspace(210,500,10), p $\rightarrow$ 200

distance $\rightarrow$ Inverse\_Fisher

Initialisation $\rightarrow$ shrinkage

plot\_cost $\rightarrow$ 0


## Figure 4
Script $\rightarrow$ \texttt{ML\_applications.m}

mu2 $\rightarrow$ mu1+80/p

application $\rightarrow$ LDA

covariance1 $\rightarrow$ Wishart/Wishart/toeplitz/

covariance2 $\rightarrow$ Wishart/toeplitz/toeplitz/
dataset $\rightarrow$ synthetic/synthetic/synthetic/eeg

p $\rightarrow$ linspace(500,200,10), n $\rightarrow$ 512
eeg dataset $\rightarrow$ $n\_train=500$, $n\_test=1000$,p=100


## Figure 5 
Script $\rightarrow$ \texttt{ML\_applications.m}

mu2 $\rightarrow$ mu1+1/p except for the third one which is mu1+80/p

application $\rightarrow$ QDA

covariance1 $\rightarrow$ Wishart/Wishart/toeplitz/

covariance2 $\rightarrow$ Wishart/toeplitz/toeplitz/

dataset $\rightarrow$ synthetic/synthetic/synthetic/eeg

p $\rightarrow$ linspace(500,200,10), n $\rightarrow$ 512

eeg dataset $\rightarrow$ $n\_train=150$, $n\_test=1000$,p=100

Please contact tiomoko_malik@yahoo.fr for further details
