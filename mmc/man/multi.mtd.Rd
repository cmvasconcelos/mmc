\name{multi.mtd}
\alias{multi.mtd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Markov Chains
}
\description{
This function estimates Ching (2002) model for Multivariare Markov Chains with Berchtold (2001) optimization algorithm.
}
\usage{
multi.mtd(y, deltaStop=0.0001, is_constrained=TRUE, delta=0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Matrix of categorical data sequences.
}
\item{deltaStop}{Value of stop criteria.}
\item{is_constrained}{Logical value, indicating if the estimation is constrained.}
\item{delta}{Quantity to increase/decrease in each parameter for each iteration of the algorithm.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Ching, W. K., E. S. Fung, and M. K. Ng (2002). 'A multivariate Markov chain modelfor categorical data sequences and its applications in demand predictions'. In: IMAJournal of Management Mathematics13(3), pp. 187-199.DOI:10.1093/imaman/13.3.187.


Berchtold,  A.  (2001).  'Estimation  in  the  Mixture  Transition  Distribution  Model'.  In:Journal of Time Series Analysis22(4), pp. 379-397.DOI:https://doi.org/10.1111/1467-9892.00231
}
\author{
Carolina Vasconcelos
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
%\examples{
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
