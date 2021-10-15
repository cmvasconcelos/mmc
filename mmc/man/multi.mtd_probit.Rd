\name{multi.mtd_probit}
\alias{multi.mtd_probit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unconstrained Multivariate Markov Chains model
}
\description{
Estimation of Nicolau (2014) model for Multivariate Markov chains.
}
\usage{
multi.mtd_probit(x, initial, nummethod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix of categorical data sequences.
}
\item{initial}{Initial values}
\item{nummethod}{Optimization method. For more details see maxLik()}
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
Nicolau, J. (2014). 'A new model for multivariate markov chains'. In:Scandinavian Journal of Statistics41(4), pp. 1124-1135.ISSN: 14679469.DOI:10.1111/sjos.12087.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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
