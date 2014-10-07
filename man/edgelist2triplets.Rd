\name{edgelist2triplets}
\alias{edgelist2triplets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identfication of the regulatory triplets given a gene regulatory network edge list 
}
\description{
This function obtains the regulatory triplets (RF1, RF2, target) from the input gene regulatory network, grn. 
}
\usage{
edgelist2triplets(grn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grn}{
Edge list of gene regulatory network, grn: 1st column is regulatory factor (RF), and 2nd column is target. 
}
}
\details{
Input the edge list of gene regulatory network, and output the regulatory triplets of the network.
}
\value{
The output table contains the regulatory triplets of gene regulatory network
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Daifeng Wang, Mark Gerstein
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
