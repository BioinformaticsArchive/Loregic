\name{loregic}
\alias{loregic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of logic gate consistency scores for all regulatory triplets of a gene regulatory network
}
\description{
This function calculates the consistency scores of 16 two-input-one-output logic gates for all regulatory triplets in inputted gene regulatory networks.
}
\usage{
loregic(grn, bdata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grn}{
Edge list of gene regulatory network, grn: 1st column is regulatory factor (RF), and 2nd column is target. Note that users can directly input triplets for grn with 1st column being RF1, 2nd column being RF2, and 3rd colun being target.
}
  \item{bdata}{
Binarized gene expression data matrix (genes by conditions) for genes in inputted gene regulatory network, grn.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The output table contains the regulatory triplets and their consistency scores across 16 logic gates
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
