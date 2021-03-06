% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_database.R
\name{read_database}
\alias{read_database}
\title{Read Database}
\usage{
read_database(path)
}
\arguments{
\item{path}{path to saved sqlite database}
}
\value{
An RSQLite database to memory
}
\description{
This function reads an sqlite database from disk to memory
}
\examples{
\dontrun{
DBs <- read_database(path = "database.sqlite")
}
}
