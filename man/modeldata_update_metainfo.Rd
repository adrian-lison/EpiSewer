% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sewer_modeldata.R
\name{modeldata_update_metainfo}
\alias{modeldata_update_metainfo}
\title{Update meta information in modeldata based on available variables}
\usage{
modeldata_update_metainfo(modeldata)
}
\description{
This update function is for all meta information that depends on
modeldata from several functions. There are also functions which add meta
information directly, in particular when this meta information cannot be
solely inferred from modeldata. This function is designed such that calling
it will never do harm to the modeldata object and not throw errors if
something is missing in the modeldata object.
}
\keyword{internal}
