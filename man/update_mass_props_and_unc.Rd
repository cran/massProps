% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/massprops.R
\name{update_mass_props_and_unc}
\alias{update_mass_props_and_unc}
\title{Update mass properties and uncertainties}
\usage{
update_mass_props_and_unc(
  df,
  target,
  sources,
  override = set_poi_conv_from_target
)
}
\arguments{
\item{df}{A data frame  with (at least) these columns: \code{id}, \code{mass}, \code{Cx},
\code{Cy}, \code{Cz}, \code{Ixx}, \code{Iyy}, \code{Izz}, \code{Ixy}, \code{Ixz}, \code{Iyz}, \code{POIconv}, \code{Ipoint},
\code{sigma_mass}, \code{sigma_Cx}, \code{sigma_Cy}, \code{sigma_Cz}, \code{sigma_Ixx}, \code{sigma_Iyy},
\code{sigma_Izz}, \code{sigma_Ixy}, \code{sigma_Ixz}, \code{sigma_Iyz}.}

\item{target}{The \code{id} value of the target row.}

\item{sources}{List of \code{id} values of the of the source rows.}

\item{override}{An override function, called as override(df, target, value). The default override sets the POI sign convention
of a computed aggregate to the \code{POIconv} column of the target row in the data frame.}
}
\value{
The updated data frame.
}
\description{
\code{update_mass_props_and_unc()} updates mass properties and uncertainties
for a specified target row from
specified source rows in a data frame.
}
\examples{
leaves <- list("Widget", "2nd Part")
update_mass_props_and_unc(sawe_table, "Combined", leaves)

}
