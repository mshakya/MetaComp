% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_gottcha_assignment.R
\name{load_gottcha_assignment}
\alias{load_gottcha_assignment}
\title{Efficiently loads a GOTTCHA (or other EDGE-like taxonomic assignment) table from a file.
An assumption has been made -- since GOTTCHA/EDGE tables are generated in an automated fashion,
they should be properly formatted -- thus the code doesn't check for any inconsistencies except
for the very file existence. This implementation fully relies on the read.table function
from data.table package gaining performance over traditional R techniques.}
\usage{
load_gottcha_assignment(filepath)
}
\arguments{
\item{filepath}{A path to EDGE-generated tab-delimeted GOTTCHA taxonomy assignment file.}
}
\value{
a data frame containing four columns: TAXA, LEVEL, COUNT, and ABUNDANCE, representing
        taxonomically anchored sequences from the sample.
}
\description{
Efficiently loads a GOTTCHA (or other EDGE-like taxonomic assignment) table from a file.
An assumption has been made -- since GOTTCHA/EDGE tables are generated in an automated fashion,
they should be properly formatted -- thus the code doesn't check for any inconsistencies except
for the very file existence. This implementation fully relies on the read.table function
from data.table package gaining performance over traditional R techniques.
}

