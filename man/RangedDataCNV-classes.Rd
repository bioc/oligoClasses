\name{RangedData-classes}
\Rdversion{1.1}
\docType{class}
\alias{RangedDataCopyNumber-class}
\alias{RangedDataCNV-class}
\alias{RangedDataCBS-class}
\alias{RangedDataHMM-class}
\alias{RangedDataCopyNumber}
\alias{coverage2,RangedDataCNV-method}
\alias{coverage2}
\alias{mean,RangedDataCBS-method}
\alias{featuresInRange,SnpSet,RangedDataCNV-method}
\alias{sampleNames,RangedDataCNV-method}
\alias{state,RangedDataCNV-method}
\alias{state}
\alias{todf,RangedDataCNV,ANY-method}
\alias{coerce,RangedData,RangedDataCBS-method}
\alias{coerce,RangedData,RangedDataHMM-method}
\alias{coerce,RangedDataHMM,GRangesList-method}
\alias{coerce,RangedDataCNV,GRangesList-method}
\alias{sampleNames<-,RangedDataCNV,character-method}
%\alias{trioNames,RangedDataCNV-method}

\title{Classes in MinimumDistance for data on ranges}
\description{

  \code{RangedDataCNV} is a class extending the virtual class
  \code{RangedDataCopyNumber}. \code{RangedDataCopyNumber} extends
  \code{\linkS4class{RangedData}}.

  \code{RangedDataCBS} extends \code{\linkS4class{RangedDataCNV}} and is
  useful for storing ranges from segmentation algorithms such as
  circular binary segmentation.  In particular, the columns 'chrom',
  'id', and 'num.mark' are required for instances of the class.

  \code{RangedDataHMM} extends \code{\linkS4class{RangedDataHMM}} and is
  useful for storing ranges from hidden Markov models such as PennCNV or
  VanillaICE.  A column labeled 'state' is required for the class.

}
\section{Objects from the Class}{

See \code{RangedDataCBS} and \code{RangedDataHMM} constructors.

}
\section{Slots}{
  \describe{
    \item{\code{ranges}:}{Object of class \code{"RangesList"} ~~ }
    \item{\code{values}:}{Object of class \code{"SplitDataFrameList"} ~~ }
    \item{\code{elementType}:}{Object of class \code{"character"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataTableORNULL"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{RangedDataCNV}"}, directly.
Class \code{"\linkS4class{RangedDataCopyNumber}"}, by class
"RangedDataCNV", distance 2.
%Class \code{RangedDataCopyNumber}, by class "RangedDataCNV", distance 2.
Class \code{"\linkS4class{RangedData}"}, by class "RangedDataCNV", distance 3.
Class \code{"\linkS4class{DataTable}"}, by class "RangedDataCNV", distance 4.
Class \code{"\linkS4class{List}"}, by class "RangedDataCNV", distance 4.
Class \code{"\linkS4class{DataTableORNULL}"}, by class "RangedDataCNV", distance 5.
Class \code{"\linkS4class{Vector}"}, by class "RangedDataCNV", distance 5.
Class \code{"\linkS4class{Annotated}"}, by class "RangedDataCNV", distance 6.
}

\section{Accessors for \code{RangedDataCNV}}{
  In the following code descriptions, \code{object} is a class derived
  from \code{RangedData}.
  \describe{

    \item{}{

      \code{as(object, "RangedDataCBS")}, \code{as(object, "RangedDataHMM")}:

      coerce an object of class \code{RangedData} to
      \code{RangedDataCBS} or \code{RangedDataHMM}, respectively.

    }

    \item{}{

      \code{as(object, "GRangesList")}:

      coerce an object of class \code{RangedDataHMM} or
      \code{RangedDataCNV} to a \code{GRangesList}.

    }

    \item{}{\code{coverage2(object)}: numeric vector indicating the
      number of markers in each range.  The length of the vector is the
      same as \code{nrow(object)}.}

    \item{}{\code{mean(object)}: returns the mean copy number for each
    segment in \code{object}.}

    \item{}{\code{sampleNames(object)}: character vector of sample ids
    for each range.  The length of the vector is the same as
    \code{nrow(object)}.}

    \item{}{\code{sampleNames(object) <- value}: Assign the
    \code{character} fector \code{value} to the sample identifiers in a
    \code{RangedDataCNV}-derived object.}

    \item{}{\code{state(object)}: For objects of class
  \code{RangedDataHMM} derived by fitting a hidden Markov model, the
  function \code{state} returns the index of the hidden state.}

}
}

\section{Accessors for \code{RangedDataCBS} objects}{
  \describe{
    In the following code descriptions, \code{object} belongs to class
    \code{RangedDataCBS}
    \item{}{\code{mean(object)}: numeric vector of segment means}
  }
}

\section{Accessors for \code{RangedDataHMM} objects}{
  \describe{
    In the following code descriptions, \code{object} belongs to class
    \code{RangedDataCBS}

    \item{}{\code{state(object)}: character vector of hidden states from
      the HMM}

  }
}

\section{Coercion to data.frame}{
  \describe{

    \item{}{\code{todf(object)}: coersion of a \code{RangedDataCNV} to a
    \code{data.frame}.  Not intended to be called directly by the user. }

  }
}

\author{
R. Scharpf
}

\seealso{

  See \code{\linkS4class{RangedData}} for a detailed description and
  additional methods.  See the constructors \code{\link{RangedDataCNV}},
  \code{\link{RangedDataCBS}}, and \code{\link{RangedDataHMM}} for
  generating an object of the respective class.

}
\examples{showClass("RangedDataCNV")}
\keyword{classes}
