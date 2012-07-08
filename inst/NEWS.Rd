\name{NEWS}
\title{News for Package 'oligoClasses'}
\section{Changes in version 1.18}{
  \subsection{New class: GenomeAnnotatedDataFrame}{
  \itemize{
    \item GenomeAnnotatedDataFrame extends AnnotatedDataFrame but requires
    columns 'isSnp', 'chromosome', and 'position'.

    \item Accessors for the SNP indicator, chromosome, and position are
    defined for the GenomeAnnotatedDataFrame class

    \item The assayData elements must all be integers. Copy number or
      relative copy number should be multipled by 100 and stored as an
      integer. Use integerMatrix(x, 100), where x is a matrix of copy number
      on the original scale.  B allele frequencies should be multipled by
      1000 and saved as an integer. Use integerMatrix(b, 1000), where 'b' is
      a matrix of B allele frequencies.

      \item The featureData in these classes should have the class
      "GenomeAnnotatedDataFrame".

    }
  }
}


\section{Changes in version 1.14}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item ocSamples() and ocProbesets() propagate the value they are
      set to to a cluster (if the cluster is set);
    }
  }
}

\section{Changes in version 1.9}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Major refactoring of classes used by SNP/CNP algorithms.
      \enumerate{
	\item New class: CNSet, AlleleSet, SnpSuperSet
	\item Removed classes: Snp*QSet, Snp*CallSet, Snp*CallSetPlus
	\item Genotyping algorithms will use SnpSet instead of SnpCallSet
      }
      \item Replaced thetaA, thetaB, senseThetaA, antisenseThetaA,
      senseThetaB and antisenseThetaB methods by allele.
      \item Exported DBPDInfo, ExpressionPDInfo, FeatureSet, TilingPDInfo
      \item Added ff utilities used by both oligo and crlmm packages
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item bothStrands() methods for AlleleSet and SnpFeatureSet.
    }
  }
}
