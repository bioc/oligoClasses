\name{NEWS}
\title{News for Package 'oligoClasses'}

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
