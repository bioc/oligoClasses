\name{checkExists}
\alias{checkExists}
\title{
	Checks to see whether an object exists and, if not, executes the appropriate function.
}
\description{

	Only loads an object if the object name is not in the global
	environment.  If not in the global environment and the file
	exists, the object is loaded (by default). If the file does
	not exist, the function FUN is run.

}

\usage{
checkExists(.name, .path = ".", .FUN, .FUN2, .save.it=TRUE, .load.it, ...)
}


\arguments{
  \item{.name}{
  Character string giving name of object in global environment
}
  \item{.path}{
  Path to where the object is saved.
}
  \item{.FUN}{
  Function to be executed if <name> is not in the global environment and the file does not exist.
}
  \item{.FUN2}{

  Not currently used.
	  
}
 \item{.save.it}{

 Logical.  Whether to save the object to the directory indicaged by
 \code{path}. This argument is ignored if the object was loaded from
 file or already exists in the .GlobalEnv.

}

\item{.load.it}{

Logical.  If load.it is TRUE, we try to load the object from the
indicated \code{path}.  The returned object will replace the object in
the .GlobalEnv unless the object is bound to a different name (symbol)
when the function is executed.

}

\item{\dots}{

	Additional arguments passed to FUN.
  
}
}

\value{

	Could be anything -- depends on what FUN, FUN2 perform.  

	Future versions could return a 0 or 1 indicating whether the
	function performed as expected.

}

\author{
R. Scharpf
}

\examples{
path <- tempdir()
dir.create(path)
x <- 3+6
x <- checkExists("x", .path=path, .FUN=function(y, z) y+z, y=3, z=6)
rm(x)
x <- checkExists("x", .path=path, .FUN=function(y, z) y+z, y=3, z=6)
rm(x)
x <- checkExists("x", .path=path, .FUN=function(y, z) y+z, y=3, z=6)
rm(x)
##now there is a file called x.rda in tempdir().  The file will be loaded
x <- checkExists("x", .path=path, .FUN=function(y, z) y+z, y=3, z=6)
rm(x)
unlink(path, recursive=TRUE)
}
\keyword{manip}

