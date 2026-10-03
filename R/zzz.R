.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "PhysioNetPhysiology ", utils::packageVersion(pkgname),
    " - Time Delay Stability & organ-interaction networks")
}
