.onAttach <- function(lib, pkg)
{
  #cat("\n \n Welcome to nontarget version 1.2 \n \n");
  library.dynam("nontarget", pkg, lib);
}

packageStartupMessage
{
  cat("\n \n Welcome to nontarget version 1.2 \n \n");
  #library.dynam("nontarget", pkg, lib);
}
