#' getComArgs function
#'
#' function to (do something)
#'
#' @param defaultArgs [value]. Default is NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getComArgs <- function(defaultArgs = NULL) {
  defaults <- !is.null(defaultArgs)
  args <- commandArgs(TRUE)
  isAssn <- grepl("=", args)
  userArgs <- args[isAssn]
  needEval <- grepl("\\(|\\)|\\:|\\'", userArgs) 
  argSplit <- strsplit(userArgs, "=")
  argList <- lapply(argSplit, "[[", 2)
  names(argList) <- lapply(argSplit, "[[", 1)
  argList[needEval] <- lapply(argList[needEval], function(x) eval(parse(text = x)))
  argList[!needEval] <- lapply(argList[!needEval], function(x) strsplit(x, ",")[[1]])
  argList[!needEval] <- type.convert(argList[!needEval], as.is = TRUE)
  print(argList)
  if (defaults){
    defaultArgs[names(argList)] <- argList
    return(defaultArgs)
  } else {
    return(argList)
  }
}
