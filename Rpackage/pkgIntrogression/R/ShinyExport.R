IntrogressionShiny <- function() {
  appDir <- system.file("ShinyApp", package = "pkgIntrogression")
  if(appDir == "") {
    stop("Could not find directory. Try re-installing `pkgIntrogression`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
