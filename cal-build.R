args <- commandArgs(trailingOnly = TRUE)
action <- if (length(args) == 0L) "cal" else args[[1L]]

if (identical(action, "clean")) {
  bookdown::clean_book(clean = TRUE)
} else {
  if (file.exists("_main.Rmd")) unlink("_main.Rmd")
  bookdown::render_book(
    paste0(action, ".Rmd"),
    "bookdown::pdf_book",
    config_file = paste0(action, ".yml")
  )
}
