`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
repo_lib <- normalizePath(
  file.path(dirname(dirname(script_path)), ".r_libs"),
  winslash = "/",
  mustWork = TRUE
)

.libPaths(c(repo_lib, .libPaths()))

cat("libPaths=", paste(.libPaths(), collapse = " | "), "\n")

load_error <- function(pkg) {
  tryCatch(
    {
      loadNamespace(pkg)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

bslib_error <- load_error("bslib")
shiny_error <- load_error("shiny")
sass_error <- load_error("sass")
httpuv_error <- load_error("httpuv")

cat("bslib=", bslib_error %||% "ok", "\n")
cat("sass=", sass_error %||% "ok", "\n")
cat("shiny=", shiny_error %||% "ok", "\n")
cat("httpuv=", httpuv_error %||% "ok", "\n")

if (!is.null(bslib_error) || !is.null(shiny_error)) {
  stop(
    paste(
      "Windows GUI runtime stack is still incomplete in this environment.",
      paste0("bslib: ", bslib_error %||% "ok"),
      paste0("shiny: ", shiny_error %||% "ok"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

theme <- bslib::bs_theme(version = 5, bootswatch = "minty")
page <- bslib::page_fluid(
  theme = theme,
  bslib::card(
    bslib::card_header("Header"),
    "Body"
  ),
  htmltools::tags$p("hello")
)

rendered <- htmltools::renderTags(page)
cat("html_chars=", nchar(rendered$html), "\n")
cat("dependencies=", length(rendered$dependencies), "\n")

output_file <- tempfile(fileext = ".html")
htmltools::save_html(page, output_file)
cat("saved_html=", file.exists(output_file), "\n")
cat("output_file=", output_file, "\n")

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 10:56
# - date: 2026-07-06
# - prompt_used: "Set up a quick test to see whether the earlier Windows bslib GUI problem is still present in this environment."
