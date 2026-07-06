test_that("Windows GUI runtime can load and render a minimal bslib page", {
  repo_lib <- normalizePath(
    testthat::test_path("..", "..", ".r_libs"),
    winslash = "/",
    mustWork = TRUE
  )
  old_libpaths <- .libPaths()
  on.exit(.libPaths(old_libpaths), add = TRUE)
  .libPaths(c(repo_lib, old_libpaths))

  load_error <- function(pkg) {
    tryCatch(
      {
        loadNamespace(pkg)
        NULL
      },
      error = function(e) conditionMessage(e)
    )
  }

  or_else <- function(x, fallback) {
    if (is.null(x)) {
      fallback
    } else {
      x
    }
  }

  bslib_error <- load_error("bslib")
  shiny_error <- load_error("shiny")

  if (!is.null(bslib_error) || !is.null(shiny_error)) {
    testthat::fail(
      paste(
        "Windows GUI runtime stack is still incomplete in this environment.",
        paste0("bslib: ", or_else(bslib_error, "ok")),
        paste0("shiny: ", or_else(shiny_error, "ok")),
        sep = "\n"
      )
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
  expect_gt(nchar(rendered$html), 0)
  expect_gt(length(rendered$dependencies), 0)

  output_file <- tempfile(fileext = ".html")
  htmltools::save_html(page, output_file)
  expect_true(file.exists(output_file))
})

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 10:56
# - date: 2026-07-06
# - prompt_used: "Set up a quick test to see whether the earlier Windows bslib GUI problem is still present in this environment."
