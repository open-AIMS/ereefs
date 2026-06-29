#' ereefs package
#'
#' @description A package for accessing, extracting, and visualising eReefs and
#'   other CSIRO-EMS model output. The current development workflow is designed
#'   around local NetCDF files, live OPeNDAP datasets, and THREDDS catalog
#'   services, with support for both curvilinear EMS grids and centre-only
#'   regular grids.
#'
#' @name ereefs-package
#' @aliases ereefs
#'
#' @references
#' https://www.ereefs.org.au/about/
#' https://research.csiro.au/ereefs/welcome-to-ereefs/about-ereefs/
#' https://ereefs.aims.gov.au/ereefs-aims
#' https://www.barrierreef.org/what-we-do/projects/eReefs
#'
#' @importFrom chron chron
#' @importFrom cowplot plot_grid
#' @importFrom dplyr arrange between bind_rows coalesce distinct
#'   filter group_by inner_join left_join mutate pull relocate rename right_join
#'   row_number rowwise select slice summarise transmute ungroup
#' @importFrom grDevices col2rgb rgb
#' @importFrom jsonlite fromJSON
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom utils tail
#' @keywords internal
"_PACKAGE"
#
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:18
# - date: 2026-04-28
# - prompt_used: "Do a final documentation polish pass so the help text reads consistently with the current tidy-dev workflow."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:46
# - date: 2026-06-29
# - prompt_used: "Save the current instruction block locally, then finish GitHub issues #24 and #25 by tightening argument documentation and making package imports/dependencies explicit."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:52
# - date: 2026-06-29
# - prompt_used: "Correct the package-level import sources, regenerate NAMESPACE/man pages, and verify the documentation matches the active refactored functions."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:24
# - date: 2026-06-29
# - prompt_used: "Resolve the remaining namespace/import warnings reported by R CMD check for issues #24 and #25."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:29
# - date: 2026-06-29
# - prompt_used: "Add the remaining explicit package imports used by plotting helpers so the package check is cleaner for issues #24 and #25."
