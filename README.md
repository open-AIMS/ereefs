<!-- README.md is generated from README.Rmd. Please edit that file -->

# ereefs <img src="man/figures/logo.png" width = 180 alt="eReefs Logo" align="right" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- [![R build status](https://github.com/open-AIMS/ereefs/actions/workflows/R-CMD-check.yaml/badge.svg?branch=dev)](https://github.com/open-AIMS/ereefs/actions) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/open-aims/ereefs/branch/master/graph/badge.svg)](https://app.codecov.io/gh/open-aims/ereefs?branch=master) -->
<!-- ![pkgdown](https://github.com/open-AIMS/ereefs/workflows/pkgdown/badge.svg) -->
[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
[![Ask Us Anything
!](https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg)](https://github.com/open-AIMS/ereefs/issues/new)
![Open Source
Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
<!-- badges: end -->

## Overview

**What is eReefs?**

[eReefs](https://www.ereefs.org.au/about/) is a project that combines
government commitment to reef protection, world-class science innovation
and contributions from leading Australian businesses. It is a
collaborative information system created (the [Great Barrier Reef
Foundation](https://www.barrierreef.org/),
[CSIRO](https://www.csiro.au/), the [Australian Institute of Marine
Science](https://www.aims.gov.au/), [Bureau of
Meteorology](https://www.bom.gov.au/), and [Queensland
Government](https://www.qld.gov.au/)) that provides a picture of what is
currently happening on the reef, and what will likely happen in the
future.

Focused on the protection and preservation of the iconic Great Barrier
Reef, it forms the first step in building comprehensive coastal
information systems for Australia. Using the latest technologies to
collate data, and new and integrated modelling, eReefs will produce
powerful visualisation, communication and reporting tools. It will
provide for the Reef information akin to that provided by the Bureau of
Meteorology for weather. This information will benefit government
agencies, Reef managers, policy makers, researchers, industry and local
communities.

**What does the `ereefs` R package do?**

The `ereefs` R package provides easy access to
[eReefs](https://www.ereefs.org.au/about/) and other
[CSIRO-EMS](https://research.csiro.au/cem/software/ems/) output files.
It is designed to assist R users who need more customised access to
eReefs data. This includes things like:

-   Accessing data from versions of the eReefs model that are not
    available through the web-based data service;

-   Accessing data for less commonly-used variables that are not
    included in the web-based data service or visualisation portal;

-   Generating customised maps or animations for a specified region at a
    specified depth over a specified period of time, including (for
    example) true colour maps and animations from modelled optical data,
    maps that combine two or more eReefs variables (e.g., to get total
    zooplankton concentrations), maps with labelled points of interest
    such as Marine Monitoring Program
    ([MMP](https://www2.gbrmpa.gov.au/our-work/programs-and-projects/marine-monitoring-program))
    sampling locations, and maps with outlines of reef areas shown. Maps
    are created as [`ggplot2`](https://ggplot2.tidyverse.org/) figures
    and can be further adjusted as required by R users familiar with
    that package.

-   Extracting and visualising data in different ways, for example:

    -   Taking a vertical profile of a variable over the depth of the
        water column and displaying it either as a single profile at a
        point in time or as a depth-vs-time contour plot;  
    -   Taking a two-dimensional vertical slice through the
        three-dimensional model data;  
    -   Calculating vertically integrated results (i.e., the average
        value of a variable over the depth of the water column rather
        than the value at a particular depth);  
    -   Extracting data along the path of a boat or glider.

**How does it compare to the existing web portal?**

If you want to extract time-series of the most commonly needed physical
and chemical variables from the main versions of the eReefs model
outputs at a fixed set of geolocations at a fixed depth from mean sea
level, the [web-based eReefs data extraction
tool](https://extraction.ereefs.aims.gov.au/) is usually faster and
easier to use. If you want animations of the most commonly-needed
variables at the surface, and do not need them to be customised, the
[AIMS eReefs visualisation
portal](https://ereefs.aims.gov.au/ereefs-aims) offers pre-generated
animations at several scales in an easy-to-navigate format.

## Installation

To install the latest release from GitHub use

    if (!requireNamespace("remotes")) {
      install.packages("remotes")
    }
    remotes::install_github("open-aims/ereefs")

The current development version can be downloaded from GitHub via

    if (!requireNamespace("remotes")) {
      install.packages("remotes")
    }
    remotes::install_github("open-aims/ereefs", ref = "dev")

## Usage

Usage and further information about `ereefs` can be seen on the [project
page](https://open-aims.github.io/ereefs/) and the
[vignettes](https://open-aims.github.io/ereefs/articles/). Help files
for the individual functions can be found on the [reference
page](https://open-aims.github.io/ereefs/reference/).

## Further Information

`ereefs` is provided by the [Australian Institute of Marine
Science](https://www.aims.gov.au) under the MIT License
([MIT](https://opensource.org/licenses/MIT)).

## Other related content

-   [Access](https://ereefs.aims.gov.au/ereefs-aims) to visualisations
    that will give you an idea of the types of information that are in
    the eReefs models.

-   [Snippets](https://ereefs.aims.gov.au/ereefs-aims/help) of code for
    interactive with the eReefs data, including in
    [python](https://ereefs.aims.gov.au/ereefs-aims/help/how-to-plot-aims-ereefs-data-with-python).

-   [CSIRO public data
    service](https://dapds00.nci.org.au/thredds/catalogs/fx3/catalog.html)
    providing access to the raw model data.

-   [Access](https://thredds.ereefs.aims.gov.au/thredds/s3catalogue/aims-ereefs-public-prod/derived/ncaggregate/ereefs/catalog.html)
    to the derived aggregated AIMS data used in the aggregate products
    on the AIMS visualisation portal. In general AIMS only keeps the
    important variables, half the depths, and aggregations (raw hourly
    hydro -&gt; daily, monthly, annual; raw daily BGC -&gt; monthly,
    annual). AIMS derived products, along with temporal aggregation are
    also resampled onto a regular grid to allow these files to be
    processed in a GIS application. The raw model data is on a
    curvilinear grid which requires special consideration when
    manipulating in code.

-   The [GBRF
    website](https://www.barrierreef.org/what-we-do/projects/eReefs),
    with videos and descriptions by comms professionals.
