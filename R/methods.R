#' Show
#'
#' @importFrom methods getMethod
#'
#' @keywords internal
#' @noRd
setMethod(
    f = "show",
    signature = "SpatialExperiment",
    definition = function(object) {
        if (
            isTRUE(x = getOption(x = "restore_SpatialExperiment_show", default = FALSE))
        ) {
            f <- methods::getMethod(
                f = "show",
                signature = "SummarizedExperiment",
                where = asNamespace(ns = "SummarizedExperiment")
            )
            f(object = object)

        } else {  
            print(object)
        }
    }
)

setClass("tidySpatialExperiment", contains = "SpatialExperiment")

#' Extract and join information for features.
#'
#' @description join_features() extracts and joins information for specified features
#'
#' @importFrom ttservice join_features
#'
#' @name join_features
#' @rdname join_features
#'
#' @param .data A SpatialExperiment object
#' @param features A vector of feature identifiers to join
#' @param all If TRUE return all
#' @param exclude_zeros If TRUE exclude zero values
#' @param shape Format of the returned table "long" or "wide"
#' @param ... Parameters to pass to join wide, i.e. assay name to extract feature abundance from and gene prefix, for shape="wide"
#'
#' @details This function extracts information for specified features and returns the information in either long or wide format.
#'
#' @return An object containing the information.for the specified features
#'
#' @examples
#' example(read10xVisium)
#' spe |>
#'     join_features(features = "ENSMUSG00000025900")
NULL

#' join_features
#'
#' @keywords internal
#' @noRd
setMethod("join_features", "SpatialExperiment", function(.data, features = NULL, all = FALSE,
                                                        exclude_zeros = FALSE, shape = "long", 
                                                        ...) {
      
      # Define unbound variable
      index <- NULL

      # Shape is long
      if (shape == "long") {
          
          # Print message about changing data type
          message(
              "tidySpatialExperiment says: A data frame is returned for independent data 
              analysis."
          )
        
          # Join feature abundance with colData by index
          .data |>
              colData() |>
              tibble::as_tibble(rownames = c_(.data)$name) |>
              tibble::rowid_to_column("index") |>
              dplyr::mutate(index = as.character(index)) |>
              dplyr::left_join(
                  get_abundance_sc_long(
                      .data = .data,
                      features = features,
                      all = all,
                      exclude_zeros = exclude_zeros
                  ),
                  by = "index"
              ) |>
              dplyr::mutate("index" = NULL)
      
      # Shape is wide
      } else {
        colData(.data) <- 
          .data |>
              colData() |>
              tibble::as_tibble(rownames = c_(.data)$name) |>
              tibble::rowid_to_column("index") |>
              dplyr::mutate(index = as.character(index)) |>
              left_join(
                  get_abundance_sc_wide(
                      .data = .data,
                      features = features,
                      all = all, ...
                  ),
                  by = "index") |>
              dplyr::mutate("index" = NULL) |>
              as_meta_data(.data)
          .data
      } 
})

#' Aggregate cells
#'
#' @description Combine cells into groups based on shared variables and aggregate feature counts.
#'
#' @importFrom rlang enquo
#' @importFrom tibble enframe
#' @importFrom Matrix rowSums
#' @importFrom dplyr full_join
#'
#' @name aggregate_cells
#' @rdname aggregate_cells
#' 
#' @param .data A tidySpatialExperiment object
#' @param .sample A vector of variables by which cells are aggregated
#' @param slot The slot to which the function is applied
#' @param assays The assay to which the function is applied
#' @param aggregation_function The method of cell-feature value aggregation
#' 
#' @return A SummarizedExperiment object
#' 
#' @examples 
#' example(read10xVisium)
#' spe |>
#'     aggregate_cells(sample_id, assays = "counts")
#'
#' @export
aggregate_cells <- function(.data, .sample = NULL, slot = "data", assays = NULL, 
                            aggregation_function = rowSums) {
  
    # Declare unbound variables
    feature <- NULL
    .sample <- enquo(.sample)
    
    # Subset only wanted assays
    if (!is.null(assays)) {
        .data@assays@data <- .data@assays@data[assays]
    }
    
    .data |>
      
        nest(data = -!!.sample) |>
        mutate(.aggregated_cells = as.integer(map(data, ~ ncol(.x)))) |>
        mutate(data = map(data, ~ 
            
            # loop over assays
            map2(
                as.list(assays(.x)), names(.x@assays),
                
                # Get counts
                ~  .x |>
                  aggregation_function(na.rm = TRUE) |>
                  enframe(
                      name  = "feature",
                      value = sprintf("%s", .y)
                  ) |>
                  mutate(feature = as.character(feature)) 
            ) |>
            reduce(function(...) dplyr::full_join(..., by = c("feature")))

        )) |>
        left_join(
            .data |> 
            as_tibble() |> 
            subset(!!.sample),
        by = quo_names(.sample)
        ) |>
        unnest(data) |>
    
    drop_class("tidySpatialExperiment_nested") |> 
    
    as_SummarizedExperiment(
        .sample = !!.sample, .transcript = feature, .abundance = !!as.symbol(names(.data@assays))
    )
}

#' Rectangle Gating Function
#'
#' @description Determines whether points specified by spatial coordinates are within a defined rectangle.
#'
#' @importFrom dplyr mutate
#'
#' @name rectangle
#' @rdname rectangle
#' 
#' @param spatial_coord1 Numeric vector for x-coordinates (e.g., array_col)
#' @param spatial_coord2 Numeric vector for y-coordinates (e.g., array_row)
#' @param center Numeric vector of length 2 specifying the center of the rectangle (x, y)
#' @param height The height of the rectangle
#' @param width The width of the rectangle
#' 
#' @return Logical vector indicating points within the rectangle
#' 
#' @examples 
#' example(read10xVisium)
#' spe |>
#'     mutate(in_rectangle = rectangle(
#'       array_col, array_row, center = c(50, 50), height = 20, width = 10)
#'       )
#'
#' @export
rectangle <- function(spatial_coord1, spatial_coord2, center, height, width) {
    x_min <- center[1] - width / 2
    x_max <- center[1] + width / 2
    y_min <- center[2] - height / 2
    y_max <- center[2] + height / 2

    within_x <- spatial_coord1 >= x_min & spatial_coord1 <= x_max
    within_y <- spatial_coord2 >= y_min & spatial_coord2 <= y_max

    return(within_x & within_y)
}

#' Ellipse Gating Function
#'
#' @name ellipse
#' @rdname ellipse
#' @description Function to create an ellipse gate in a SpatialExperiment object
#' @param spatial_coord1 Numeric vector for x-coordinates
#' @param spatial_coord2 Numeric vector for y-coordinates
#' @param center Numeric vector (length 2) for ellipse center (x, y)
#' @param axes_lengths Numeric vector (length 2) for the lengths of the major and minor axes of the
#' ellipse
#' @return Logical vector indicating points within the ellipse
#' @examples
#' example(read10xVisium)
#' spe |>
#'     mutate(in_ellipse = ellipse(
#'         array_col, array_row, center = c(50, 50), axes_lengths = c(20, 10))
#'     )
#'     
#' @export
ellipse <- function(spatial_coord1, spatial_coord2, center, axes_lengths) {
    # axes_lengths should be a vector of length 2: [major_axis, minor_axis]

    # Scaling factor to normalize the ellipse to a unit circle
    scale_x <- 1 / axes_lengths[1]
    scale_y <- 1 / axes_lengths[2]

    # Normalized coordinates relative to ellipse center
    normalized_x <- (spatial_coord1 - center[1]) * scale_x
    normalized_y <- (spatial_coord2 - center[2]) * scale_y

    # Check if points are within the unit circle (ellipse after normalization)
    within_ellipse <- (normalized_x^2 + normalized_y^2) <= 1

    return(within_ellipse)
}

#' Gate spatial
#' 
#' Interactively gate points by their location in space, with image data overlaid. 
#' 
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom magick image_read
#' @importFrom magick image_info
#' @importFrom plotly raster2uri
#' @importFrom plotly ggplotly
#' @importFrom plotly layout
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 coord_cartesian
#' @importFrom SpatialExperiment imgData
#' @importFrom SpatialExperiment imgRaster
#' @importFrom SpatialExperiment imgSource
#' @importFrom shiny shinyApp
#' @importFrom shiny runApp
#' @importFrom tidygate ui
#' @importFrom tidygate server
#' @examples
#' \dontrun{
#' example(read10xVisium)
#' spe |>
#'     gate_spatial()
#' }
#' @param spe A SpatialExperiment object
#' @param image_index The image to display if multiple are stored within the provided
#' SpatialExperiment object. 
#' @param alpha The opacity of points, with 1 being completely opaque and 0 being completely
#' transparent.
#' @param size The size of points.
#' @return The input SpatialExperiment object with a new column `.selected`, with TRUE for elements 
#' inside gate points and FALSE for elements outside gate points. A record of the selected points is 
#' stored in `tidygate_env$select_data` and a record of the gates is stored in 
#' `tidygate_env$brush_data`.
#' @export
gate_spatial <-
  
  function(spe, image_index = 1, alpha = 1, size = 1) {
    
    message("tidySpatialExperiment says: this feature is in early development and may undergo changes or contain bugs.")
    
    # Create tibble with necessary information
    data <- 
      tibble::tibble(
        dimension_x = 
          spe |>
          pull(pxl_col_in_fullres), 
        dimension_y = 
          spe |>
          pull(pxl_row_in_fullres)
      ) |>
      dplyr::mutate(.key = dplyr::row_number()) |>
      dplyr::mutate(.selected = FALSE)
    
    # Prepare spatial information
    image <- 
      SpatialExperiment::imgData(spe)[1, ]@listData$data[[1]] |>
      SpatialExperiment::imgSource() |>
      magick::image_read()
    
    image_x_size <- 
      magick::image_info(image)$width /
      SpatialExperiment::imgData(spe)[image_index, ]@listData$scaleFactor
    
    image_y_size <- 
      magick::image_info(image)$height /
      SpatialExperiment::imgData(spe)[image_index, ]@listData$scaleFactor
    
    image_uri <- 
      SpatialExperiment::imgData(spe)[image_index, ]@listData$data[[1]] |>
      SpatialExperiment::imgRaster() |>
      plotly::raster2uri()
    
    # Create plot 
    spatial_plot <-
      data |>
      ggplot2::ggplot(ggplot2::aes(x = dimension_x, y = dimension_y, key = .key)) +
      ggplot2::geom_point(alpha = alpha, size = size) +
      ggplot2::coord_fixed(
        xlim =  c(0, image_x_size), 
        ylim = rev(c(0, image_y_size)), 
        expand = FALSE,
        ratio = 1
      )

    spatial_plot <- 
      spatial_plot |>
      plotly::ggplotly() |>
      plotly::layout(images = list(
        list(
          source = image_uri,
          xref = "x",
          yref = "y",
          x = 0,
          y = 0,
          sizex = image_x_size,
          sizey = image_y_size,
          sizing = "stretch",
          opacity = 1,
          layer = "below"
        )
      ))
    
    # Create environment and save input variables
    tidygate_env <<- rlang::env()
    tidygate_env$input_data <- data
    tidygate_env$dimension_x <- data$dimension_x
    tidygate_env$dimension_y <- data$dimension_y
    tidygate_env$custom_plot <- spatial_plot
    
    # Launch Shiny App
    app <- shiny::shinyApp(tidygate::ui, tidygate::server)
    shiny::runApp(app, port = 1234)
    
    spe$.gate <- tidygate_env$select_data$.selected 
    
    return(spe)
  }
