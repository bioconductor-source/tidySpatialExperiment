#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr mutate
#' @importFrom utils tail
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySpatialExperiment
#' @param features A character
#' @param all A boolean
#' @param ... Parameters to pass to join wide, i.e. assay name to extract feature abundance from
#'
#' @importFrom SummarizedExperiment assays
#'
#' @return A tidySpatialExperiment object
#'
#' @noRd
get_abundance_sc_wide <- function(.data, features = NULL, all = FALSE, 
                                  assay = SummarizedExperiment::assays(.data) |> as.list() |> 
                                  tail(1) |> names(),  prefix = "") {

    # Define unbound variable
    index <- NULL
  
    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(features) &
            all == FALSE
    ) {
        stop(
            "Your object does not contain variable feature labels,
            feature argument is empty and all arguments are set to FALSE.
            Either:
            1. use detect_variable_features() to select variable feature
            2. pass an array of feature names
            3. set all=TRUE (this will output a very large object, does your computer have enough 
            RAM?)"
        )
    }

    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(features) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } else {
        variable_genes <- NULL
    }

    # Just grub last assay
    assays(.data) |>
        as.list() |>
        pluck(assay) |>
        when(
            variable_genes |>
              is.null() |>
              lapply(`!`) |>
              unlist() ~ 
                (.)[variable_genes, , drop = FALSE],
            features |>
                is.null() |>
                lapply(`!`) |>
                unlist() ~ 
                (.)[features, , drop = FALSE],
            ~ stop(
                "It is not convenient to extract all genes, you should have either variable features
                or feature list to extract"
            )
        ) |>
        as.matrix() |>
        t() |>
        as_tibble() |>
      
        # Add index for joining
        tibble::rowid_to_column("index") |>
        dplyr::mutate(index = as.character(index))
}

#' Bind columns without checking for duplicate sample_ids
#' 
#' @importFrom SingleCellExperiment int_metadata 
#' @importFrom SpatialExperiment imgData
#' @importFrom BiocGenerics rbind cbind
#' @importFrom methods callNextMethod
#' 
#' @keywords internal
#' @noRd
setMethod("cbind", "SpatialExperiment", function(..., deparse.level = 1) {
  
    args <- list(...)
    
    # Bind SpatialExperiment objects
    out <- do.call(
        methods::callNextMethod, 
        c(args, list(deparse.level=1)))
    
    # Merge 'imgData'
    if (!is.null(imgData(args[[1]]))) { 
        newimgdata <- do.call(rbind, lapply(args, SpatialExperiment::imgData))
        SingleCellExperiment::int_metadata(out)[names(SingleCellExperiment::int_metadata(out)) == 
                                                "imgData"] <- NULL
        SingleCellExperiment::int_metadata(out)$imgData <- newimgdata
    } 
    return(out)
})

#' get abundance long
#'
#' @keywords internal
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom purrr when
#' @importFrom purrr map2
#' @importFrom SummarizedExperiment assays
#'
#' @param .data A tidySpatialExperiment
#' @param features A character
#' @param all A boolean
#' @param exclude_zeros A boolean
#'
#' @return A tidySpatialExperiment object
#'
#' @noRd
get_abundance_sc_long <- function(.data, features = NULL, all = FALSE, exclude_zeros = FALSE) {

    # Define unbound variable
    .feature <- NULL
    
    # For SCE there is not filed for variable features
    variable_feature <- c()

    # Check if output would be too big without forcing
    if (
        length(variable_feature) == 0 &
            is.null(features) &
            all == FALSE
    ) {
        stop(
            "Your object does not contain variable feature labels,
            feature argument is empty and all arguments are set to FALSE.
            Either:
            1. use detect_variable_features() to select variable feature
            2. pass an array of feature names
            3. set all=TRUE (this will output a very large object, does your computer have enough 
            RAM?)"
        )
    }

    # Get variable features if existing
    if (
        length(variable_feature) > 0 &
            is.null(features) &
            all == FALSE
    ) {
        variable_genes <- variable_feature
    } 
    else {
        variable_genes <- NULL
    }

    assay_names <- assays(.data) |> names()

    # Check that I have assay manes
    if(length(assay_names) == 0)
        stop(
            "tidySpatialExperiment says: there are no assays names in the source SpatialExperiment."
        )

    assays(.data) |>
        as.list() |>
        
        # Take active assay
        map2(assay_names, function(x, y) {
            x <- 
                x |>
                    when(
                        variable_genes |> 
                            is.null() |> 
                            lapply(`!`) ~ 
                            x[variable_genes, , drop = FALSE],
                        features |>
                            is.null() |> 
                            lapply(`!`) |> 
                            unlist() ~ 
                            x[toupper(rownames(x)) %in% toupper(features), , drop = FALSE],
                        all ~ x,
                        ~ stop(
                            "It is not convenient to extract all genes, you should have either 
                            variable features or feature list to extract"
                        )
                    ) |>
                    as.matrix() |>
                    DataFrame() |>
                    tibble::as_tibble(rownames = ".feature")
            
            # Rename columns with index to handle duplicate cell names
            colnames(x) <- c(".feature", 2:length(colnames(x))-1)
            
            x |>
                tidyr::pivot_longer(
                    cols = -.feature,
                    names_to = "index",
                    values_to = ".abundance" |> paste(y, sep = "_"),
                    values_drop_na = TRUE
                )
        }) |>
        reduce(function(...) full_join(..., by = c(".feature", c_(.data)$name)))
}

#' @importFrom S4Vectors DataFrame
#'
#' @keywords internal
#'
#' @param .data A tibble
#' @param SpatialExperiment_object A tidySpatialExperiment
#'
#' @noRd
as_meta_data <- function(.data, SpatialExperiment_object) {

    col_to_exclude <-

        # special_datasets_to_tibble(SpatialExperiment_object) |>
        get_special_columns(SpatialExperiment_object) |>
  
        # I need this in case we have multiple reduced dimension data frames with overlapping names 
        # of the columns. For example multiple PCA versions
        vctrs::vec_as_names(repair = "unique") |>

    # To avoid name change by the bind_cols of as_tibble
    trick_to_avoid_renaming_of_already_unique_columns_by_dplyr()
    
    col_names <- 
        colnames(.data)[!colnames(.data) %in% col_to_exclude]
    
    .data_df <-
        .data |>
        select(col_names) |>
        data.frame()
    
    # Select row names and change to dataframe class, allowing for duplicate values
    row_names <- 
        .data_df |> 
        pull(!!c_(SpatialExperiment_object)$symbol)
    
    .data_df <-
        .data_df |>
        DataFrame()
    
    # Set new rownames and remove column of origin
    rownames(.data_df) <- row_names
    .data_df <- .data_df[, !names(.data_df) == c_(SpatialExperiment_object)$symbol, drop = FALSE]
    
    .data_df
}

#' @importFrom purrr map_chr
#'
#' @keywords internal
#'
#' @param SpatialExperiment_object A tidySpatialExperiment
#'
#' @noRd
#'
get_special_columns <- function(SpatialExperiment_object) {
    SpatialExperiment_object |>
        get_special_datasets() |>
        map(~ .x |> colnames()) |>
        unlist() |>
        as.character()
}

#' @importFrom SpatialExperiment spatialCoords
#'
#' @keywords internal
#' @noRd
get_special_datasets <- function(SpatialExperiment_object, n_dimensions_to_return = Inf) {
    rd <- SpatialExperiment_object@int_colData@listData$reducedDims
    
    reduced_dimensions <-
        map2(rd |> as.list(), names(rd), ~ {
            mat <- .x[, seq_len(min(n_dimensions_to_return, ncol(.x))), drop = FALSE]
            
            # Set names as SCE is much less constrained and there could be missing names
            if (length(colnames(mat)) == 0) colnames(mat) <- sprintf("%s%s", .y, seq_len(ncol(mat)))
            mat
      })
    
    spatial_coordinates <- 
        SpatialExperiment::spatialCoords(SpatialExperiment_object)
    
    list(reduced_dimensions, spatial_coordinates)
}

#' @importFrom dplyr select
#' @importFrom rlang expr
#' @importFrom tidyselect eval_select
#' @importFrom tidyselect all_of
#' 
#' @keywords internal
#' @noRd
select_helper <- function(.data, ...) {
  
  loc <- tidyselect::eval_select(rlang::expr(c(...)), .data)
  dplyr::select(.data, tidyselect::all_of(loc))
}

# Key column names
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
ping_old_special_column_into_metadata <- 
    tidySingleCellExperiment:::ping_old_special_column_into_metadata

#' @importFrom S4Vectors metadata
c_ <-  tidySingleCellExperiment:::c_

# Greater than
gt <- function(a, b) {
  a > b
}

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr <- tidySingleCellExperiment:::add_attr

#' Subset columns
#'
#' @keywords internal
#' @noRd
#' 
#' @importFrom rlang enquo
#' 
#' @param .data A tibble
#' @param .column A vector of column names
#'
#' @return A tibble
subset <- tidySingleCellExperiment:::subset

#' Add class to abject
#'
#'
#' @keywords internal
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class <- tidySingleCellExperiment:::add_class

#' Remove class to abject
#'
#' @keywords internal
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
#' @keywords internal
drop_class <- tidySingleCellExperiment:::drop_class

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- tidySingleCellExperiment:::quo_names

#' as_SummarizedExperiment
#'
#' @keywords internal
#' @noRd
#' 
#' @description as_SummarizedExperiment creates a `SummarizedExperiment` object from a `tbl` or 
#' `tidybulk` tbl formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#' @importFrom utils data
#' @importFrom tidyr pivot_longer
#'
#' @param .data A tibble
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @return A `SummarizedExperiment` object
as_SummarizedExperiment <- tidySingleCellExperiment:::as_SummarizedExperiment

# This function is used for the change of special sample column to .sample
# Check if "sample" is included in the query and is not part of any other existing annotation
#' @importFrom stringr str_detect
#' @importFrom stringr regex
is_sample_feature_deprecated_used <- tidySingleCellExperiment:::is_sample_feature_deprecated_used

#' @importFrom stringr str_replace_all
trick_to_avoid_renaming_of_already_unique_columns_by_dplyr <- 
    tidySingleCellExperiment:::trick_to_avoid_renaming_of_already_unique_columns_by_dplyr

# Inherited functions with no documentation
get_special_column_name_cell <- tidySingleCellExperiment:::get_special_column_name_cell
get_needed_columns <- tidySingleCellExperiment:::get_needed_columns
get_special_column_name_symbol <- tidySingleCellExperiment:::get_special_column_name_symbol
special_datasets_to_tibble <- tidySingleCellExperiment:::special_datasets_to_tibble
prepend <- tidySingleCellExperiment:::prepend

# Define messages
data_frame_returned_message <- 
    "tidySpatialExperiment says: A data frame is returned for independent data analysis."
duplicated_cell_names <- 
    "tidySpatialExperiment says: This operation lead to duplicated cell names. A data frame is 
    returned for independent data analysis."
