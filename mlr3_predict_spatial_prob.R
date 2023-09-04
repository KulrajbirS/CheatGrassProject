predict_spatial_prob = function(newdata, learner, chunksize = 200L, format = "terra", filename = NULL) {
  library(checkmate)
  task = as_task_unsupervised(newdata)
  assert_multi_class(task$backend, c("DataBackendRaster", "DataBackendVector"))
  assert_learner(learner)
  
  if (test_class(task$backend, "DataBackendRaster")) {
    assert_number(chunksize)
    assert_choice(format, c("terra", "raster", "stars"))
    filename = filename %??% tempfile(fileext = ".tif")
    assert_path_for_output(filename)
    
    stack = task$backend$stack
    start_time = proc.time()[3]
    learner_start = learner$clone()
    learner = switch(learner$task_type,
                     "classif" = mlr3spatial:::LearnerClassifSpatial$new(learner),
                     "regr" = mlr3spatial:::LearnerRegrSpatial$new(learner))
    
    # calculate block size
    bs = block_size(stack, chunksize)
    
    # initialize target raster
    target_raster = terra::rast(terra::ext(stack), resolution = terra::res(stack), 
                                crs = terra::crs(stack), nlyrs = length(learner$learner$state$train_task$class_names))
    terra::writeStart(target_raster, filename = filename, overwrite = TRUE, datatype = "FLT8S",
                      names = learner$learner$state$train_task$class_names)
    
    mlr3:::lg$info("Start raster prediction")
    mlr3:::lg$info("Prediction is executed with a chunksize of %s Megabytes, %i chunk(s) in total, %i values per chunk",
            as.character(chunksize), length(bs$cells_seq), ceiling(terra::ncell(task$backend$stack) / length(bs$cells_seq)))
    
    mlr3misc::pmap(list(bs$cells_seq, bs$cells_to_read, seq_along(bs$cells_seq)), function(cells_seq, cells_to_read, n) {
      
      stack = task$backend$stack
      pred = learner$predict(task, row_ids = cells_seq:((cells_seq + cells_to_read - 1)))
      
      # The pred$prob object does not contain rows with NA values, those need
      # to be added in and that can be accomplished by inserting those values
      # into an empty vector at the non-NA positions gleaned from the responses.
      non_na_pos = which(!is.nan(pred$response))
      v = matrix(rep(NaN, ncol(pred$prob) * length(pred$response)), 
                 ncol = ncol(pred$prob), dimnames = list(c(), colnames(pred$prob)))
      for(i in 1:ncol(pred$prob)) {
        v[non_na_pos, i] = pred$prob[, i]
      }
      terra::writeValues(x = target_raster, v = v,
                         start = terra::rowFromCell(stack, cells_seq), # start row number
                         nrows = terra::rowFromCell(stack, cells_to_read)) # how many rows
      mlr3:::lg$info("Chunk %i of %i finished", n, length(bs$cells_seq))
    })
    
    terra::writeStop(target_raster)
    mlr3:::lg$info("Finished raster prediction in %i seconds", as.integer(proc.time()[3] - start_time))
    
    if (learner_start$task_type == "classif" && learner_start$predict_type == "response") {
      levels = learner$learner$state$train_task$levels()[[learner$learner$state$train_task$target_names]]
      value = data.table(ID = seq_along(levels), categories = levels)
      target_raster = terra::categories(target_raster, value = value, index = 2)
    }
    # target_raster = set_names(target_raster, learner$learner$state$train_task$target_names)
    
    switch(format,
           "terra" = target_raster,
           "stars" = stars::st_as_stars(target_raster),
           "raster" = as(target_raster, "Raster")
    )
  } else {
    assert_string(format, "sf")
    if (!is.null(filename)) assert_path_for_output(filename)
    pred = learner$predict(task)
    vector = set_names(sf::st_as_sf(data.frame(pred$response, task$backend$sfc)), c(learner$state$train_task$target_names, "geometry"))
    
    if (!is.null(filename)) sf::st_write(vector, filename, quiet = TRUE)
    vector
  }
}