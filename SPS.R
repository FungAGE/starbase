library(systemPipeShiny)
spsInit(app_path="/home/adrian/Systematics/Starship_Database",
  project_name="starbase",
  database_name="starbase.db",overwrite=T)
sps_tmp_dir <- tempdir()
spsInit(app_path = sps_tmp_dir, change_wd = FALSE, project_name = "SPSProject")
sps_dir <- file.path(sps_tmp_dir, "SPSProject")