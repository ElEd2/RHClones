# library(tidyverse)

#' Read xls file and process
#'
#' @export
#' @param file_name Name of xls file.
#'
read_raw_data <- function(file_name){
  # file_name = "counts/FMR1_data.xlsx"
  # file_out  = "FMR1_data_proc.Rdata"
  # file_out_all  = "FMR1_data_all_patch.csv"
  raw_file  = readxl::read_xlsx(file_name)

  # Sep col, rename a few, remove repeat entries
  proc_file = raw_file %>%
    tidyr::separate(`Achiever ID`, c("Block.ID", "Section")) %>%
    dplyr::rename(monoclonal        = `Total WPC`, partial = PPC, num_crypts = NCrypts,
           SEX               = Sex,
           DW.Patient.number = `Tissue Bank ID`,
           Patch_amount      = `Patch Amount`,
           Patch_size        = `Patch Size`) %>%
    dplyr::distinct()

  # Filter on QC
  proc_file = proc_file %>%
    dplyr::filter(`Stain Quality` >= 1.5) %>%
    dplyr::mutate(remove_me = stringr::str_detect(Comments, stringr::fixed("DO NOT", ignore_case = T)),
           remove_me = dplyr::if_else(is.na(remove_me), F, remove_me)) %>%
    dplyr::filter(!remove_me) %>%
    dplyr::select(-remove_me)

  # Total mutant crypts
  proc_file = proc_file %>%
    dplyr::group_by(Block.ID, Section) %>%
    dplyr::mutate(tot_mut = as.numeric(Patch_size%*%Patch_amount)) %>%
    dplyr::ungroup()

  # proc_file %>% write_csv(file_out_all)
  proc_file_patch = proc_file

  # create 10 or more category, Spread pacthes, get rid of zeros
  patch_sizes_str =  paste0("fullx", 0:11)
  proc_file = proc_file %>%  dplyr::mutate(Patch_size = dplyr::if_else(Patch_size > 10, 11, Patch_size)) %>%
    dplyr::mutate(Patch_size = factor(paste( paste0("fullx", Patch_size)), levels = patch_sizes_str))

  # If multiple larger than 10 patches need to merge
  proc_file = proc_file %>% dplyr::group_by(dplyr::across(-Patch_amount)) %>%
    dplyr::summarise(Patch_amount = sum(Patch_amount)) %>% dplyr::ungroup()

  # dodgy workaround
  dummy_tbl = dplyr::tibble(Block.ID = "Ed_added", Patch_size = factor(patch_sizes_str, levels = patch_sizes_str), Patch_amount = 0)

  proc_file = proc_file %>% dplyr::bind_rows(dummy_tbl) %>% tidyr::spread(Patch_size, Patch_amount, fill = 0) %>%
    dplyr::filter(Block.ID != "Ed_added") %>% dplyr::select(-fullx0) %>% dplyr::rename(fullxmore10 = fullx11, full = fullx1)

  # Choose a slide if multiple slides per-block
  proc_file = proc_file %>% dplyr::group_by(Block.ID) %>% dplyr::top_n(1, `Stain Quality`) %>% dplyr::ungroup()

  num_pat    = proc_file %>% dplyr::pull(DW.Patient.number) %>% dplyr::n_distinct()
  num_slides = nrow(proc_file)
  print(paste0("Slides =  ", num_slides, " Patients = ", num_pat))
  negative_full_clones = proc_file %>% dplyr::filter(full < 0)
  print(paste0("Num neg full error ", nrow(negative_full_clones)))

  list(slides = proc_file, patch = proc_file_patch)
  # save(proc_file, file = file_out)
}
# print(out_vals$missing_pat %>% select(Section, Block.ID, `stain quality`) %>% as.data.frame())
# print(out_vals$clone_counts %>% filter(is.na(num_crypts)) %>% select(Section, Block.ID, `stain quality`) %>% as.data.frame())

