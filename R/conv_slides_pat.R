#' Process tibble with slide level info and sum to get patient level info
#'
#' @export
#' @param full_part_all tibble with slide data.
#' @return tibble with patient level data
#'
convertSlidesToPatient = function(full_part_all){

  if(!("tot_mut_crypts" %in% colnames(full_part_all))) full_part_all = full_part_all %>% dplyr::mutate(tot_mut_crypts = NA)

  # 1. Filter so that 1 slide per block
  filterSlidesPerBlock = full_part_all    %>%
    dplyr::group_by(Block.ID)                    %>%
    dplyr::filter(num_crypts == max(num_crypts)) %>% dplyr::ungroup()
  # 2. Collapse crypts and counts for blocks from same patient
  PerPat_data = filterSlidesPerBlock %>%
    dplyr::group_by(DW.Patient.number, Age, SEX) %>%
    dplyr::summarise( num_crypts = sum(num_crypts, na.rm = T),
               partial    = sum(partial, na.rm = T),
               monoclonal = sum(monoclonal, na.rm = T), full = sum(full, na.rm = T),
               fullx2 = sum(fullx2, na.rm = T), fullx3 = sum(fullx3, na.rm = T), fullx4  = sum(fullx4, na.rm = T),
               fullx5 = sum(fullx5, na.rm = T), fullx6 = sum(fullx6, na.rm = T), fullx7  = sum(fullx7, na.rm = T),
               fullx8 = sum(fullx8, na.rm = T), fullx9 = sum(fullx9, na.rm = T), fullx10 = sum(fullx10, na.rm = T),
               fullxmore10 = sum(fullxmore10, na.rm = T),
               tot_mut_crypts = sum(tot_mut_crypts, na.rm = T)) %>%
    dplyr::ungroup()
  # 3. Remove cases without crypts or patient
  PerPat_data = PerPat_data %>% dplyr::filter(!is.na(num_crypts), !is.na(Age))
  PerPat_data
}
