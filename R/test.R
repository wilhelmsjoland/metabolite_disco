#  C15H10O5 apigenin
#  C15H10O6 luteolin


Rdisop::getMass(Rdisop::getMolecule(aglycone_form))

try <- hits %>%
  dplyr::filter(chem_change == "± O")

try$feat1
try$feat2

try[1,]

possible.adducts %>%
  dplyr::filter(feature == "FT02089")

anno %>%
  dplyr::filter(peak_id == "FT03779") %>%
  view()

possible.adducts %>%
  dplyr::filter(feature == "FT03779")

test <- predictBiotransfAdductsSubset(
  data = possible.adducts,
  biotransf.data = bio.transf2, # bio.transf
  tolerance_ppm = 1000, # glycoside_ppm
  feat_filt = c(pot.glycosides)
)

test2 <- test %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = all_of(c("feat1", "feat2")),
      .fns = ~ .x %in% c("FT03779", "FT02089")
    )
  )

# calc obs ppm
test2 %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = dplyr::all_of(c("feat1", "feat2")),
      .fns = ~ .x %in% c("FT03779", "FT02089")
    )
  ) %>%
  dplyr::filter(feat1 != feat2) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(obs_ppm = num_to_ppm(abs(delta_mass - obs_delta_mass))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(obs_ppm) %>%
  dplyr::select(-name) %>%
  view()

#  C15H10O5 apigenin
#  C15H10O6 luteolin

x <- getMass(getMolecule("C15H10O5"))
y <- getMass(getMolecule("C15H10O6"))

num_to_ppm(abs(x - y))

test %>%
  dplyr::filter(
    dplyr::if_all(
      .cols = all_of(c("feat1", "feat2")),
      .fns = ~ .x %in% c("FT02089", "FT03779")
    )
  ) %>%
  view()

# MAKE SURE I CALCULATE THE NEUTRAL MASSES FOR EVERY ADDUCT FOR EVERY FEATURE

xx <- possible.adducts %>%
  dplyr::filter(feature == "FT02089") %>%
  dplyr::mutate(group = paste0(feature, "_", adduct))

yy <- possible.adducts %>%
  dplyr::filter(feature == "FT03779") %>%
  dplyr::mutate(group = paste0(feature, "_", adduct))

tibble(
  api = xx$group,
  lut = yy$group
) %>%
combn()


test.ppm <- getMass(getMolecule("O")) + c(-ppm_to_num(5000), +ppm_to_num(5000))

try <- combn(
  c(xx$group, yy$group),
  m = 2
) %>%
t() %>%
tibble::as_tibble(.) %>%
dplyr::mutate(
  V1_name = stringr::str_split_i(V1, "_", 1),
  V1_adduct = stringr::str_split_i(V1, "_", 2),
  V2_name = stringr::str_split_i(V2, "_", 1),
  V2_adduct = stringr::str_split_i(V2, "_", 2)
) %>%
dplyr::filter(V1_name != V2_name) %>%
dplyr::left_join(
  x = .,
  y = xx %>%
    dplyr::select(feature, adduct, mass),
  by = c(
    "V1_name" = "feature",
    "V1_adduct" = "adduct"
    )
) %>%
dplyr::left_join(
  x = .,
  y = yy %>%
    dplyr::select(feature, adduct, mass),
  by = c(
    "V2_name" = "feature",
    "V2_adduct" = "adduct"
    )
) %>%
rowwise() %>%
dplyr::mutate(abs_diff = abs(mass.x - mass.y)) %>%
ungroup() 

# 43 & 62
num_to_ppm(abs(try[62,]$abs_diff - getMass(getMolecule("O"))))

try %>%
  dplyr::filter(between(abs_diff, test.ppm[1], test.ppm[2]))

anno %>%
  dplyr::filter(peak_id == "FT03779") %>%
  dplyr::filter(grepl("Lut", target_name)) %>% view

test <- mass2mz(getMass(getMolecule("C15H10O6")), adduct = adducts(polarity = "negative")) %>% t() %>% as_tibble(., rownames ="adduct")
test

lut.adduct <- mass2mz(getMass(getMolecule("C15H10O6")), adduct = "[M+C2H3O2]-")

try <- mass2mz(getMass(getMolecule("C15H10O5")), adduct = adducts(polarity = "negative")) %>% 
t() %>%
tibble::as_tibble(., rownames = "adduct") %>%
rowwise() %>%
dplyr::mutate(diff = abs(lut.adduct[1] - V1))

difference.is <- try[9,]$diff

abs(difference.is - theory.difference.is)

getMass(getMolecule("O")) + 0.1

# real diffs
test <- try %>%
  dplyr::filter(between(abs_diff, getMass(getMolecule("O")) - 0.1, getMass(getMolecule("O")) + 0.1))

# theory diff
theory.difference.is <- abs(getMass(getMolecule("C15H10O5")) - getMass(getMolecule("C15H10O6")))

num_to_ppm(abs(test[1,]$abs_diff - theory.difference.is))

anno.test <- anno %>%
  dplyr::filter(peak_id == "FT03779") %>%
  dplyr::filter(grepl("Lut", target_name))

num_to_ppm(abs(anno.test[1,][["mz"]] - lut.adduct))

test <- target_df %>% 
  as_tibble(., rownames = "id") %>%
  dplyr::filter(name == "Luteolin")

# the difference between massbank and expected ppm is between 0.02 to 191262
min(num_to_ppm(abs(test[["exactmass"]] - getMass(getMolecule("C15H10O6")))))
max(num_to_ppm(abs(test[["exactmass"]] - getMass(getMolecule("C15H10O6")))))

# luteolin mass
getMass(getMolecule("C15H10O6"))

# differences in mass vs massbank in ppm (0.02 to 191262)
min(num_to_ppm(abs(test[["exactmass"]] - getMass(getMolecule("C15H10O6")))))
max(num_to_ppm(abs(test[["exactmass"]] - getMass(getMolecule("C15H10O6")))))

lut.observed.mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature == "FT03779") %>%
  dplyr::pull(mzmed)

my.lut.mass <- mass2mz(
  getMass(getMolecule("C15H10O6")),
  adduct = "[M+C2H3O2]-"
  ) %>%
  .[1]
  
api.observed.mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature == "FT02089") %>%
  dplyr::pull(mzmed)

my.api.mass <- mass2mz(
  getMass(getMolecule("C15H10O5")),
  adduct = "[M+C2H3O2]-"
  ) %>%
  .[1]

lut.observed.mz
api.observed.mz

# raw m/z
abs(api.observed.mz - lut.observed.mz)

# known apigenin mass vs luteolin mass calculated via mz and adduct
abs(my.api.mass - my.lut.mass)

abs(getMass(getMolecule("C15H10O6")) - my.lut.mass)

num_to_ppm(abs(calculateMass("C15H10O6") - getMass(getMolecule("C15H10O6"))))

formula2mz("C15H10O6", adduct = "[M+C2H3O2]-")



theory.api.mass <- calculateMass("C15H10O5")[[1]]
theory.lut.mass <- calculateMass("C15H10O6")[[1]]

abs(theory.api.mass - theory.lut.mass)





observed.api.mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature == "FT02089") %>%
  dplyr::pull(mzmed)

observed.lut.mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature == "FT03779") %>%
  dplyr::pull(mzmed)

observed.api.mz
observed.lut.mz

# apigenin
xx <- mz2mass(observed.api.mz, adduct = "[M-H]-")[[1]]
# luteolin
yy <- mz2mass(observed.lut.mz, adduct = "[M+C2H3O2]-")[[1]]

observed.diff <- abs(xx - yy)

num_to_ppm(abs(calculateMass("O") - observed.diff))

observed.api.mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature == "FT02089") %>%
  dplyr::pull(mzmed)

observed.lut.mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature == "FT03779") %>%
  dplyr::pull(mzmed)

api <- mz2mass(
  observed.api.mz,
  adduct = adducts(polarity = "negative")
  ) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("api_mass" = V1) %>%
  dplyr::mutate(adduct = paste0("apigenin_", adduct))

lut <- mz2mass(
  observed.lut.mz,
  adduct = adducts(polarity = "negative")
  ) %>%
  t() %>%
  tibble::as_tibble(., rownames = "adduct") %>%
  dplyr::rename("lut_mass" = V1) %>%
  dplyr::mutate(adduct = paste0("luteolin_", adduct))

api
lut

test <- combn(
  c(api$adduct, lut$adduct),
  m = 2
  ) %>%
  t() %>%
  tibble::as_tibble(.) %>%
  dplyr::mutate(
    V1_name = stringr::str_split_i(V1, "_", 1),
    V1_adduct = stringr::str_split_i(V1, "_", 2),
    V2_name = stringr::str_split_i(V2, "_", 1),
    V2_adduct = stringr::str_split_i(V2, "_", 2)
  ) %>%
  dplyr::filter(V1_name != V2_name) %>%
  dplyr::select(-c("V1_name", "V2_name", "V1_adduct", "V2_adduct")) %>%
  dplyr::left_join(
    x = .,
    y = api,
    by = c("V1" = "adduct")
  ) %>%
  dplyr::left_join(
    x = .,
    y = lut,
    by = c("V2" = "adduct")
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    diff = abs(api_mass - lut_mass)
  ) %>%
  dplyr::filter(between(diff, 14.5, 16.8))

exp.diff <- abs(calculateMass("C15H10O5")[[1]] - calculateMass("C15H10O6")[[1]])

test
num_to_ppm(abs(test[1,][["diff"]] - exp.diff))
num_to_ppm(abs(test[2,][["diff"]] - exp.diff))
num_to_ppm(abs(test[3,][["diff"]] - exp.diff))

num_to_ppm(abs(test[1,][["diff"]] - exp.diff))
num_to_ppm(abs(test[2,][["diff"]] - exp.diff))
num_to_ppm(abs(test[3,][["diff"]] - exp.diff))

try <- abs(test[1,][["diff"]] - calculateMass("O"))
num_to_ppm(try[[1]])

massbank.luts <- target_df %>%
  as_tibble(.) %>%
  dplyr::filter(name == "Luteolin") %>%
  pull(exactmass)

for (i in massbank.luts) {
  tmp <- abs(calculateMass("C15H10O6")[[1]] - i)
  tmp2 <- num_to_ppm(tmp)
  print(tmp2)
}

target_df %>%
  as_tibble(.) %>% 
  dplyr::filter(name == "Luteolin") %>% view

anno %>%
  dplyr::filter(target_name == "Luteolin") %>%
  view

num_to_ppm(abs(calculateMass("C15H10O6")[[1]] - 286.04773803))

featureDefinitions(xchr9) %>%
  as_tibble(., rownames = "feature") %>% 
  dplyr::filter(
    between(
      mzmed,
      c(calculateMass("C15H10O6")[[1]] - ppm_to_num(15)),
      c(calculateMass("C15H1006")[[1]] + ppm_to_num(15))
    ))

calculateMass("C15H10O6")[[1]] + c(+ppm_to_num(100), -ppm_to_num(100))
