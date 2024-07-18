library(tidyverse)
library(tidymodels)
library(ggbeeswarm)
library(ggtext)
library(conflicted)
library(themis)
library(patchwork)
library(here)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

# functions ---------------------------------------------------------------

calc_chord <- function(df) {
  
  df <- df %>% 
    # group_by(plate, well, particle) %>%
    mutate(min = min(frame),
           max = max(frame)) %>% 
    filter(frame == min | frame == max) %>% 
    arrange(frame) %>% 
    mutate(chord = sqrt((lead(x) - x)^2 + (lead(y) - y)^2)) %>% 
    drop_na(chord)
  
  return(df$chord[1])
  
}


# import ------------------------------------------------------------------

flow_files <- tibble(file = list.files(here(), "*_flow.csv", recursive = TRUE)) %>% 
  mutate(experiment_date = str_extract(file, '202[0-9]{5}'),
         plate = str_extract(file, "202[0-9]{5}-p[0-9]{2}-[A-Z]{2,3}"))

tracking_files <- tibble(file = list.files(here(), "*_tracking.csv", recursive = TRUE))


# tidy --------------------------------------------------------------------
                     
flow_data <- flow_files %>% 
  mutate(data = map(file, read_csv)) %>% 
  unnest(data) %>% 
  mutate(normalized_flow = optical_flow / worm_area) %>% 
  pivot_longer(cols = optical_flow:normalized_flow) %>% 
  mutate(conc = case_when(
    conc == '0.1p' ~ 0.0001,
    TRUE ~ as.numeric(str_remove(conc, 'uM'))
  ))

metadata <- flow_data %>% 
  select(experiment_date:timepoint) %>% 
  distinct()

tracking_data <- tracking_files %>% 
  mutate(df = map(file, read_csv)) %>% 
  unnest(cols = c(df)) %>% 
  mutate(
    plate = str_extract(file, "202[0-9]{5}-p[0-9]{2}-[A-Z]{2,3}"),
    plate_number = str_extract(plate, 'p[0-9]{2}')
  ) %>% 
  arrange(plate, well, particle, frame) %>% 
  nest(.by = c(plate, plate_number, well, particle)) %>% 
  mutate(chord = unlist(map(data, calc_chord))) %>% 
  unnest(c(data)) %>% 
  group_by(plate, plate_number, well, particle) %>%
  mutate(
    diff_x = c(0, diff(x)),
    diff_y = c(0, diff(y)),
    change_x = case_when(
      sign(diff_x) != sign(lead(diff_x)) & diff_x != 0 ~ TRUE,
      TRUE ~ FALSE),
    change_y = case_when(
      sign(diff_y) != sign(lead(diff_y)) & diff_y != 0 ~ TRUE,
      TRUE ~ FALSE),
    dist = sqrt((lead(x) - x)^2 + (lead(y) - y)^2),
    w = atan2(lead(y) - y, lead(x) - x) - atan2(y - lag(y), x - lag(x))
  ) %>% 
  summarise(
    n_frames = n(),
    total_distance = sum(dist, na.rm = TRUE),
    total_changes = sum(change_x + change_y),
    time = (max(frame) - min(frame)) / 16,
    velocity = total_distance / time,
    rcd = total_changes / time, 
    mean_w = mean(w, na.rm = TRUE),
    sd_w = sd(w, na.rm = TRUE),
    chord = mean(chord),
    ave_mass = mean(mass, na.rm = TRUE),
    ave_size = mean(size, na.rm = TRUE),
    ave_ecc = mean(ecc, na.rm = TRUE),
    ave_signal = mean(signal, na.rm = TRUE),
    ave_rawmass = mean(raw_mass, na.rm = TRUE),
    # ave_ep = mean(ep, na.rm = TRUE)
  ) %>%
  mutate(
    arc_chord_ratio = total_distance / chord
  ) %>% 
  drop_na()


# build model -------------------------------------------------------------

features <- read_rds(here('classification', 'feature_summary.rds'))

classes <- read_csv(here('classification', 'particle_annotations.csv')) %>% 
  left_join(features)

classes %>% 
  pivot_longer(where(is.numeric)) %>% 
  ggplot() +
  geom_quasirandom(aes(x = 1, y = value, color = class)) +
  facet_wrap(vars(name), scales = 'free_y') +
  theme_minimal()

model_data <- classes %>% 
  select(-plate, -well, -plate_number) %>% 
  mutate(class = factor(class)) %>% 
  drop_na()

set.seed(123)
data_split <- initial_split(model_data,
                            strata = class)
train_data <- training(data_split)
test_data <- testing(data_split)

set.seed(234)
folds <- vfold_cv(train_data,
                  v = 10,
                  strata = class)

cores <- parallel::detectCores()
rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune()) %>%
  set_engine('ranger', num.threads = cores) %>%
  set_mode('classification')

recipe <-
  recipe(class ~ ., data = model_data) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors()) %>%
  step_corr(all_numeric_predictors(), threshold = .5) %>%
  step_smote(class)

prep <- prep(recipe)
juice <- juice(prep)

rf_workflow <-
  workflow() %>%
  add_model(rand_forest_ranger_spec) %>%
  add_recipe(recipe)

rf_grid <- grid_regular(finalize(mtry(), model_data),
                        min_n())

rf_tune <-
  rf_workflow %>%
  tune_grid(
    resamples = folds,
    grid = rf_grid,
    control = control_grid(save_pred = TRUE,
                           verbose = TRUE),
    metrics = metric_set(roc_auc, sens))

# extract the best decision model
best_rf <- rf_tune %>%
  select_best(metric = "roc_auc")

# print metrics
(rf_metrics <- rf_tune %>% 
    collect_metrics() %>% 
    semi_join(best_rf) %>% 
    select(.metric:.config) %>% 
    mutate(model = 'Random forest'))

# finalize the wf with the best model
rf_workflow <-
  rf_workflow %>%
  finalize_workflow(best_rf)

# generate predictions on the hold-out test data
rf_auc <-
  rf_tune %>%
  collect_predictions(parameters = best_rf) %>%
  roc_curve(.pred_bad, truth = class) %>%
  mutate(model = "Random forest")

rf_auc %>%
  autoplot()

mtry <- best_rf$mtry
min_n <- best_rf$min_n

last_mod <-
  rand_forest(mtry = mtry,
              min_n = min_n) %>%
  set_engine('ranger', num.threads = cores, importance = "impurity") %>%
  set_mode('classification')

last_workflow <-
  rf_workflow %>%
  update_model(last_mod)

set.seed(345)
last_fit <-
  last_workflow %>%
  last_fit(data_split, 
           metrics = metric_set(roc_auc, sens))

collect_metrics(last_fit)

last_fit %>%
  extract_fit_engine() %>% 
  vip::vip() +
  theme_minimal()

(final_auc <-
    last_fit %>%
    collect_predictions() %>%
    roc_curve(.pred_bad, truth = class) %>%
    autoplot())

last_fit %>%
  collect_predictions() %>%
  conf_mat(truth = class, estimate = .pred_class) %>%
  autoplot()

final_wf <- last_fit %>% 
  extract_workflow()


# filter ------------------------------------------------------------------

post_filter <- augment(final_wf, tracking_data) %>% 
  filter(.pred_class == 'good',
         time > 5
  )

tracking_well <- post_filter %>% 
  select(-contains('pred')) %>% 
  select(plate:particle, sd_w, velocity, arc_chord_ratio, rcd) %>% 
  pivot_longer(c(sd_w, velocity, arc_chord_ratio, rcd)) %>% 
  group_by(plate, well, name) %>% 
  summarise(value = median(value, na.rm = TRUE)) %>% 
  left_join(metadata)

# plot --------------------------------------------------------------------

moving_worms <- post_filter %>% 
  left_join(select(metadata, -plate)) %>% 
  group_by(species, stage, treatment, conc, timepoint, plate, well) %>% 
  tally() %>% 
  rename(moving_worms = n)

final_df <- post_filter %>% 
  left_join(metadata) %>% 
  mutate(replicate = case_when(
    str_detect(well, 'A') ~ '1',
    str_detect(well, 'B') ~ '2',
    str_detect(well, 'C') ~ '3'
  )) %>% 
  select(-mean_w, -total_changes) %>% 
  pivot_longer(c(total_distance, sd_w, velocity, arc_chord_ratio, rcd)) %>% 
  mutate(
    value = case_when(
      name %in% c('arc_chord_ratio', 'total_distance', 'velocity') ~ log10(value),
      name == 'rcd' ~ log10(value + 1),
      TRUE ~ value
    ),
    name = case_when(
      name == 'arc_chord_ratio' ~ 'Log<sub>10</sub>(Tortuosity)',
      name == 'total_distance' ~ 'Log<sub>10</sub>(Distance)',
      name == 'velocity' ~ 'Log<sub>10</sub>(Velocity)',
      name == 'sd_w' ~ 'Angular<br>Velocity (SD)',
      name == 'rcd' ~ 'Log<sub>10</sub>(RCD + 1)',
    )
  ) %>% 
  group_by(plate, well, species, stage, treatment, conc, timepoint, name) %>%
  summarise(well_mean = mean(value)) %>%
  bind_rows(select(flow_data, -file, -experiment_date, well_mean = value)) %>% 
  bind_rows(select(moving_worms, well_mean = moving_worms) %>% mutate(name = 'Number of Tracks')) %>% 
  filter(!name %in% c('worm_area', 'normalized_flow')) %>% 
  # different fps on the 4/18 run
  mutate(
    well_mean = case_when(
      name == 'optical_flow' & plate == '20240418-p01-RVH' ~ well_mean * 30,
      TRUE ~ well_mean),
    name = case_when(
      name == 'optical_flow' ~ 'Total Motility',
      TRUE ~ name
    )) 

stats <- final_df %>% 
  group_by(name) %>%
  group_nest() %>% 
  mutate(
    aov = map(data, ~aov(.x$well_mean ~ .x$conc)),
    tidy = map(aov, broom::tidy),
  ) %>% 
  unnest(cols = tidy) %>% 
  drop_na() %>% 
  select(name, sumsq:p.value)

labels <- c(
  'Log<sub>10</sub>(Tortuosity)' = 'Log<sub>10</sub>(Tortuosity)<br>p = 0.0006',
  'Log<sub>10</sub>(Distance)' = 'Log<sub>10</sub>(Distance)<br>p = 0.002',
  'Log<sub>10</sub>(Velocity)' = 'Log<sub>10</sub>(Velocity)<br>p = 0.0002',
  'Angular<br>Velocity (SD)' = 'Angular<br>Velocity (SD)<br>p = 0.002',
  'Log<sub>10</sub>(RCD + 1)' = 'Log<sub>10</sub>(RCD + 1)<br>p = 0.0008',
  'Number of Tracks' = 'Number of Tracks<br>p = 0.054',
  'Total Motility' = 'Total Motility<br>p = 0.097'
)

(track_plot <-  final_df %>% 
    ggplot() +
    ggbeeswarm::geom_quasirandom(aes(x = conc, y = well_mean,),
                                 size = 0.5, alpha = 0.5
    ) +
    stat_summary(aes(x = conc, y = well_mean),
                 color = 'indianred',
                 size = 0.25, alpha = 1) +
    facet_wrap(facets = vars(name),
               scales = 'free_y',
               labeller = as_labeller(labels)) +
    scale_x_log10(
      breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10),
      labels = c('∅', '0.001', '0.01', '0.1', '1', '10')) +
    labs(x = 'PZQ concentration (µM)',
         y = 'Value',
         color = 'Replicate') +
    cowplot::theme_half_open() +
    theme(
      strip.text = ggtext::element_markdown(size = 9),
      strip.background = element_rect(fill = 'white'),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9)
    ) +
    NULL)


# example tracks ----------------------------------------------------------

ex_tracks <- tracking_files %>% 
  filter(str_detect(file, '20240307-p01-RVH')) %>% 
  mutate(df = map(file, read_csv)) %>% 
  unnest(cols = c(df)) %>% 
  mutate(
    plate = str_extract(file, "202[0-9]{5}-p[0-9]{2}-[A-Z]{2,3}"),
    plate_number = str_extract(plate, 'p[0-9]{2}')
  ) %>% 
  left_join(metadata) %>% 
  arrange(plate, well, particle, frame) %>% 
  filter(str_detect(well, '^A')) %>% 
  left_join(select(post_filter, .pred_class:particle)) %>% 
  filter(.pred_class == 'good')

(plot_tracks <- ex_tracks %>% 
    mutate(conc = case_when(
      conc == 1e-04 ~ "\U2205",
      conc == 0.001 ~ '0.001 µM',
      conc == 0.01 ~ '0.01 µM',
      conc == 0.1 ~ '0.1 µM',
      conc == 1 ~ '1 µM',
      conc == 10 ~ '10 µM',
    )) %>% 
    ggplot() +
    ggforce::geom_circle(
      data = tibble(x0 = 1616 / 2, y0 = 1616 / 2, r = 1616 / 2),
      aes(x0 = x0, y0 = y0, r = r), 
      inherit.aes = FALSE) +
    geom_path(
      aes(x = x, y = y, color = frame / 15, group = as.factor(particle))
    ) +
    facet_wrap(facets = vars(conc), ncol = 6, strip.position = "bottom") +
    scale_x_continuous(limits = c(0, 1616)) +
    scale_y_continuous(limits = c(0, 1616)) +
    labs(color = 'Time (s)') +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position.inside = c(1, 1),
      legend.key.size = unit(0.3, 'cm'),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6)
    ) +
    NULL)


# final plot --------------------------------------------------------------

final <- cowplot::plot_grid(plot_tracks, track_plot,
                            nrow = 2, rel_heights = c(.25, 1),
                            labels = 'AUTO', label_size = 12)

ggsave(here('pzq', 'Fig4.pdf'), final, width = 5, height = 6, device = cairo_pdf)
ggsave(here('pzq', 'Fig4.png'), final, width = 5, height = 6)
