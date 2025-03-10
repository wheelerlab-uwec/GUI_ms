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

quick_nest <- function(df, nest_cols) {
  nested_df <- vctrs::vec_split(
    x = df[setdiff(colnames(df), nest_cols)],
    by = df[nest_cols]
  )
  
  nested_df <- vctrs::vec_cbind(nested_df$key, new_tibble(list(data = nested_df$val)))
}

calc_features <- function(df) {
  print(str_c(df$row[[1]], df$col[[1]]))
  df <- df |>
    arrange(frame, .by_group = TRUE) |>
    mutate(
      diff_x = c(0, diff(x)),
      diff_y = c(0, diff(y)),
      change_x = case_when(
        sign(diff_x) != sign(lead(diff_x)) & diff_x != 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      change_y = case_when(
        sign(diff_y) != sign(lead(diff_y)) & diff_y != 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      total_changes = sum(change_x) + sum(change_y),
      dist = sqrt((lead(x) - x)^2 + (lead(y) - y)^2),
      w = atan2(lead(y) - y, lead(x) - x) - atan2(y - lag(y), x - lag(x)),
      sd_w = sd(w, na.rm = TRUE),
      total_dist = sum(dist, na.rm = TRUE)
    ) |>
    slice(c(1, n())) |>
    mutate(
      chord = sqrt((lead(x) - x)^2 + (lead(y) - y)^2),
      # 8 FPS
      time = (max(frame) - min(frame)) / 8,
      vel = total_dist / time,
      rcd = total_changes / time
    ) |>
    select(chord, total_dist, time, vel, sd_w, rcd) |>
    drop_na() |>
    distinct()
  
  return(df)
}


# import and tidy ---------------------------------------------------------

flow_data <- read_csv(here('pzq', 'data', '20250305-p01-NJW_flow.csv')) |> 
  select(-worm_area) |> 
  pivot_longer(cols = optical_flow) |> 
  select(batch, well, concentrations:value) |> 
  rename(well_mean = value)

tracking_data <- read_csv(here('pzq', 'data', '20250305-p01-NJW_tidy.csv')) |> 
  arrange(batch, well, particle, frame)

nest_cols <- c("batch", "well", "particle")
nested <- quick_nest(tracking_data, nest_cols) |> 
  mutate(
    features = map(data, calc_features),
    sd_x = map_dbl(data, ~ sd(.x$x, na.rm = TRUE)),
    sd_y = map_dbl(data, ~ sd(.x$y, na.rm = TRUE))
  )

write_rds(nested, here('pzq', 'data', 'nested.rds'))
nested <- read_rds(here('pzq', 'data', 'nested.rds'))

filtered <- nested |> 
  unnest(c(features)) |>
  unnest(c(data)) |> 
  group_by(batch, well, particle) |>
  arrange(frame, .by_group = TRUE) |>
  group_by(batch, well) |>
  filter(
    vel > 2.5,
    chord > 20,
    sd_x > 5 & sd_y > 5,
    signal > 5,
    # only those inside a circle/well
    (x - 1616 / 2)^2 + (y - 1616 / 2)^2 <= (1500 / 2)^2
  )

filtered_groups <- filtered |> 
  group_by(batch, well, particle) |> 
  summarise(max_mass = max(mass),
            min_mass = min(mass)) |> 
  mutate(diff_mass = max_mass - min_mass) |> 
  select(-max_mass, -min_mass)

more_filtered <- filtered |> 
  left_join(filtered_groups) |> 
  # filter by mass
  filter(diff_mass > 600) 

write_rds(more_filtered, here('pzq', 'data', 'more_filtered.rds'))


# plot --------------------------------------------------------------------

final_df <- more_filtered %>% 
  select(batch, well, treatments, concentrations, particle, chord:rcd) |> 
  distinct() |> 
  mutate(arc_chord_ratio = total_dist / chord) |> 
  select(-chord) |> 
  pivot_longer(c(total_dist, sd_w, vel, arc_chord_ratio, rcd)) %>% 
  mutate(
    value = case_when(
      name %in% c('arc_chord_ratio', 'total_dist', 'vel') ~ log10(value),
      name == 'rcd' ~ log10(value + 1),
      TRUE ~ value
    )
  ) %>% 
  group_by(batch, well, treatments, concentrations, name) %>%
  summarise(well_mean = mean(value)) |> 
  bind_rows(flow_data) |> 
  # different fps on the 4/18 run
  mutate(
    well_mean = case_when(
      name == 'optical_flow' & batch == '20240418' ~ well_mean * 30,
      TRUE ~ well_mean)) |> 
  mutate(
    name = case_when(
      name == 'arc_chord_ratio' ~ 'Log<sub>10</sub>(Tortuosity)',
      name == 'total_dist' ~ 'Log<sub>10</sub>(Distance)',
      name == 'vel' ~ 'Log<sub>10</sub>(Velocity)',
      name == 'sd_w' ~ 'Angular<br>Velocity (SD)',
      name == 'rcd' ~ 'Log<sub>10</sub>(RCD + 1)',
      name == 'optical_flow' ~ 'Total Motility'
    )
  ) 

stats <- final_df %>% 
  group_by(name) %>%
  group_nest() %>% 
  mutate(
    aov = map(data, ~aov(.x$well_mean ~ .x$concentrations)),
    tidy = map(aov, broom::tidy),
  ) %>% 
  unnest(cols = tidy) %>% 
  drop_na() %>% 
  select(name, sumsq:p.value)

labels <- c(
  'Log<sub>10</sub>(Tortuosity)' = 'Log<sub>10</sub>(Tortuosity)<br>p = 0.001',
  'Log<sub>10</sub>(Distance)' = 'Log<sub>10</sub>(Distance)<br>p = 0.0006',
  'Log<sub>10</sub>(Velocity)' = 'Log<sub>10</sub>(Velocity)<br>p = 0.9',
  'Angular<br>Velocity (SD)' = 'Angular<br>Velocity (SD)<br>p = 0.05',
  'Log<sub>10</sub>(RCD + 1)' = 'Log<sub>10</sub>(RCD + 1)<br>p = 0.04',
  # 'Number of Tracks' = 'Number of Tracks<br>p = 0.054',
  'Total Motility' = 'Total Motility<br>p = 0.4'
)

(track_plot <-  final_df %>% 
    mutate(concentrations = case_when(
      concentrations == '0.1p' ~ 0.0001,
      concentrations == '10uM' ~ 10,
      concentrations == '1uM' ~ 1,
      concentrations == '0.1uM' ~ 0.1,
      concentrations == '0.01uM' ~ 0.01,
      concentrations == '0.001uM' ~ 0.001,
    )) |> 
    ggplot() +
    ggbeeswarm::geom_quasirandom(aes(x = concentrations, y = well_mean),
                                 size = 0.5, alpha = 0.5
    ) +
    stat_summary(aes(x = concentrations, y = well_mean),
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

ex_tracks <- more_filtered %>% 
  filter(batch == '20240307') %>% 
  arrange(well, particle, frame) %>% 
  filter(str_detect(well, '^A'))


(plot_tracks <- ex_tracks %>% 
    mutate(conc = case_when(
      concentrations == '0.1p' ~ "\u2205",
      concentrations == '0.001uM' ~ '0.001 µM',
      concentrations == '0.01uM' ~ '0.01 µM',
      concentrations == '0.1uM' ~ '0.1 µM',
      concentrations == '1uM' ~ '1 µM',
      concentrations == '10uM' ~ '10 µM',
    )) %>% 
    ggplot() +
    ggforce::geom_circle(
      data = tibble(x0 = 1616 / 2, y0 = 1616 / 2, r = 1500 / 2),
      aes(x0 = x0, y0 = y0, r = r), 
      inherit.aes = FALSE) +
    geom_path(
      aes(x = x, y = y, color = frame / 15, group = as.factor(particle))
    ) +
    facet_wrap(facets = vars(conc), ncol = 6, strip.position = "bottom") +
    scale_x_continuous(limits = c(0, 1600)) +
    scale_y_continuous(limits = c(0, 1600)) +
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
                            nrow = 2, rel_heights = c(.3, 1),
                            labels = 'AUTO', label_size = 12)

ggsave(here('pzq', 'Fig4.pdf'), final, width = 5, height = 6, device = cairo_pdf)
ggsave(here('pzq', 'Fig4.png'), final, width = 5, height = 6)
