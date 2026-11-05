
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(grid)

# ------------------------------------------
# 1. Read panels A-C
# ------------------------------------------

dat_abc <- all_panels_combined

dat_abc <- dat_abc %>%
  select(Panel, Combined_Sample_Size, Method, Power)

# ------------------------------------------
# 2. Read and reshape panel D
# ------------------------------------------

dat_d <- panel_D %>%
  dplyr::select(-any_of("...1")) %>%
  pivot_longer(
    cols = -SampleSize,
    names_to = "Method",
    values_to = "Power"
  ) %>%
  mutate(
    Method = recode(
      Method,
      glmpwr = "Logistic Regression, 7-Day PPA",
      any1glmpwr = "Logistic Regression, Any Cigarette Free Days",
      nbpwr = "Negative Binomial, Number of Cigarette Free Days",
      AGpwr = "Andersen-Gill Model, Recurrent Event",
      cox = "Stratified Cox Regression, Time-to-Event"
    ),
    Panel = "d",
    Combined_Sample_Size = SampleSize
  ) %>%
  dplyr::select(Panel, Combined_Sample_Size, Method, Power)

# ------------------------------------------
# 3. Combine all four panels
# ------------------------------------------

dat_all <- bind_rows(dat_abc, dat_d)

# Set panel and method order
method_order <- c(
  "Logistic Regression, 7-Day PPA",
  "Andersen-Gill Model, Recurrent Event",
  "Logistic Regression, Any Cigarette Free Days",
  "Stratified Cox Regression, Time-to-Event",
  "Negative Binomial, Number of Cigarette Free Days"
)

dat_all <- dat_all %>%
  mutate(
    Panel = factor(
      Panel,
      levels = c("a", "b", "c", "d")
    ),
    Method = factor(
      Method,
      levels = method_order
    )
  )

# ------------------------------------------
# 4. Define black-and-white aesthetics
# ------------------------------------------

line_types <- c(
  "solid",
  "dashed",
  "dotted",
  "dotdash",
  "longdash"
)

point_shapes <- c(
  16,  # Circle
  17,  # Triangle
  15,  # Square
  18,  # Diamond
  3    # Plus
)

# ------------------------------------------
# 5. Generate combined figure
# ------------------------------------------

p <- ggplot(
  dat_all,
  aes(
    x = Combined_Sample_Size,
    y = Power,
    group = Method,
    linetype = Method,
    shape = Method
  )
) +
  
  geom_line(
    color = "black",
    linewidth = 0.8
  ) +
  
  geom_point(
    color = "black",
    size = 2.3,
    stroke = 0.8
  ) +
  
  scale_linetype_manual(
    values = setNames(line_types, method_order)
  ) +
  
  scale_shape_manual(
    values = setNames(point_shapes, method_order)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  facet_wrap(
    ~Panel,
    ncol = 2,
    scales = "free_x",
    labeller = labeller(
      Panel = c(
        a = "(a)",
        b = "(b)",
        c = "(c)",
        d = "(d)"
      )
    )
  ) +
  
  labs(
    x = "Combined Sample Size",
    y = "Power"
  ) +
  
  # Two-row legend
  guides(
    linetype = guide_legend(
      nrow = 2,
      byrow = TRUE
    ),
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    panel.grid.minor = element_blank(),
    
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    ),
    
    axis.title = element_text(size = 12),
    
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "white",
      color = "black"
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.position = "bottom",
    
    legend.title = element_blank(),
    
    legend.text = element_text(size = 9),
    
    legend.key.width = unit(1.1, "cm"),
    
    legend.box = "vertical",
    
    legend.spacing.y = unit(0.1, "cm"),
    
    plot.margin = margin(10, 10, 10, 10)
  )

# Display figure
print(p)

# ------------------------------------------
# 6. Save publication-quality figures
# ------------------------------------------

ggsave(
  "all_four_panels_black_white.png",
  p,
  width = 11,
  height = 9,
  dpi = 300
)

ggsave(
  "all_four_panels_black_white.pdf",
  p,
  width = 11,
  height = 9
)
