# ----------------------------
# 1) Paquetes necesarios
# ----------------------------
packages <- c("tidyverse")  # 'tidyverse' incluye dplyr, readr, ggplot2, etc.
to_install <- setdiff(packages, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
lapply(packages, library, character.only = TRUE)


datosA <- read_csv("C:/Users/retrechero/PyCharmMiscProject/gaia_DR3_zoneA.csv") # Zona A (plano galáctico)
datosB <- read_csv("C:/Users/retrechero/PyCharmMiscProject/gaia_DR3_zoneB.csv") # Zona B (polo galáctico)


# 2) comprobación rápido de datos correctos
datosA %>%
  summarise(n = n(),
            b_min = min(b, na.rm=TRUE),
            b_q05 = quantile(b, .05, na.rm=TRUE),
            b_med = median(b, na.rm=TRUE),
            b_q95 = quantile(b, .95, na.rm=TRUE),
            b_max = max(b, na.rm=TRUE),
            frac_cerca_plano = mean(abs(b) <= 5, na.rm=TRUE))  # % cerca del plano

datosB %>%
  summarise(n = n(),
            b_min = min(b, na.rm=TRUE),
            b_q05 = quantile(b, .05, na.rm=TRUE),
            b_med = median(b, na.rm=TRUE),
            b_q95 = quantile(b, .95, na.rm=TRUE),
            b_max = max(b, na.rm=TRUE),
            frac_polo = mean(b >= 80, na.rm=TRUE))             # % en el polo N

A <- datosA
B <- datosB

# verificamos dimensiones y nombres de columnas
dim(A); dim(B)
names(A); names(B)


# ----------------------------
# 4) Función de preparación
# ----------------------------
prep <- function(df, zone_lab){
  # 4.1 Estandariza nombres
  df <- df %>%
    dplyr::rename(
      source_id         = dplyr::any_of("source_id"),
      ra                = dplyr::any_of("ra"),
      dec               = dplyr::any_of("dec"),
      l                 = dplyr::any_of("l"),
      b                 = dplyr::any_of("b"),
      phot_g_mean_mag   = dplyr::any_of("phot_g_mean_mag"),
      bp_rp             = dplyr::any_of("bp_rp"),
      parallax          = dplyr::any_of("parallax"),
      parallax_error    = dplyr::any_of("parallax_error"),
      random_index      = dplyr::any_of("random_index")
    )
  
  # 4.2 Crea razón señal-ruido de la paralaje
  df <- df %>%
    dplyr::mutate(
      parallax_over_error = parallax / parallax_error,
      zone = zone_lab
    )
  
  # 4.3 Filtro suave de color BP-RP:
  df <- df %>%
    dplyr::filter(is.na(bp_rp) | (bp_rp > -2 & bp_rp < 6))
  
  # 4.4 Elimina 'random_index' (ya no es útil tras la descarga)
  df <- df %>% dplyr::select(-dplyr::any_of("random_index"))
  
  return(df)
}

# ----------------------------
# 5) Preparar A y B y unir
# ----------------------------
A_clean <- prep(A, "A")
B_clean <- prep(B, "B")

# Verificaciones rápidas
dim(A_clean); dim(B_clean)
summary(A_clean$phot_g_mean_mag)
summary(B_clean$phot_g_mean_mag)
summary(A_clean$bp_rp)
summary(B_clean$bp_rp)

# Data frame combinado para análisis comparativo
AB <- dplyr::bind_rows(A_clean, B_clean)

# ----------------------------
# 6) Checks mínimos de calidad
# ----------------------------

# 6.1 Variables clave presentes
vars_required <- c("phot_g_mean_mag", "bp_rp", "parallax", "parallax_error", "zone")
setdiff(vars_required, names(AB))  # debería devolver character(0)

# 6.2 Rangos razonables
summary(AB$phot_g_mean_mag)
summary(AB$bp_rp)

has_lb  <- all(c("l","b") %in% names(AB))
has_radec <- all(c("ra","dec") %in% names(AB))
cat("Tiene l/b:", has_lb, " | Tiene ra/dec:", has_radec, "\n")


if (has_lb) {
  AB %>%
    dplyr::group_by(zone) %>%
    dplyr::summarise(
      n = dplyr::n(),
      b_min = min(b, na.rm=TRUE),
      b_q05 = quantile(b, .05, na.rm=TRUE),
      b_med = median(b, na.rm=TRUE),
      b_q95 = quantile(b, .95, na.rm=TRUE),
      b_max = max(b, na.rm=TRUE),
      frac_plano = mean(abs(b) <= 5, na.rm=TRUE),  # debería ser alto en Zona A
      frac_polo  = mean(b >= 80, na.rm=TRUE)       # debería ser alto en Zona B
    ) %>% print(n=Inf)
}

#########################
# SECCIÓN 1: Análisis exploratorio de datos (EDA)
# ------------------------------------------------
# Objetivo:
# - Comparar las distribuciones de magnitud (phot_g_mean_mag)
#   y color (bp_rp) entre las zonas A (plano galáctico) y B (polo galáctico).
# - Resumir medias, dispersiones y tamaños muestrales por zona
#   contexto para las pruebas de hipótesis posteriores.
#########################

# 2.1 Resúmenes numéricos por zona
summary_by_zone <- AB %>%
  group_by(zone) %>%
  summarise(
    n = n(),
    g_mean = mean(phot_g_mean_mag, na.rm = TRUE),
    g_sd   = sd(phot_g_mean_mag,   na.rm = TRUE),
    g_q25  = quantile(phot_g_mean_mag, 0.25, na.rm = TRUE),
    g_med  = median(phot_g_mean_mag,   na.rm = TRUE),
    g_q75  = quantile(phot_g_mean_mag, 0.75, na.rm = TRUE),
    color_mean = mean(bp_rp, na.rm = TRUE),
    color_sd   = sd(bp_rp,   na.rm = TRUE),
    color_q25  = quantile(bp_rp, 0.25, na.rm = TRUE),
    color_med  = median(bp_rp,   na.rm = TRUE),
    color_q75  = quantile(bp_rp, 0.75, na.rm = TRUE),
    parallaxError_mean = mean(parallax_over_error, na.rm = TRUE),
    parallaxError_sd = sd(parallax_over_error, na.rm = TRUE),
    parallaxError_med = median(parallax_over_error, na.rm = TRUE),
    parallaxError_q25 = quantile(parallax_over_error, 0,25, na.rm = TRUE),
    parallaxError_q275 = quantile(parallax_over_error, 0,75, na.rm = TRUE)
    
  )
summary_by_zone

view(AB)

# ----------------------------
# Histograma de magnitud G por zona
# ----------------------------
ggplot(AB, aes(phot_g_mean_mag, fill = zone)) +
  geom_histogram(alpha = 0.45, bins = 80, position = "identity") +
  labs(title = "Distribución de magnitud G por zona",
       x = "phot_g_mean_mag (G) — menores = más brillantes",
       y = "conteo") +
  theme_minimal()

# 2.4. Boxplots comparativos (mediana y dispersión)
# visualizamos medianas (línea central) y dispersión (caja y “bigotes”) por zona.
# Esto ayuda a anticipar las pruebas de varianza.
# Diferencia de medianas entre zonas.

# Diferencias de anchura (IQR): si A es más ancha en BP−RP, 
# su varianza podría ser mayor
# ----------------------------
# Boxplot de magnitud G por zona
# ----------------------------
ggplot(AB, aes(zone, phot_g_mean_mag, fill = zone)) +
  geom_boxplot(outlier.alpha = 0.1) +
  labs(title = "Magnitud G por zona (boxplot)",
       x = "Zona", y = "phot_g_mean_mag") +
  theme_minimal()

# ----------------------------
# Boxplot de color BP−RP por zona
# ----------------------------
ggplot(AB, aes(zone, bp_rp, fill = zone)) +
  geom_boxplot(outlier.alpha = 0.1) +
  labs(title = "BP−RP por zona (boxplot)",
       x = "Zona", y = "bp_rp") +
  theme_minimal()


# Distribución de las Variables
plot_var_dist <- function(df, var){
  ggplot(df %>% filter(!is.na(.data[[var]])),
         aes(x = .data[[var]], fill = zone)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 80, alpha = 0.45, position = "identity") +
    geom_density(alpha = 0.8) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste("Distribución original de", var, "por zona"),
      x = var, y = "Densidad", fill = "Zona"
    ) +
    theme_minimal(base_size = 11)
}

plot_var_dist(AB, "phot_g_mean_mag")
plot_var_dist(AB, "bp_rp")

p99 <- quantile(AB$parallax_over_error, 0.99, na.rm = TRUE)

ggplot(AB %>% filter(parallax_over_error <= p99),
       aes(x = parallax_over_error, fill = zone)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, alpha = 0.45, position = "identity") +
  geom_density(alpha = 0.7) +
  labs(title="Distribución original limitada (p99) de parallax_over_error",
       x="parallax_over_error", y="Densidad") +
  theme_minimal()

write.csv(AB, "AB.csv", row.names = FALSE)
getwd()


############################################
# SECCIÓN 2: Inferencia estadística
# ------------------------------------------
# Demostraciones numéricas de:
# - Teorema del límite central (TLC)
# - Consistencia de la media muestral
# - Suficiencia bajo normalidad
# - Intervalos de confianza para la media
# - Teorema de Rao–Blackwell
###########################################

# 1. Ilustración del Teorema del límite Central
# ==============================
# TLC: panel con Z
# ==============================
comparar_clt <- function(df, var, zone_sel = c("A","B"),
                         ns   = c(5, 10, 30, 100, 500, 1000),
                         reps = 5000,
                         trim_q = c(0.01, 0.99),   # recorte para legibilidad
                         seed = 2025) {
  
  zone_sel <- match.arg(zone_sel)
  stopifnot(var %in% names(df))
  
  set.seed(seed)
  
  # 1) Extrae y limpia la variable en la zona elegida
  x <- df %>% filter(zone == zone_sel) %>% pull(all_of(var))
  x <- x[is.finite(x)]
  if (!is.null(trim_q)) {
    qq <- quantile(x, trim_q, na.rm = TRUE)
    x  <- x[x >= qq[1] & x <= qq[2]]
  }
  
  mu_hat  <- mean(x)
  sig_hat <- sd(x)
  
  # 2) Simula medias y estandariza para cada n
  sim_list <- lapply(ns, function(n){
    xbar <- replicate(reps, mean(sample(x, size = n, replace = TRUE)))
    Z    <- (xbar - mu_hat) / (sig_hat / sqrt(n))
    tibble(
      panel = factor(paste0("Z (n=", n, ")"), levels = paste0("Z (n=", ns, ")")),
      value = Z
    )
  })
  
  df_all  <- bind_rows(sim_list)
  df_norm <- tibble(panel = factor(paste0("Z (n=", ns, ")"),
                                   levels = levels(df_all$panel)))
  
  # 3) Gráfico final
  ggplot(df_all, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "gray75") +
    facet_grid(panel ~ ., scales = "free_x") +
    stat_function(
      data = df_norm,
      fun = dnorm, args = list(mean = 0, sd = 1),
      linewidth = 0.9, inherit.aes = FALSE, aes(x = ..x..),
      xlim = c(-4, 4)
    ) +
    labs(
      title = paste0("TLC — ", var, " (zona ", zone_sel, ")"),
      subtitle = "Distribución de Zₙ = (X̄ − μ) / (σ/√n) para distintos tamaños muestrales",
      x = "Z estandarizado", y = "densidad"
    ) +
    theme_minimal(base_size = 11)
}

# Zona A
comparar_clt(AB, var = "parallax_over_error", zone_sel = "A")
comparar_clt(AB, var = "phot_g_mean_mag", zone_sel = "A")
# Zona B
comparar_clt(AB, var = "parallax_over_error", zone_sel = "B")
comparar_clt(AB, var = "phot_g_mean_mag", zone_sel = "B")


#####
# Consistencia
####
consistencia <- function(df, var = "bp_rp", zone = "A",
                         ns = c(20, 50, 100, 200, 500, 1000),
                         reps = 1000, trim_q = NULL, seed = 123) {
  set.seed(seed)
  stopifnot(var %in% names(df))
  
  # 1) Datos base
  x <- df %>% filter(zone == zone) %>% pull(all_of(var))
  x <- x[is.finite(x)]
  if (!is.null(trim_q)) {
    q <- quantile(x, trim_q, na.rm = TRUE); x <- x[x >= q[1] & x <= q[2]]
  }
  mu_hat <- mean(x)
  
  # 2) Simulación de medias para cada n
  sim <- lapply(ns, function(n) {
    tibble(n = n,
           xbar = replicate(reps, mean(sample(x, n, replace = TRUE))))
  }) %>% bind_rows() %>%
    mutate(n = factor(n, levels = ns))
  
  # 3) Gráfico: concentración de medias alrededor de mu_hat
  ggplot(sim, aes(n, xbar)) +
    geom_violin(fill = "gray85", color = "gray40", scale = "width") +
    geom_boxplot(width = 0.15, outlier.alpha = 0.15) +
    geom_hline(yintercept = mu_hat, linetype = 2, linewidth = 0.8) +
    labs(
      title = paste0("Consistencia de la media — ", var, " (zona ", zone, ")"),
      subtitle = expression(paste("Distribución de ", bar(X)[n], " por tamaño muestral; línea punteada = ", hat(mu))),
      x = "tamaño muestral n", y = expression(bar(X)[n])
    ) +
    theme_minimal(base_size = 11)
}

consistencia(AB, var = "bp_rp", zone = "A")

#####
# Suficiencia
####
suficiencia <- function(df, var = "bp_rp", zone = "A", n = 2000, seed = 123) {
  set.seed(seed)
  
  # 1) Tomamos muestra real
  x <- df %>% filter(zone == zone) %>% pull(all_of(var))
  x <- x[is.finite(x)]
  sample_x <- sample(x, n)
  
  # 2) Estadísticos
  xbar <- mean(sample_x)
  s    <- sd(sample_x)
  
  # 3) Curva normal reconstruida
  xx <- seq(min(sample_x), max(sample_x), length.out = 300)
  dens_model <- dnorm(xx, mean = xbar, sd = s)
  df_curve <- data.frame(xx, dens_model)
  
  # 4) Gráfico
  ggplot(data.frame(x = sample_x), aes(x)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "gray80", color = "black") +
    geom_line(data = df_curve,
              aes(x = xx, y = dens_model),
              inherit.aes = FALSE,
              linewidth = 1.2, color = "red") +
    labs(
      title = paste0("Suficiencia (ilustración) — ", var, " (zona ", zone, ")"),
      subtitle = paste0("Histograma real vs Normal(", round(xbar,3), ", ", round(s,3), ")"),
      x = "valor observado", y = "densidad"
    ) +
    theme_minimal(base_size = 11)
}
# Ejemplo:
suficiencia(AB, var = "bp_rp", zone = "A")

#####
# Intervalos de Confianza
####
ci_demo <- function(df, var = "bp_rp", zone = "A",
                    n_vec = c(30, 200), conf = 0.95, reps = 400,
                    seed = 123) {
  set.seed(seed)
  stopifnot(var %in% names(df))
  
  # Datos y "verdadero" (empírico) mu_hat
  x <- df %>% filter(zone == zone) %>% pull(all_of(var))
  x <- x[is.finite(x)]
  mu_hat <- mean(x)
  
  alpha <- 1 - conf
  
  sim <- lapply(n_vec, function(n){
    replicate(reps, {
      xi <- sample(x, n, replace = TRUE)
      xb <- mean(xi); s <- sd(xi)
      tcrit <- qt(1 - alpha/2, df = n - 1)
      hw <- tcrit * s / sqrt(n)
      c(n = n, center = xb, lo = xb - hw, hi = xb + hw)
    }) |> t() |> as.data.frame()
  }) |> bind_rows()
  
  sim <- sim |>
    mutate(n = factor(n, levels = n_vec),
           cover = (lo <= mu_hat & hi >= mu_hat)) |>
    arrange(n, center) |>
    group_by(n) |>
    mutate(idx = row_number()) |>
    ungroup()
  
  # Gráfico: 400 IC por panel, línea vertical = mu_hat
  p <- ggplot(sim, aes(y = idx)) +
    geom_segment(aes(x = lo, xend = hi, yend = idx, color = cover), linewidth = 0.8) +
    geom_point(aes(x = center), size = 0.8) +
    geom_vline(xintercept = mu_hat, linetype = 2) +
    scale_color_manual(values = c("TRUE" = "seagreen4", "FALSE" = "firebrick")) +
    facet_grid(n ~ ., scales = "free_x") +
    labs(
      title = paste0("IC para la media de ", var, " — zona ", zone),
      subtitle = paste0(conf*100, "% IC t de Student; línea vertical = μ̂; verde=cubre, rojo=no"),
      x = "valor", y = "réplica"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  
  # Resumen de cobertura y ancho
  resumen <- sim |>
    summarise(
      cobertura = mean(cover),
      ancho_medio = mean(hi - lo),
      mediana_ancho = median(hi - lo)
      , .by = n)
  
  list(plot = p, resumen = resumen)
}

out <- ci_demo(AB, var = "bp_rp", zone = "A", n_vec = c(30, 200), conf = 0.95, reps = 400)
out$plot
out$resumen

# ahora ajustamos a 0.99 para ver como cambia el ancho
out2 <- ci_demo(AB, var = "bp_rp", zone = "B", n_vec = c(30, 200), conf = 0.99, reps = 400)
out2$plot
out2$resumen


#######################
# RAO - BLACKWELL
#######################
rao_blackwell_demo <- function(df, var = "bp_rp", zone = "A",
                               ns = c(10, 50, 200), reps = 2000, seed = 123){
  set.seed(seed)
  stopifnot(var %in% names(df))
  
  # "verdadero" (empírico) con todos los datos disponibles
  x <- df %>% filter(zone == zone) %>% pull(all_of(var))
  x <- x[is.finite(x)]
  mu_hat <- mean(x)
  
  # Simulación: para cada n, reps veces tomamos una muestra con reemplazo
  sim <- lapply(ns, function(n){
    replicate(reps, {
      xi <- sample(x, size = n, replace = TRUE)
      c(n = n,
        naive = xi[1],           # estimador ingenuo: 1 dato aleatorio
        rb    = mean(xi))        # Rao–Blackwell: media muestral
    }) |> t() |> as.data.frame()
  }) |> bind_rows() |>
    mutate(n = factor(n, levels = ns))
  
  # Resumen: sesgo, varianza y eficiencia relativa
  resumen <- sim |>
    summarise(
      bias_naive = mean(naive) - mu_hat,
      var_naive  = var(naive),
      bias_rb    = mean(rb)    - mu_hat,
      var_rb     = var(rb),
      eff_rel    = var_naive / var_rb      # > 1 implica mejora grande
      , .by = n)
  
  # Datos en formato largo para graficar
  long <- sim |>
    tidyr::pivot_longer(cols = c(naive, rb),
                        names_to = "estimador", values_to = "valor") |>
    mutate(estimador = factor(estimador, levels = c("naive","rb"),
                              labels = c("Ingenuo: X1", "Rao–Blackwell: media")))
  
  # Gráfico: distribución de los estimadores (violines + boxplots)
  p <- ggplot(long, aes(x = estimador, y = valor, fill = estimador)) +
    geom_violin(color = "gray35", scale = "width") +
    geom_boxplot(width = 0.15, outlier.alpha = 0.15) +
    geom_hline(yintercept = mu_hat, linetype = 2, linewidth = 0.9) +
    facet_wrap(~ n, nrow = 1) +
    scale_fill_manual(values = c("gray75","gray60")) +
    labs(
      title = paste0("Rao–Blackwell — ", var, " (zona ", zone, ")"),
      subtitle = expression(paste("Línea punteada = ", hat(mu),
                                  ".  Comparación de varianza entre estimador ingenuo y RB.")),
      x = NULL, y = "estimación de μ"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  
  list(plot = p, resumen = resumen)
}

out <- rao_blackwell_demo(AB, var = "bp_rp", zone = "A", ns = c(10,50,200), reps = 2000)
out$plot
out$resumen



######################
# SECCIÓN 3: Pruebas de hipótesis
# --------------------------------
# Comparación formal entre zonas A y B para:
#  - Magnitud G (phot_g_mean_mag)
#  - Color BP−RP (bp_rp)
#  - Calidad de paralaje (parallax_over_error)
# Incluye dos pruebas educativas para ilustrar errores Tipo I y Tipo II.
######################
set.seed(123)
library(dplyr); library(tidyr)

# Utilidades --------------------------------------------------------------
take_sample <- function(df, var, n_per_zone){
  df %>%
    filter(zone %in% c("A","B")) %>%
    select(zone, {{var}}) %>%
    drop_na() %>%
    group_by(zone) %>%
    slice_sample(n = n_per_zone, replace = FALSE) %>%
    ungroup() %>%
    mutate(zone = factor(zone, levels = c("A","B")))
}

brief_desc <- function(df, var){
  df %>% group_by(zone) %>% summarise(
    n = n(),
    mean = mean({{var}}),
    sd = sd({{var}})
  )
}

# Prueba 1 — Magnitud G: comparación de medias entre zonas (t de Welch)
# H0: μ_A = μ_B
# H1: μ_A ≠ μ_B
# Con n = 500 por zona se espera no rechazar H0, ya que la diferencia real de medias es muy pequeña.
G_mean <- take_sample(AB, phot_g_mean_mag, n_per_zone = 500)
desc_1 <- brief_desc(G_mean, phot_g_mean_mag); print(desc_1)

tt_g <- t.test(phot_g_mean_mag ~ zone, data = G_mean, var.equal = FALSE, conf.level = 0.95)
res_1 <- tibble::tibble(
  prueba   = "g: t de Welch (medias)",
  n_por_z  = 500,
  diff_A_B = with(desc_1, mean[zone=="A"] - mean[zone=="B"]),
  t_stat   = unname(tt_g$statistic),
  gl       = unname(tt_g$parameter),
  p_value  = tt_g$p.value,
  ci_lo    = tt_g$conf.int[1],
  ci_hi    = tt_g$conf.int[2]
); print(res_1)

df_ic <- data.frame(
  estimate = res_1$diff_A_B,
  lo = res_1$ci_lo,
  hi = res_1$ci_hi
)

ggplot(df_ic, aes(x = "Zonas A-B", y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .1, size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color="red") +
  labs(
    x = "",
    y = "Diferencia de medias (A − B)",
    title = "Intervalo de confianza del 95% para la diferencia de medias"
  ) +
  theme_minimal(base_size = 14)


# Prueba 2 — Magnitud G: comparación de varianzas entre zonas (F de Fisher)
# H0: σ_A^2 = σ_B^2
# H1: σ_A^2 ≠ σ_B^2
# Con n = 100 por zona la diferencia de varianzas es suficientemente grande como para esperar rechazo de H0.
G_var <- take_sample(AB, phot_g_mean_mag, n_per_zone = 100)
desc_2 <- brief_desc(G_var, phot_g_mean_mag); print(desc_2)

vf_g <- var.test(phot_g_mean_mag ~ zone, data = G_var, conf.level = 0.95)
res_2 <- tibble::tibble(
  prueba   = "g: F de Fisher (varianzas)",
  n_por_z  = 100,
  ratioVar = unname(vf_g$estimate),  # Var(A)/Var(B)
  F_stat   = unname(vf_g$statistic),
  gl1      = vf_g$parameter[1], gl2 = vf_g$parameter[2],
  p_value  = vf_g$p.value,
  ci_lo    = vf_g$conf.int[1], ci_hi = vf_g$conf.int[2]
); print(res_2)

ic  <- vf_g$conf.int
est <- unname(vf_g$estimate)  # razón Var(A)/Var(B)

plot(1, est, pch = 19,
     ylim = c(ic[1], ic[2]),
     xlab = "", ylab = "Var(A) / Var(B)",
     main = "IC 95% para la razón de varianzas")
segments(1, ic[1], 1, ic[2], lwd = 3)
abline(h = 1, lty = 2, col = "red")  # línea en 1 (igualdad de varianzas)

# Prueba 3 — Color BP−RP: comparación de medias entre zonas (t de Welch)
# H0: μ_A = μ_B
# H1: μ_A ≠ μ_B
# La diferencia real entre zonas es grande; incluso con n = 30 por zona se espera rechazo de H0.
C_mean <- take_sample(AB, bp_rp, n_per_zone = 30)
desc_3 <- brief_desc(C_mean, bp_rp); print(desc_3)

tt_c <- t.test(bp_rp ~ zone, data = C_mean, var.equal = FALSE, conf.level = 0.95)
res_3 <- tibble::tibble(
  prueba   = "Color: t de Welch (medias)",
  n_por_z  = 30,
  diff_A_B = with(desc_3, mean[zone=="A"] - mean[zone=="B"]),
  t_stat   = unname(tt_c$statistic),
  gl       = unname(tt_c$parameter),
  p_value  = tt_c$p.value,
  ci_lo    = tt_c$conf.int[1],
  ci_hi    = tt_c$conf.int[2]
); print(res_3)

df_ic <- data.frame(
  estimate = est,
  lo = ic[1],
  hi = ic[2]
)

ggplot(df_ic, aes(x = "Color A-B", y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, size = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "",
    y = "Diferencia de medias (A − B)",
    title = "Intervalo de confianza del 95% — Color BP–RP"
  ) +
  theme_minimal(base_size = 14)


# Prueba 4 — Color BP−RP: comparación de varianzas entre zonas (F de Fisher)
# H0: σ_A^2 = σ_B^2
# H1: σ_A^2 ≠ σ_B^2
# Las varianzas son muy similares. con n = 200 por zona se anticipa que no se rechace H0.
C_var <- take_sample(AB, bp_rp, n_per_zone = 200)
desc_4 <- brief_desc(C_var, bp_rp); print(desc_4)

vf_c <- var.test(bp_rp ~ zone, data = C_var, conf.level = 0.95)
res_4 <- tibble::tibble(
  prueba   = "Color: F de Fisher (varianzas)",
  n_por_z  = 200,
  ratioVar = unname(vf_c$estimate),
  F_stat   = unname(vf_c$statistic),
  gl1      = vf_c$parameter[1], gl2 = vf_c$parameter[2],
  p_value  = vf_c$p.value,
  ci_lo    = vf_c$conf.int[1], ci_hi = vf_c$conf.int[2]
); print(res_4)


df_ic <- data.frame(
  estimate = est,
  lo = ic[1],
  hi = ic[2]
)

ggplot(df_ic, aes(x = "Zonas A-B", y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, size = 1.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = "",
    y = "Razón Var(A) / Var(B)",
    title = "Intervalo de confianza del 95% — Varianzas BP–RP"
  ) +
  theme_minimal(base_size = 14)


# Prueba 5 — parallax_over_error: comparación de medias entre zonas (t de Welch)
# H0: μ_A = μ_B
# H1: μ_A ≠ μ_B
# La diferencia entre zonas es muy grande, aunque con alta dispersión; con n = 200 por zona debería observarse rechazo de H0.
P_mean <- take_sample(AB, parallax_over_error, n_per_zone = 200)
desc_5 <- brief_desc(P_mean, parallax_over_error); print(desc_5)

tt_p <- t.test(parallax_over_error ~ zone, data = P_mean, var.equal = FALSE, conf.level = 0.95)
res_5 <- tibble::tibble(
  prueba   = "Parallax_over_error: t de Welch (medias)",
  n_por_z  = 200,
  diff_A_B = with(desc_5, mean[zone=="A"] - mean[zone=="B"]),
  t_stat   = unname(tt_p$statistic),
  gl       = unname(tt_p$parameter),
  p_value  = tt_p$p.value,
  ci_lo    = tt_p$conf.int[1],
  ci_hi    = tt_p$conf.int[2]
); print(res_5)

df_ic <- data.frame(
  estimate = est,
  lo = ic[1],
  hi = ic[2]
)

ggplot(df_ic, aes(x = "A-B", y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, size = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "",
    y = "Diferencia de medias (A − B)",
    title = "IC 95% — parallax_over_error"
  ) +
  theme_minimal(base_size = 14)


# Prueba 6 — Ejemplo educativo de error Tipo I (falso positivo)
# --------------------------------------------------------------
# Variable: brillo G (phot_g_mean_mag)
# Diseño deliberadamente deficiente:
#  - Muestra muy pequeña (n = 5 por zona)
# Objetivo: mostrar que una mala especificación del modelo puede llevar a rechazar H0 sin justificación.
set.seed(4)  

# 1) Submuestra pequeña por zona
G_small <- AB %>%
  filter(zone %in% c("A","B")) %>%
  select(zone, phot_g_mean_mag) %>%
  drop_na(phot_g_mean_mag) %>%
  group_by(zone) %>%
  slice_sample(n = 5, replace = FALSE) %>%
  ungroup() %>%
  mutate(zone = factor(zone, levels = c("A","B")))

# 2) Descriptivos muestrales
desc_small <- G_small %>%
  group_by(zone) %>%
  summarise(n = n(),
            mean_g = mean(phot_g_mean_mag),
            sd_g   = sd(phot_g_mean_mag),
            .groups = "drop")
print(desc_small)

# 3) PRUEBA MAL PLANTEADA: t de Student
tt_bad <- t.test(phot_g_mean_mag ~ zone, data = G_small,
                 var.equal = TRUE,  # <- suposición errónea
                 alternative = "two.sided", conf.level = 0.95)

res_bad <- tibble(
  prueba   = "g: t de Student (var.equal=TRUE) — MAL PLANTEADA",
  n_por_z  = 5,
  diff_A_B = with(desc_small, mean_g[zone=="A"] - mean_g[zone=="B"]),
  t_stat   = unname(tt_bad$statistic),
  gl       = unname(tt_bad$parameter),
  p_value  = tt_bad$p.value,
  ci_lo    = tt_bad$conf.int[1],
  ci_hi    = tt_bad$conf.int[2]
)
print(res_bad)


############
# Prueba 7 — Ejemplo educativo de error Tipo II (falso negativo)
############
# Variable: color (bp_rp)
# Condiciones: n = 10 por zona, α = 0.01 (dos colas)
# La prueba (t de Welch) está bien planteada, pero la muestra es pequeña y el nivel de significancia muy estricto,
# por lo que la potencia es baja y es fácil no detectar una diferencia real.
set.seed(2025)
alpha <- 0.01

# 1) Submuestra pequeña para color
C_small <- AB %>%
  filter(zone %in% c("A","B")) %>%
  select(zone, bp_rp) %>%
  drop_na(bp_rp) %>%
  group_by(zone) %>%
  slice_sample(n = 10, replace = FALSE) %>%
  ungroup() %>%
  mutate(zone = factor(zone, levels = c("A","B")))

# 2) Descriptivos muestrales
desc <- C_small %>%
  group_by(zone) %>%
  summarise(n = n(),
            mean = mean(bp_rp),
            sd   = sd(bp_rp),
            .groups = "drop")
print(desc)

# 3) Prueba t de Welch con α = 0.01
tt <- t.test(bp_rp ~ zone, data = C_small,
             var.equal = FALSE,
             alternative = "two.sided",
             conf.level = 1 - alpha)

res_edu <- tibble(
  prueba   = "Color: t de Welch (educativa, n=10, α=0.01)",
  n_por_z  = 10,
  diff_A_B = with(desc, mean[zone=="A"] - mean[zone=="B"]),
  t_stat   = unname(tt$statistic),
  gl       = unname(tt$parameter),
  p_value  = tt$p.value,
  ci99_lo  = tt$conf.int[1],
  ci99_hi  = tt$conf.int[2],
  decision = ifelse(tt$p.value < alpha, "Rechaza H0 (α=0.01)", "NO rechaza H0 (α=0.01)")
)
print(res_edu)

# 4) Cálculo aproximado de n requerido para 80% de potencia
# Parámetros aproximados basados en los datos originales
alpha    = alpha
power    = 0.80
Delta    = 0.52   # diferencia empírica aproximada de medias entre zonas
sigma    = 0.65   # desviación estándar aproximada de bp_rp
z_alpha  = qnorm(1 - alpha/2)
z_beta   = qnorm(power)
n_per_group <- ((z_alpha + z_beta) * sigma / Delta)^2

n_per_group

# 5) Verificación con un tamaño muestral cercano
n <- 40
C_n <- AB %>%
  filter(zone %in% c("A","B")) %>%
  select(zone, bp_rp) %>%
  drop_na(bp_rp) %>%
  group_by(zone) %>%
  slice_sample(n = n, replace = FALSE) %>%
  ungroup() %>%
  mutate(zone = factor(zone, levels = c("A","B")))

desc <- C_n %>%
  group_by(zone) %>%
  summarise(n = n(),
            mean = mean(bp_rp),
            sd   = sd(bp_rp),
            .groups = "drop")
print(desc)

tt <- t.test(bp_rp ~ zone, data = C_n,
             var.equal = FALSE,
             alternative = "two.sided",
             conf.level = 0.99)

tibble(
  prueba   = "Color: t de Welch (n≈40, α=0.01)",
  n_por_z  = n,
  diff_A_B = with(desc, mean[zone=="A"] - mean[zone=="B"]),
  t_stat   = unname(tt$statistic),
  gl       = unname(tt$parameter),
  p_value  = tt$p.value,
  ci99_lo  = tt$conf.int[1],
  ci99_hi  = tt$conf.int[2],
  decision = ifelse(tt$p.value < 0.01, "Rechaza H0 (α=0.01)", "NO rechaza H0 (α=0.01)")
)


