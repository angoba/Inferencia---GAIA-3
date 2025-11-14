
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.animation import FFMpegWriter

AB = pd.read_csv("AB.csv")

def simular_Zn(data, variable, zone, ns, reps=5000, trim_quantiles=(0.01, 0.99), random_state=2025):
    rng = np.random.default_rng(random_state)
    x = data.loc[data["zone"] == zone, variable].astype(float)
    x = x.replace([np.inf, -np.inf], np.nan).dropna().values

    if trim_quantiles is not None:
        ql, qh = np.quantile(x, trim_quantiles)
        x = x[(x >= ql) & (x <= qh)]

    mu_hat = x.mean()
    sigma_hat = x.std(ddof=1)

    results = {}
    for n in ns:
        samples = rng.choice(x, size=(reps, n), replace=True)
        xbar = samples.mean(axis=1)
        Z = (xbar - mu_hat) / (sigma_hat / np.sqrt(n))
        results[n] = Z

    return results

# -------------------------
# Configuración TLC
# -------------------------
ns = [5, 10, 20, 30, 50, 80, 120, 200, 300, 500, 800, 1000]
reps = 5000

Z_g_A   = simular_Zn(AB, "phot_g_mean_mag",     "A", ns, reps=reps, random_state=2025)
Z_snr_A = simular_Zn(AB, "parallax_over_error", "A", ns, reps=reps, random_state=2026)
Z_snr_B = simular_Zn(AB, "parallax_over_error", "B", ns, reps=reps, random_state=2027)

# -------------------------
# Preparar figura
# -------------------------
fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
x_grid = np.linspace(-4, 4, 400)
y_norm = norm.pdf(x_grid, 0, 1)

writer = FFMpegWriter(fps=1)

with writer.saving(fig, "TLC_evolucion_video.mp4", dpi=120):
    for n in ns:
        axes[0].cla()
        axes[1].cla()
        axes[2].cla()

        # Panel 1: G Zona A
        Z1 = Z_g_A[n]
        axes[0].hist(Z1, bins=50, density=True, alpha=0.7)
        axes[0].plot(x_grid, y_norm, lw=2)
        axes[0].set_title(f"G (Zona A)\nn={n}")
        axes[0].set_xlim(-4, 4)

        # Panel 2: Parallax Zona A
        Z2 = Z_snr_A[n]
        axes[1].hist(Z2, bins=50, density=True, alpha=0.7)
        axes[1].plot(x_grid, y_norm, lw=2)
        axes[1].set_title(f"Parallax (Zona A)\nn={n}")
        axes[1].set_xlim(-4, 4)

        # Panel 3: Parallax Zona B
        Z3 = Z_snr_B[n]
        axes[2].hist(Z3, bins=50, density=True, alpha=0.7)
        axes[2].plot(x_grid, y_norm, lw=2)
        axes[2].set_title(f"Parallax (Zona B)\nn={n}")
        axes[2].set_xlim(-4, 4)

        fig.suptitle("Teorema del Límite Central — Evolución simultánea", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.92])

        writer.grab_frame()

print("Video generado: TLC_evolucion_video.mp4")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

AB = pd.read_csv("AB.csv")

# Parámetros de la simulación
variable = "bp_rp"
zone     = "A"
ns       = [20, 50, 100, 200, 500, 1000]
reps     = 1000
trim_quantiles = None

# -------------------------
# 2. Preparar datos base
# -------------------------
x = AB.loc[AB["zone"] == zone, variable].astype(float)
x = x.replace([np.inf, -np.inf], np.nan).dropna().values

if trim_quantiles is not None:
    ql, qh = np.quantile(x, trim_quantiles)
    x = x[(x >= ql) & (x <= qh)]

mu_hat = x.mean()
sigma_hat = x.std(ddof=1)

print(f"Variable: {variable}, zona: {zone}")
print(f"mu_hat = {mu_hat:.4f}, sigma_hat = {sigma_hat:.4f}, N base = {len(x)}")

# -------------------------
# 3. Simulación de medias
# -------------------------
means_by_n = {}
rmse_emp   = {}

rng = np.random.default_rng(2025)

for n in ns:
    samples = rng.choice(x, size=(reps, n), replace=True)
    xbars = samples.mean(axis=1)
    means_by_n[n] = xbars

    rmse_emp[n] = np.sqrt(np.mean((xbars - mu_hat)**2))

# RMSE teórico ~ sigma/sqrt(n)
rmse_teo = {n: sigma_hat/np.sqrt(n) for n in ns}

# -------------------------
# 4. Animación con Matplotlib
# -------------------------
fig, (ax_violin, ax_rmse) = plt.subplots(
    1, 2, figsize=(11, 4), gridspec_kw={"width_ratios": [2, 1]}
)

writer = FFMpegWriter(fps=0.4)

with writer.saving(fig, "consistencia_media_video.mp4", dpi=120):
    for k in range(len(ns)):
        current_ns = ns[:k+1]

        ax_violin.cla()
        data_violin = [means_by_n[n] for n in current_ns]
        positions   = np.arange(len(current_ns)) + 1

        parts = ax_violin.violinplot(
            data_violin,
            positions=positions,
            showmeans=False,
            showextrema=False,
            showmedians=False
        )

        for pc in parts['bodies']:
            pc.set_facecolor("lightgray")
            pc.set_edgecolor("gray")
            pc.set_alpha(0.8)

        # Boxplots por encima (mediana y cuartiles)
        ax_violin.boxplot(
            data_violin,
            positions=positions,
            widths=0.15,
            showfliers=True
        )

        # Línea horizontal en mu_hat (media "verdadera")
        ax_violin.axhline(mu_hat, linestyle="--", linewidth=1.2, color="red",
                          label=r"$\hat{\mu}$" if k == 0 else "")

        ax_violin.set_xticks(positions)
        ax_violin.set_xticklabels([str(n) for n in current_ns])
        ax_violin.set_xlabel("Tamaño muestral n")
        ax_violin.set_ylabel(r"$\bar{X}_n$")
        ax_violin.set_title(f"Consistencia de la media — {variable} (zona {zone})")

        # ----- Panel derecho: RMSE vs n -----
        ax_rmse.cla()
        # Empírico
        emp_x = np.array(current_ns)
        emp_y = np.array([rmse_emp[n] for n in current_ns])
        ax_rmse.plot(emp_x, emp_y, "o-", label="RMSE empírico")

        # Teórico
        teo_y = np.array([rmse_teo[n] for n in current_ns])
        ax_rmse.plot(emp_x, teo_y, "s--", label=r"$\sigma / \sqrt{n}$")

        ax_rmse.set_xscale("log")
        ax_rmse.set_yscale("log")
        ax_rmse.set_xlabel("n (escala log)")
        ax_rmse.set_ylabel("RMSE")
        ax_rmse.set_title("Disminución del error de la media")
        ax_rmse.legend()

        fig.suptitle("Ilustración de la consistencia de la media muestral", fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.92])

        writer.grab_frame()

plt.close(fig)
print("Video generado: consistencia_media_video.mp4")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from scipy.stats import t

def simulate_ic_for_n(x, n, conf=0.95, reps=200, rng=None):
    if rng is None:
        rng = np.random.default_rng(2025)
    alpha = 1 - conf
    mu_hat = x.mean()

    lo_list = []
    hi_list = []
    center_list = []
    cover_list = []

    for _ in range(reps):
        xi = rng.choice(x, size=n, replace=True)
        xb = xi.mean()
        s = xi.std(ddof=1)
        tcrit = t.ppf(1 - alpha/2, df=n - 1)
        hw = tcrit * s / np.sqrt(n)
        lo = xb - hw
        hi = xb + hw
        cover = (lo <= mu_hat) and (hi >= mu_hat)

        lo_list.append(lo)
        hi_list.append(hi)
        center_list.append(xb)
        cover_list.append(cover)

    return (np.array(lo_list), np.array(hi_list),
            np.array(center_list), np.array(cover_list), mu_hat)


def make_ci_video_fast(
    csv_path,
    variable="bp_rp",
    zone="A",
    n_small=30,
    n_large=200,
    conf=0.95,
    reps=200,
    filename="ic_video_fast.mp4",
    fps=0.5,
    frame_step=5
):
    # 1) Cargar datos
    AB = pd.read_csv(csv_path)
    x = AB.loc[AB["zone"] == zone, variable].astype(float)
    x = x.replace([np.inf, -np.inf], np.nan).dropna().values

    # 2) Simular IC para cada tamaño muestral
    rng = np.random.default_rng(2025)
    lo_s, hi_s, cen_s, cov_s, mu_hat = simulate_ic_for_n(
        x, n_small, conf=conf, reps=reps, rng=rng
    )
    lo_l, hi_l, cen_l, cov_l, _ = simulate_ic_for_n(
        x, n_large, conf=conf, reps=reps, rng=rng
    )

    # 3) Figura 2x2
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharex="col")
    ax_int_s, ax_cov_s = axes[0]
    ax_int_l, ax_cov_l = axes[1]

    writer = FFMpegWriter(fps=fps)
    frame_indices = list(range(frame_step, reps + 1, frame_step))

    with writer.saving(fig, filename, dpi=120):
        for k in frame_indices:
            ax_int_s.cla()
            ax_cov_s.cla()
            ax_int_l.cla()
            ax_cov_l.cla()

            # ---- n pequeño ----
            y_s = np.arange(1, k + 1)
            colors_s = np.where(cov_s[:k], "seagreen", "firebrick")
            # intervalos
            for i in range(k):
                ax_int_s.hlines(y=y_s[i],
                                xmin=lo_s[i],
                                xmax=hi_s[i],
                                colors=colors_s[i],
                                linewidth=0.8)
            ax_int_s.plot(cen_s[:k], y_s, "ko", markersize=2)
            ax_int_s.axvline(mu_hat, linestyle="--", color="black", linewidth=1.0)
            ax_int_s.set_ylabel("réplica")
            ax_int_s.set_title(f"IC (n = {n_small}, conf = {conf:.2f}) — zona {zone}")

            cov_run_s = np.cumsum(cov_s[:k]) / np.arange(1, k + 1)
            ax_cov_s.plot(np.arange(1, k + 1), cov_run_s, "-")
            ax_cov_s.axhline(conf, color="red", linestyle="--", label="nivel teórico")
            ax_cov_s.set_ylabel("cobertura acumulada")
            ax_cov_s.set_title("Cobertura (n pequeño)")
            ax_cov_s.legend(loc="lower right", fontsize=8)

            # ---- n grande ----
            y_l = np.arange(1, k + 1)
            colors_l = np.where(cov_l[:k], "seagreen", "firebrick")
            for i in range(k):
                ax_int_l.hlines(y=y_l[i],
                                xmin=lo_l[i],
                                xmax=hi_l[i],
                                colors=colors_l[i],
                                linewidth=0.8)
            ax_int_l.plot(cen_l[:k], y_l, "ko", markersize=2)
            ax_int_l.axvline(mu_hat, linestyle="--", color="black", linewidth=1.0)
            ax_int_l.set_xlabel(variable)
            ax_int_l.set_ylabel("réplica")
            ax_int_l.set_title(f"IC (n = {n_large}, conf = {conf:.2f}) — zona {zone}")

            cov_run_l = np.cumsum(cov_l[:k]) / np.arange(1, k + 1)
            ax_cov_l.plot(np.arange(1, k + 1), cov_run_l, "-")
            ax_cov_l.axhline(conf, color="red", linestyle="--", label="nivel teórico")
            ax_cov_l.set_xlabel("número de intervalos")
            ax_cov_l.set_ylabel("cobertura acumulada")
            ax_cov_l.set_title("Cobertura (n grande)")
            ax_cov_l.legend(loc="lower right", fontsize=8)

            fig.suptitle(
                f"IC para la media de {variable} — zona {zone} (conf = {conf:.2f})",
                fontsize=14
            )
            plt.tight_layout(rect=[0, 0, 1, 0.93])

            writer.grab_frame()

    plt.close(fig)
    print(f"Video generado: {filename}")


make_ci_video_fast(
    csv_path="AB.csv",
    variable="bp_rp",
    zone="A",
    n_small=30,
    n_large=200,
    conf=0.95,
    reps=500,
    filename="ic_bp_rp_zonaA_95_fast.mp4",
    fps=2,
    frame_step=7
)