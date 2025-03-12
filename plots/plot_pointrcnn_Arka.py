#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.DataFrame({
    "filtering_threshold": [00, 10, 20, 30, 40, 50, 60, 70, 80],
    "latency": [124 ,118 ,113 ,98 ,87 ,76 ,71 ,62 ,56],
    "accuracy": [89.1 ,89.1 ,89.0 ,89.0 ,88.9 ,88.8 ,85.2 ,53.8 ,32.8],
    "signal": [39.43 ,36.61 ,36.14 ,27.46 ,22.99 ,19.07 ,19.58 ,64.85 ,np.nan],
    "std_dev_signal": [18.02 ,17.56 ,37.85 ,28.61 ,15.38 ,18.78 ,24.24 ,123.13 ,111608.50],
})

# Set up the figure and axis
fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("outward", 60))  # Offset accuracy axis to the right

# Plot latency as bars
bars = ax1.bar(
    df["filtering_threshold"],
    df["latency"],
    width=6, # Bar width
    #label="Latency (ms)",
    color="#80d4e9",
    alpha=0.7,
)

# Plot accuracy as a line plot
ax3.plot(
    df["filtering_threshold"],
    df["accuracy"], "o-",
    color="#19727F"
)

# Plot signal as line with error bars
ax2.plot(
    df["filtering_threshold"],
    df["signal"], "o-",
    # yerr=df["std_dev_signal"],
    # fmt="-o",
    color="#A80C28",
)

ax1_offset = [0,0,-35,0,0,0,0,-25,0]
# Annotate values on bars
for i, bar in enumerate(bars):
    height = bar.get_height()
    offset = ax1_offset[i]
    ax1.annotate(
        f"{height}",
        xy=(bar.get_x() + bar.get_width() / 2, height),
        xytext=(0, offset),  # Offset for better visibility
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=16,
        color="black",
    )

# Annotate values on line plots
for i, txt in enumerate(df["signal"]):
    ax2.annotate(
        f"{txt:.2f}",
        (df["filtering_threshold"][i], df["signal"][i]),
        textcoords="offset points",
        xytext=(0, 10),
        ha="center",
        fontsize=16,
        color="#A80C28",
    )
for i, txt in enumerate(df["accuracy"]):
    ax3.annotate(
        f"{txt}",
        (df["filtering_threshold"][i], df["accuracy"][i]),
        textcoords="offset points",
        xytext=(0, -18),
        ha="center",
        fontsize=16,
        color="#0A6A47",
    )

# Labels and legend
ax1.set_xlabel("Filtering Threshold (%)", fontsize=17)
ax1.set_ylabel("Latency (ms)", color="black", fontsize=17)
ax2.set_ylabel("Average RMSE", color="#d20a2e", fontsize=17)
ax3.set_ylabel("Accuracy (%)", color="#19727F", fontsize=17)
ax1.set_xticks(df["filtering_threshold"])
ax1.tick_params(axis='both', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.set_ylim(10, 70)
ax3.tick_params(axis='y', labelsize=16)
ax3.set_ylim(0, 100)  
plt.subplots_adjust(hspace=0.9)

# plt.title("POINTRCNN: Latency, Accuracy, and RMSE tradeoff")
plt.grid(True, which="both", linestyle="--", linewidth=0.9)
plt.tight_layout()
plt.show()
