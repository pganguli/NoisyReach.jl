#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

df = pd.DataFrame({
    "filtering_threshold": [0, 10, 20, 30, 40, 50, 60, 70, 80],
    "latency": [143, 138, 130, 118, 111, 92, 85, 67, 48],
    "accuracy": [90.4, 90.4, 90.3, 90.2, 90.2, 89.8, 86.5, 57.3, 34.2],
    "signal": [52.16, 48.11, 45.79, 39.83, 36.62, 27.12, 23.93, 57.74, 183.26],
    "std_dev_signal": [28.05, 27.83, 27.40, 24.48, 28.78, 22.70, 19.48, 79.31, 416.60],
})

# Set up the figure and axis
fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()
ax3 = ax1.twinx()
ax3.spines["right"].set_position(("outward", 60))  # Offset accuracy axis to the right

# Plot latency as bars with increased width and lighter color
bars = ax1.bar(
    df["filtering_threshold"],
    df["latency"],
    width=6,  # Thicker bars
    # color="#6acde5",
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

ax1_offset = [0,0,-30,0,0,0,0,0,-30]
# Annotate values on bars with more space
for i, bar in enumerate(bars):
    height = bar.get_height()
    offset = ax1_offset[i]
    ax1.annotate(
        f"{height}",
        xy=(bar.get_x() + bar.get_width() / 2, height),
        xytext=(0, offset),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=16,
        color="black",
    )

# Annotate values on line plots with adjusted spacing
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
ax2.set_ylim(0, 200)
ax3.tick_params(axis='y', labelsize=16)
ax3.set_ylim(0, 100)  
plt.subplots_adjust(hspace=0.9)

# plt.title("CASA: Latency, Accuracy, and RMSE tradeoff", fontsize=17)
plt.grid(True, which="both", linestyle="--", linewidth=0.9)
plt.tight_layout()
plt.show()
