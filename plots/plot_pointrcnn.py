#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

df = pd.DataFrame({
    "filtering_threshold": [00, 10, 20, 30, 40, 50, 60, 70, 80],
    "latency": [124 ,118 ,113 ,98 ,87 ,76 ,71 ,62 ,56],
    "accuracy": [89.1 ,89.1 ,89.0 ,89.0 ,88.9 ,88.8 ,85.2 ,53.8 ,32.8],
    "signal": [39.43 ,36.61 ,36.14 ,27.46 ,22.99 ,19.07 ,19.58 ,64.85 ,9272.56],
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
    3, # Bar width
    #label="Latency (ms)",
    color="b",
    alpha=0.6,
)

# Plot accuracy as a line plot
ax3.plot(
    df["filtering_threshold"],
    df["accuracy"], "go-",
    #label="Accuracy (%)"
)

# Plot signal as line with error bars
ax2.errorbar(
    df["filtering_threshold"],
    df["signal"],
    yerr=df["std_dev_signal"],
    fmt="-o",
    color="r",
    #label="RMSE",
)

# Annotate values on bars
for bar in bars:
    height = bar.get_height()
    ax1.annotate(
        f"{height}",
        xy=(bar.get_x() + bar.get_width() / 2, height),
        xytext=(0, 5),  # Offset for better visibility
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=10,
        color="black",
    )

# Annotate values on line plots
for i, txt in enumerate(df["signal"]):
    ax2.annotate(
        f"{txt:.2f}",
        (df["filtering_threshold"][i], df["signal"][i]),
        textcoords="offset points",
        xytext=(0, 8),
        ha="center",
        fontsize=10,
        color="red",
    )
for i, txt in enumerate(df["accuracy"]):
    ax3.annotate(
        f"{txt}",
        (df["filtering_threshold"][i], df["accuracy"][i]),
        textcoords="offset points",
        xytext=(0, -12),
        ha="center",
        fontsize=10,
        color="green",
    )

# Labels and legend
ax1.set_xlabel("Filtering Threshold (%)")
ax1.set_ylabel("Latency (ms)", color="b")
ax2.set_ylabel("RMSE", color="r")
ax3.set_ylabel("Accuracy (%)", color="g")
ax1.set_xticks(df["filtering_threshold"])
#ax1.legend(loc="upper left")
#ax2.legend(loc="upper right")
#ax3.legend(loc="lower right")

plt.title("POINTRCNN: Latency, Accuracy, and RMSE tradeoff")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()
