#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

df = pd.DataFrame({
    "filtering_threshold": [00, 10, 20, 30, 40, 50, 60, 70, 80],
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

# Plot latency as bars
bars = ax1.bar(
    df["filtering_threshold"],
    df["latency"],
    3, # Bar width
    label="Latency (ms)",
    color="b",
    alpha=0.6,
)

# Plot accuracy as a line plot
ax3.plot(df["filtering_threshold"], df["accuracy"], "go-", label="Accuracy (%)")

# Plot signal as line with error bars
ax2.errorbar(
    df["filtering_threshold"],
    df["signal"],
    yerr=df["std_dev_signal"],
    fmt="-o",
    color="r",
    label="Signal",
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
ax2.set_ylabel("Signal (unitless)", color="r")
ax3.set_ylabel("Accuracy (%)", color="g")
ax1.set_xticks(df["filtering_threshold"])
ax1.legend(loc="upper left")
ax2.legend(loc="upper right")
ax3.legend(loc="lower right")

plt.title("Merged Chart of Latency, Accuracy, and Signal")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()
