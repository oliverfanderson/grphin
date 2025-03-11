from matplotlib import pyplot as plt
import numpy as np


def plot_runtime_stats():
    species = ["bsub", "cerevisiae", "drerio", "elegans", "fly"]
    runtime_full_algorithm_data = [6.132, 4507.445, 407.541, 303.814, 340.379]
    runtime_triplet_iteration_data = [0.461, 194.703, 51.803, 24.740, 29.512]

    stacked_bar_data = {
        "below": runtime_full_algorithm_data,
        "above": runtime_triplet_iteration_data,
    }

    fig, ax = plt.subplots()
    width = 0.5
    bottom = np.zeros(len(species))

    for stack, data in stacked_bar_data.items():
        p = ax.bar(species, data, width=width, label=stack, bottom=bottom)
        bottom += data

    plt.show()
