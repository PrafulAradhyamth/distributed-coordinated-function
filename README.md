# Wi-Fi Saturation Throughput Analysis: Simulation vs. Analytical Model

---

This project provides a MATLAB script to analyze the **saturation throughput** of an IEEE 802.11 Distributed Coordination Function (DCF) wireless network. It offers two distinct approaches for calculating throughput:

1.  **Discrete-Event Simulation:** A simplified event-driven model that simulates the backoff process, transmissions, and collisions.
2.  **Bianchi's Analytical Model:** A theoretical model based on mathematical equations to predict network performance.

By combining both methods, the script allows for direct comparison and validation of the simulation results against the widely accepted analytical model.

## Table of Contents

* [Features](#features)
* [Getting Started](#getting-started)
    * [Prerequisites](#prerequisites)
    * [Running the Script](#running-the-script)
* [How it Works](#how-it-works)
    * [System Parameters](#system-parameters)
    * [Discrete-Event Simulation](#discrete-event-simulation)
    * [Bianchi's Analytical Model](#bianchis-analytical-model)
* [Interpreting the Results](#interpreting-the-results)
* [Parameters](#parameters)

---

## Features

* Calculates Wi-Fi DCF saturation throughput.
* Implements both a discrete-event simulation and Bianchi's analytical model.
* Allows for easy configuration of key Wi-Fi parameters (e.g., contention window, packet sizes, interframe spaces).
* Supports multiple simulation scenarios (different Contention Window settings).
* Generates a single plot for direct comparison of simulation and analytical results.
* Provides clear console output for ongoing progress.

---

## Getting Started

Follow these instructions to get a copy of the project up and running on your local machine.

### Prerequisites

* **MATLAB:** You need a working installation of MATLAB (R2018a or newer is recommended).

### Running the Script

1.  **Save the Code:**
    * Copy the provided MATLAB code (the combined script) into a file named `wifi_throughput_analysis.m`.
    * Save this file in your preferred working directory.
2.  **Open MATLAB:** Launch your MATLAB application.
3.  **Navigate:** In the MATLAB command window or file browser, navigate to the directory where you saved `wifi_throughput_analysis.m`.
4.  **Run:** Execute the script by typing the following in the MATLAB Command Window and pressing Enter:
    ```matlab
    wifi_throughput_analysis
    ```
    Alternatively, you can click the "Run" button in the MATLAB editor.

The script will execute, print progress to the command window, and then display a plot showing the saturation throughput curves for both simulation and the analytical model across different numbers of stations.

---

## How it Works

The script operates in two main parts for each configured scenario:

### System Parameters

The initial section defines all the core parameters of the Wi-Fi network based on typical 802.11 specifications, including:

* **`payload`**: Size of the data portion of a packet in bits.
* **`MACheader`**: Size of the MAC header in bits.
* **`PHYheader`**: Size of the Physical Layer header in bits.
* **`ack`**: Size of the Acknowledgment (ACK) frame in bits.
* **`ProDelay`**: Propagation delay in microseconds ($\mu s$).
* **`SIFS`**: Short Interframe Space in $\mu s$.
* **`SlotTime`**: Fundamental time unit for backoff in $\mu s$.
* **`DIFS`**: DCF Interframe Space in $\mu s$.

From these, the script calculates the durations of successful transmissions (`Ts`) and collisions (`Tc`) in terms of `SlotTime` units, which are crucial for Bianchi's model.

### Discrete-Event Simulation

This part of the code simulates the behavior of Wi-Fi stations over time:

* **Initialization:** Each station is assigned a random backoff time within its initial Contention Window (`CW`).
* **Event Loop:** The simulation progresses by finding the station with the minimum backoff time. This time represents an idle period on the channel.
* **Collision Detection:** If multiple stations reach a backoff of zero at the same time, a **collision** occurs.
    * Colliding stations double their `CW` (up to a maximum defined by `m`) and select a new random backoff time.
* **Successful Transmission:** If only one station transmits, it's a **success**.
    * The transmitting station resets its `CW` to the initial value and selects a new random backoff.
* **Time Advancement:** Simulation time advances by the duration of the idle period plus the time taken for the transmission (or collision).
* **Throughput Calculation:** The simulation runs until a fixed number of successful packets are transmitted, and the total time taken is used to calculate the throughput.

### Bianchi's Analytical Model

This part implements the well-known Bianchi's analytical model for 802.11 DCF saturation throughput:

* **Fixed Point Iteration:** The core of this model is finding the **conditional collision probability (`p`)** by solving a non-linear equation using `fzero`. This `p` represents the probability that a transmitted packet experiences a collision, given that at least one station is transmitting.
* **Transmission Probability (`$\tau$`):** Once `p` is determined, the probability that a station transmits in a randomly chosen slot time (`$\tau$`, or `tae` in the code) is calculated.
* **System State Probabilities:** Based on `$\tau$` and the number of stations, the probabilities of an idle slot, a successful transmission, or a collision in any given slot are derived.
* **Average Slot Time:** The average duration of a slot is calculated as a weighted sum of idle time, successful transmission time, and collision time.
* **Throughput Formula:** The final saturation throughput is then calculated using the average number of payload bits transmitted per average slot time.

---

## Interpreting the Results

The generated plot will visually compare the throughput obtained from the simulation and the analytical model for varying numbers of stations.

* **Simulation (Solid Lines):** These lines represent a more realistic (though simplified) approximation of network behavior. You might observe some minor fluctuations due to the random nature of the simulation.
* **Analytical Model (Dashed Lines):** These lines represent the theoretical maximum throughput under ideal conditions. They typically show a smoother curve and might be slightly higher than the simulation results, as the analytical model often makes simplifying assumptions.

Differences between the simulation and analytical model can highlight the impact of factors not fully captured by the analytical model (e.g., more complex backoff interactions, non-ideal channel conditions if they were to be introduced in a more advanced simulation).

---

## Parameters

You can modify the following parameters at the beginning of the `wifi_throughput_analysis.m` script to explore different scenarios:

* `payload`, `MACheader`, `PHYheader`, `ack`, `ProDelay`, `SIFS`, `SlotTime`, `DIFS`: Network timing and sizing parameters.
* `W_scenarios`: A vector defining different initial Contention Window (`CW_min + 1`) values to test.
* `m_scenarios`: A vector defining the corresponding maximum backoff stages for each `W_scenario`.
* `num_stations_range`: The range of the number of stations for which to calculate throughput.
* `Spkt_sim_target`: The number of successful packets to simulate for each run (determines simulation length).
