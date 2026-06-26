import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------
# CONFIG
# --------------------------------------------------

filepath = "/Users/saenz-arevalo/ILL/CSPEC/CT2_april_2026/neutron_detectors/data/"
run_number = 46067
filename = filepath+f"0{run_number}.csv"

N_META = 10
WF_SIZE = 1800
TOTAL_WF = WF_SIZE * 2  # 3600 expected
N_SAMPLES_BASELINE = 100  # first 100 samples for baseline estimation
ANALYSIS_CHANNEL = [27, 28, 29, 30, 31]  # example channel to analyze

# --------------------------------------------------
# STORAGE
# --------------------------------------------------

metadata_list = []
wf1_list = []
wf2_list = []

bad_rows = 0

print("Reading:", filename)

# --------------------------------------------------
# PARSING LOOP (robust ASCII ingestion)
# --------------------------------------------------

with open(filename, "r") as f:

    for line_number, line in enumerate(f, start=1):

        fields = line.strip().split(";")

        # remove trailing empty field from final ;
        if fields and fields[-1] == "":
            fields = fields[:-1]

        # basic sanity check
        if len(fields) < N_META:
            bad_rows += 1
            continue

        meta_fields = fields[:N_META]
        wf_fields = fields[N_META:]

        # waveform validation
        if len(wf_fields) != TOTAL_WF:
            bad_rows += 1
            continue

        try:
            meta = np.array(meta_fields, dtype=np.int64)
            wf = np.array(wf_fields, dtype=np.int16)

        except ValueError:
            bad_rows += 1
            continue

        metadata_list.append(meta)
        wf1_list.append(wf[:WF_SIZE])
        wf2_list.append(wf[WF_SIZE:])

# --------------------------------------------------
# BUILD STRUCTURED EVENT ARRAY
# --------------------------------------------------

n_events = len(metadata_list)

events = np.zeros(
    n_events,
    dtype=[
        ("meta", np.int64, (N_META,)),
        ("wf1", np.int16, (WF_SIZE,)),
        ("wf2", np.int16, (WF_SIZE,))
    ]
)

events["meta"] = np.array(metadata_list)
events["wf1"] = np.array(wf1_list, dtype=np.int16)
events["wf2"] = np.array(wf2_list, dtype=np.int16)

print("\n--- SUMMARY ---")
print("Loaded events:", n_events)
print("Rejected rows:", bad_rows)
print("Event dtype:", events.dtype)

# --------------------------------------------------
# Plotting histogram of counts per channel for diagnostics
# --------------------------------------------------
channels = events["meta"][:, 3]
plt.figure(figsize=(8,5))
plt.hist(channels, bins=np.arange(-0.5, 32.5, 1), alpha=0.7)
plt.xlabel("Channel")
plt.ylabel("Count")
plt.title("Channel Distribution")
plt.grid()
plt.show()

# --------------------------------------------------
# METADATA CUT EXAMPLE
# --------------------------------------------------

# meta is composed as follows:
# [0] timestamp
# [1] acq_board (0 @ CT2)
# [2] adc_borad (0 @ CT2)
# [3] channel
# [4] event_type (6 @ CT2 for QDIV)
# [5] A
# [6] A+B
# [7] flags (not useful for now)
# [8] position = A/(A+B)*2**16 (16-bit fixed point)
# [9] number of samples (should be TOTAL_WF)

# Example cuts:
# channel == 3 AND event_type == 6


# create multiple channel masks
masks = []
for ch in ANALYSIS_CHANNEL:
    mask = (events["meta"][:, 3] == ch) & (events["meta"][:, 4] == 6)
    masks.append(mask)

# create different event subsets for each channel
cut_events_list = []
for mask in masks:
    cut_events = events[mask]
    cut_events_list.append(cut_events)
    print(f"Channel {cut_events['meta'][0, 3]}: {len(cut_events)} events after cut")



# --------------------------------------------------
# Plotting first 16 raw waveforms to check data quality 
# --------------------------------------------------

for cut_events in cut_events_list:
    # 16 subplots in a 4x4 grid
    plt.figure(figsize=(12, 10))
    # main title with channel info
    plt.suptitle(f"Raw Waveforms for Channel {cut_events['meta'][0, 3]}", fontsize=16)
    for i in range(min(16, len(cut_events))):
        # subplot
        plt.subplot(4, 4, i+1)
        plt.plot(cut_events["wf1"][i], label="WF1")
        plt.plot(cut_events["wf2"][i], label="WF2")
        plt.legend()
        if i >= 12:  # only label x-axis on bottom row
            plt.xlabel("Sample")
        if i % 4 == 0:  # only label y-axis on left column
            plt.ylabel("ADC")
        plt.grid()

    plt.show()

# --------------------------------------------------
# BASELINE CORRECTION
# --------------------------------------------------

def baseline_correct(wf):
    baseline = wf[:, :N_SAMPLES_BASELINE].mean(axis=1)
    return wf - baseline[:, None]

# Plotting histogram of baseline values for diagnostics
baseline_values = []
for cut_events in cut_events_list:
     baseline_values.append(cut_events["wf1"][:, :N_SAMPLES_BASELINE].mean(axis=1))
     baseline_values.append(cut_events["wf2"][:, :N_SAMPLES_BASELINE].mean(axis=1))

# plot baseline distributions for analysis channels
plt.figure(figsize=(15,8))
plt.suptitle("Baseline ADC Value Distribution", fontsize=16)
for i, cut_events in enumerate(cut_events_list):
    plt.subplot(2, 3, i+1)
    plt.hist(baseline_values[2*i], bins=np.arange(0, 200, 5), alpha=0.7, label="WF1")
    plt.hist(baseline_values[2*i+1], bins=np.arange(0, 200, 5), alpha=0.7, label="WF2")
    plt.xlabel("ADC")
    plt.ylabel("Count")
    plt.title(f"Ch. {cut_events['meta'][0, 3]}")
    plt.legend()
    plt.grid()
plt.show()

# --------------------------------------------------
# Plotting histogram of A+B values for diagnostics
# --------------------------------------------------
plt.figure(figsize=(15,10))
plt.suptitle("A+B ADC Value Distribution", fontsize=16)

linestyles = ['-', '--', '-.', ':']

for i, cut_events in enumerate(cut_events_list):

    a_plus_b_values = cut_events["meta"][:,6]

    n_bins = 200 # number of bins for histogram 
    counts, bins = np.histogram(
        a_plus_b_values,
        
        bins=np.arange(0, 50000, int(50000/n_bins) )
    )

    # normalize per channel
    counts = counts / counts.sum()
    channel = cut_events["meta"][0,3]

    plt.step(
        bins[:-1],
        counts,
        where="post",
        label=f"Ch. {channel}",
        linestyle=linestyles[i % len(linestyles)],
        linewidth=2
    )

    print()

plt.xlabel("A+B")
plt.ylabel("Normalized counts")
plt.grid()
plt.legend()
plt.show()


wf1_corr_list = []
wf2_corr_list = []
for cut_events in cut_events_list:  
    wf1_corr = baseline_correct(cut_events["wf1"])
    wf2_corr = baseline_correct(cut_events["wf2"])
    wf1_corr_list.append(wf1_corr)
    wf2_corr_list.append(wf2_corr)



# --------------------------------------------------
# AVERAGE WAVEFORMS
# --------------------------------------------------
avg_wf1_list = []
avg_wf2_list = []
for wf1_corr, wf2_corr in zip(wf1_corr_list, wf2_corr_list):
    avg_wf1 = wf1_corr.mean(axis=0)
    avg_wf2 = wf2_corr.mean(axis=0)
    avg_wf1_list.append(avg_wf1)
    avg_wf2_list.append(avg_wf2)



# plotting average waveforms 
plt.figure(figsize=(12,8))
plt.suptitle("Average Waveforms After Baseline Correction", fontsize=16)

colors = np.vstack([
    plt.cm.tab20(np.linspace(0,1,20)),
    plt.cm.Dark2(np.linspace(0,1,8)),
    plt.cm.Set1(np.linspace(0,1,9))
])

for i, (avg_wf1, avg_wf2) in enumerate(zip(avg_wf1_list, avg_wf2_list)):

    channel = cut_events_list[i]["meta"][0, 3]
    color = colors[i % len(colors)]

    # WF1
    plt.plot(
        avg_wf1,
        color=color,
        linestyle="-",
        linewidth=2,
        label=f"Ch.{channel} WF1"
    )

    # WF2
    plt.plot(
        avg_wf2,
        color=color,
        linestyle="--",
        linewidth=2,
        label=f"Ch.{channel} WF2"
    )

plt.xlabel("Sample")
plt.ylabel("ADC (baseline corrected)")
plt.grid()

plt.legend(
    ncol=2,
    fontsize=8
)

plt.show()

# --------------------------------------------------
# TIME OVER THRESHOLD
# --------------------------------------------------

def tot(wf, thr):
    return np.sum(wf > thr, axis=1)

thresholds = np.arange(50, 500, 10)

tot1_list = []
tot2_list = []

for wf1_corr, wf2_corr in zip(wf1_corr_list, wf2_corr_list):

    tot1_mean = []
    tot2_mean = []

    for thr in thresholds:
        tot1 = tot(wf1_corr, thr)
        tot2 = tot(wf2_corr, thr)
        tot1_mean.append(tot1.mean())
        tot2_mean.append(tot2.mean())
    
    tot1_list.append(tot1_mean)
    tot2_list.append(tot2_mean)


# plotting mean ToT vs threshold
plt.figure(figsize=(12,8))
plt.suptitle("Mean Time Over Threshold vs Threshold", fontsize=16)

# Build a large pool of very different colors
colors = np.vstack([
    plt.cm.tab20(np.linspace(0,1,20)),
    plt.cm.Dark2(np.linspace(0,1,8)),
    plt.cm.Set1(np.linspace(0,1,9))
])

for i, (tot1_mean, tot2_mean) in enumerate(zip(tot1_list, tot2_list)):

    channel = cut_events_list[i]["meta"][0,3]

    color = colors[i % len(colors)]

    # A waveform
    plt.plot(
        thresholds,
        tot1_mean,
        color=color,
        linestyle="-",
        linewidth=2,
        label=f"Ch.{channel} A"
    )

    # B waveform
    plt.plot(
        thresholds,
        tot2_mean,
        color=color,
        linestyle="--",
        linewidth=2,
        label=f"Ch.{channel} B"
    )

plt.xlabel("Threshold")
plt.ylabel("Mean ToT")
plt.grid()

plt.legend(
    ncol=2,
    fontsize=8
)

plt.show()

# --------------------------------------------------
# OPTIONAL: QUICK DIAGNOSTICS
# --------------------------------------------------

print("\n--- DIAGNOSTICS ---")
print("WF1 shape:", events["wf1"].shape)
print("WF2 shape:", events["wf2"].shape)
print("Metadata shape:", events["meta"].shape)

start_sample = 350
# Get RMS from each waveform starting from sample start_sample for diagnostics
rms_list_wf1 = []
rms_list_wf2 = []
for wf1_corr, wf2_corr in zip(wf1_corr_list, wf2_corr_list):
    rms1 = np.sqrt(np.mean(wf1_corr[:, start_sample:]**2, axis=1))
    rms2 = np.sqrt(np.mean(wf2_corr[:, start_sample:]**2, axis=1))
    rms_list_wf1.append(rms1)
    rms_list_wf2.append(rms2)

# Plotting histogram of RMS values for diagnostics
plt.figure(figsize=(15,8))
plt.suptitle(f"RMS Value Distribution (from sample {start_sample} onwards)", fontsize=16)
for i, cut_events in enumerate(cut_events_list):
    plt.subplot(2, 3, i+1)
    plt.hist(rms_list_wf1[i], bins=np.arange(0, 10, 0.2), alpha=0.7, label="WF1")
    plt.hist(rms_list_wf2[i], bins=np.arange(0, 10, 0.2), alpha=0.7, label="WF2")
    plt.xlabel("RMS")
    plt.ylabel("Count")
    plt.title(f"Ch. {cut_events['meta'][0, 3]}")
    plt.legend()
    plt.grid()
plt.show() 


# Get maximum amplitude of average waveforms for each analysis channel
max_amp_list = []
for avg_wf1, avg_wf2 in zip(avg_wf1_list, avg_wf2_list):
    max_amp1 = avg_wf1.max()
    max_amp2 = avg_wf2.max()
    max_amp_list.append((max_amp1, max_amp2))

# normalize average waveforms by their maximum amplitude for better comparison
avg_wf1_norm_list = []
avg_wf2_norm_list = []
for avg_wf1, avg_wf2, (max_amp1, max_amp2) in zip(avg_wf1_list, avg_wf2_list, max_amp_list):
    avg_wf1_norm = avg_wf1 / max_amp1 if max_amp1 != 0 else avg_wf1
    avg_wf2_norm = avg_wf2 / max_amp2 if max_amp2 != 0 else avg_wf2
    avg_wf1_norm_list.append(avg_wf1_norm)
    avg_wf2_norm_list.append(avg_wf2_norm)


# plot average waveform after normalization to check shape differences
plt.figure(figsize=(12,8))
plt.suptitle("Normalized Average Waveforms", fontsize=16)
colors = np.vstack([
    plt.cm.tab20(np.linspace(0,1,20)),
    plt.cm.Dark2(np.linspace(0,1,8)),
    plt.cm.Set1(np.linspace(0,1,9))
])

for i, (avg_wf1_norm, avg_wf2_norm) in enumerate(zip(avg_wf1_norm_list, avg_wf2_norm_list)):
    
    channel = cut_events_list[i]["meta"][0, 3]
    color = colors[i % len(colors)]

    # WF1
    plt.plot(
        avg_wf1_norm,
        color=color,
        linestyle="-",
        linewidth=2,
        label=f"Ch.{channel} WF1"
    )

    # WF2
    plt.plot(
        avg_wf2_norm,
        color=color,
        linestyle="--",
        linewidth=2,
        label=f"Ch.{channel} WF2"
    )
plt.xlabel("Sample")
plt.ylabel("Normalized ADC")
plt.grid()
plt.legend(
    ncol=2,
    fontsize=8
)
plt.show()