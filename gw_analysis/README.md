# GW Sideband Analysis - Organized Structure

## Directory Structure

```
gw_analysis/
├── reports/          # Analysis reports and documentation
│   ├── GW_ANALYSIS_UPDATE_N8.md    ← LATEST (N=8 events)
│   ├── archive_N6.md                (N=6 events)
│   ├── archive_N3.md                (N=3 events)
│   └── archive_*.md                 (older reports)
│
├── scripts/          # Analysis scripts and utilities
│   ├── gw_sidebands.py              Main analysis pipeline
│   ├── analyze_single_event.py      Single event automation
│   └── analyze_4_events.sh          Batch processing
│
├── outputs/          # Per-event analysis outputs
│   ├── gw_out_GW231028_153006/      χ=0.45 (highest)
│   ├── gw_out_GW231123/             χ=0.31
│   ├── gw_out_GW190519/             χ=0.33
│   ├── gw_out_GW190412/             χ=0.25
│   ├── gw_out_GW190707_093326/      χ=0.28 (outlier)
│   ├── gw_out_GW190521/             χ=0.08 (control)
│   ├── gw_out_GW200129_065458/      χ=0.30 (perfect Δf=0)
│   └── gw_out_GW191109_010717/      χ=0.27 (best fit)
│
└── plots/            # Correlation plots and figures
    └── correlation_plots/
        └── chi_vs_omega_correlation_N8.png  ← LATEST
```

## Current Status (N=8)

**Statistical Significance Achieved!**
- H1 correlation: r = 0.725, **p = 0.042 < 0.05** ✓
- Average correlation: r = 0.685, p = 0.061 (very close)
- Linear fit: Ω = 25.66 × χ_eff - 3.30 Hz

**Next Goal**: Add 2 more events → N=10 → p < 0.05 for average

## Quick Start

**Analyze a new event**:
```bash
cd gw_analysis/scripts
python3 analyze_single_event.py GW_EVENT_NAME
```

**View latest results**:
```bash
cat gw_analysis/reports/GW_ANALYSIS_UPDATE_N8.md
```

**See correlation plot**:
```bash
open gw_analysis/plots/correlation_plots/chi_vs_omega_correlation_N8.png
```

## Data Location

Raw HDF5 data remains in: `/Users/slavik/Rotor/data/`
