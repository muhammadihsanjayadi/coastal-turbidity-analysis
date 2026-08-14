# Coastal Turbidity Analysis

Automated ArcPy workflow for analysing coastal turbidity change from Sentinel-2 imagery.

This project examines changes in the optical turbidity signal of coastal waters adjacent to the Indonesia Morowali Industrial Park in Central Sulawesi, Indonesia. Sentinel-2 Multispectral Instrument imagery from 2 April 2016 and 6 April 2021 is processed through an automated ArcPy workflow to identify persistent water and calculate multi-date Normalized Difference Turbidity Index (NDTI) change.

## Final Map

![Coastal turbidity change near the Indonesia Morowali Industrial Park](figures/morowali_turbidity_change_map.png)

## Workflow

1. Applies Sentinel-2 Scene Classification Layer (SCL) masking to remove unsuitable observations.
2. Calculates the Normalized Difference Water Index (NDWI) to delineate water.
3. Intersects water masks from the two observation dates to identify persistent water.
4. Calculates the Normalized Difference Turbidity Index (NDTI) within identified water pixels.
5. Calculates NDTI change as later-date NDTI minus earlier-date NDTI.

The analysis identifies widespread positive NDTI change, with the strongest increases concentrated along sections of shoreline adjacent to the industrial park.

## Repository Structure

```text
coastal-turbidity-analysis/
├── src/
│   └── coastal_turbidity_analysis.py
├── figures/
│   └── morowali_turbidity_change_map.png
├── docs/
│   └── project_overview.pdf
├── .gitignore
└── README.md
