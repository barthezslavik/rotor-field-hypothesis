#!/bin/bash
# Download GW event data from GWOSC

EVENT="GW231028_153006"
GPS=1382542224.3

echo "Downloading data for $EVENT (GPS: $GPS)"
echo "This will download H1 and L1 strain data around the event..."

# Download 32 seconds around the event (4096 Hz data)
wget -nc "https://gwosc.org/eventapi/json/$EVENT/v1/H-H1_GWOSC_4KHZ_R1-1382542216-32.hdf5" -P data/
wget -nc "https://gwosc.org/eventapi/json/$EVENT/v1/L-L1_GWOSC_4KHZ_R1-1382542216-32.hdf5" -P data/

echo "Download complete! Files saved to data/"
echo ""
echo "To analyze for sidebands, run:"
echo "python gw_sidebands.py \\"
echo "  --h1 data/H-H1_GWOSC_4KHZ_R1-1382542216-32.hdf5 \\"
echo "  --l1 data/L-L1_GWOSC_4KHZ_R1-1382542216-32.hdf5 \\"
echo "  --gps $GPS \\"
echo "  --pre 10 --post 10"
