#!/bin/bash
echo "===== Analyzing 4 events sequentially ====="

for event in GW200129_065458 GW190620_030421 GW190630_185205 GW191109_010717; do
    echo ""
    echo "Starting $event..."
    python3 analyze_single_event.py $event
    echo "✓ $event complete"
done

echo ""
echo "===== ALL 4 EVENTS COMPLETE ====="
