#!/bin/bash
set -e

# Ensure we're in the right directory
cd /private/groups/corbettlab/alex/tomG


echo "Running panmap with profiling..."

# Set up signal handling to ensure we get profiling data
cleanup() {
    echo "Caught signal - generating profile report..."
    cd /private/groups/corbettlab/alex/tomG
    gprof build/bin/panmap gmon.out > panmap_profile.txt
    echo "Profile saved to panmap_profile.txt"
    
    # Generate call graph visualization if gprof2dot is available
    if command -v gprof2dot &> /dev/null && command -v dot &> /dev/null; then
        gprof build/bin/panmap | gprof2dot | dot -Tpng -o panmap_callgraph.png
        gprof build/bin/panmap | gprof2dot | dot -Tsvg -o panmap_callgraph.svg
        echo "Call graph visualization saved as panmap_callgraph.png/svg"
    fi
    
    exit 1
}

trap cleanup SIGINT SIGTERM

# Run with a timeout to ensure we get profiling data after a fixed time
# Adjust the timeout value as needed (currently 60 seconds)
timeout --foreground 60s build/bin/panmap -f rsv_4K.panman || true

# Generate the profiling report
gprof build/bin/panmap gmon.out > panmap_profile.txt
echo "Profile saved to panmap_profile.txt"

# Generate call graph visualization if gprof2dot is available
if command -v gprof2dot &> /dev/null && command -v dot &> /dev/null; then
    gprof build/bin/panmap | gprof2dot | dot -Tpng -o panmap_callgraph.png
    gprof build/bin/panmap | gprof2dot | dot -Tsvg -o panmap_callgraph.svg
    echo "Call graph visualization saved as panmap_callgraph.png/svg"
fi

# Display the top functions by time
echo -e "\nTop functions by time:"
head -n 20 panmap_profile.txt

echo -e "\nDone. Full profile is in panmap_profile.txt"
