#!/bin/bash
set -e

echo "Detecting platform..."

# Get user platform
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="linux"
    echo "Detected: Linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="darwin"
    echo "Detected: macOS/Darwin"
else
    echo "Unsupported platform: $OSTYPE"
    exit 1
fi

# Add "#" in front of "platform: linux/amd64" in the docker-compose.yml file if linux
# Comment out platform line on Linux, uncomment on Mac
if [[ "$PLATFORM" == "linux" ]]; then
    echo "Configuring for Linux (native x86_64)..."
    sed -i 's/^\([[:space:]]*\)platform: linux\/amd64/#\1platform: linux\/amd64/' docker-compose.yml
    echo "✓ Platform line commented out for native performance"
else
    echo "Configuring for macOS (emulation required for InterProScan)..."
    sed -i '' 's/^[[:space:]]*#[[:space:]]*platform: linux\/amd64/    platform: linux\/amd64/' docker-compose.yml
    echo "✓ Platform line enabled for x86_64 emulation"
fi

echo "docker-compose.yml configured successfully!"
