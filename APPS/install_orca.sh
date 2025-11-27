#!/usr/bin/env bash
#
# ORCA 5.x Installation Script
# ORCA is free for academic use but requires manual download
#

set -e

ORCA_DIR="/Users/stevens/Dropbox/Backplane/APPS/orca"
ORCA_VERSION="5.0.4"

echo "========================================"
echo "ORCA 5.x Installation Script"
echo "========================================"
echo

# Check if ORCA already installed
if [ -f "$ORCA_DIR/orca" ]; then
    echo "✓ ORCA already installed at: $ORCA_DIR"
    $ORCA_DIR/orca --version 2>&1 | head -10
    exit 0
fi

echo "ORCA requires manual download (free for academic use)"
echo
echo "Download instructions:"
echo "1. Visit: https://orcaforum.kofo.mpg.de/"
echo "2. Register for free academic account"
echo "3. Download: ORCA 5.0.4 macOS ARM64 (Apple Silicon)"
echo "4. Place downloaded file in: /Users/stevens/Dropbox/Backplane/APPS/"
echo

# Check for downloaded file
ORCA_ARCHIVE=$(ls /Users/stevens/Dropbox/Backplane/APPS/orca_*.tar.* 2>/dev/null | head -1)

if [ -z "$ORCA_ARCHIVE" ]; then
    echo "❌ ORCA archive not found in APPS directory"
    echo
    echo "Expected filename pattern: orca_5_0_4_*.tar.xz or orca_5_0_4_*.tar.gz"
    echo
    echo "After downloading, run this script again:"
    echo "  bash install_orca.sh"
    exit 1
fi

echo "Found ORCA archive: $ORCA_ARCHIVE"
echo "Extracting..."

# Extract archive
cd /Users/stevens/Dropbox/Backplane/APPS/
if [[ "$ORCA_ARCHIVE" == *.tar.xz ]]; then
    tar -xf "$ORCA_ARCHIVE"
elif [[ "$ORCA_ARCHIVE" == *.tar.gz ]]; then
    tar -xzf "$ORCA_ARCHIVE"
fi

# Find extracted directory
ORCA_EXTRACTED=$(ls -d orca_*_*/ 2>/dev/null | head -1)
if [ -z "$ORCA_EXTRACTED" ]; then
    echo "❌ Could not find extracted ORCA directory"
    exit 1
fi

# Rename to standard location
mv "$ORCA_EXTRACTED" orca

echo "✓ ORCA extracted to: $ORCA_DIR"

# Make binaries executable
chmod +x $ORCA_DIR/orca
chmod +x $ORCA_DIR/orca_*

# Test ORCA
echo
echo "Testing ORCA installation..."
$ORCA_DIR/orca --version 2>&1 | head -10

echo
echo "========================================"
echo "✅ ORCA installation complete!"
echo "========================================"
echo
echo "ORCA binary: $ORCA_DIR/orca"
echo "Version: $($ORCA_DIR/orca --version 2>&1 | grep -i version | head -1)"
echo
echo "Add to PATH (optional):"
echo "  export PATH=\$PATH:$ORCA_DIR"
echo
