#!/bin/bash
# setup_conda.sh
# Install Miniforge in $HOME and configure it for Snakemake
# Run this once to set up a working conda installation

set -euo pipefail

echo "=========================================="
echo "Miniforge Installation Script"
echo "=========================================="
echo ""

# Check if conda already exists and works
if command -v conda >/dev/null 2>&1; then
  if python3 -c "import conda" 2>/dev/null; then
    echo "[INFO] Working conda found: $(which conda)"
    echo "[INFO] Conda version: $(conda --version)"
    echo "[INFO] You may not need to install Miniforge."
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      exit 0
    fi
  else
    echo "[WARNING] Conda found but Python module is broken (likely Spack conda)"
    echo "[INFO] Will install Miniforge to fix this"
  fi
fi

# Installation directory
INSTALL_DIR="${HOME}/miniforge3"
INSTALLER="${HOME}/Miniforge3.sh"

echo "[INFO] Installation directory: $INSTALL_DIR"
echo "[INFO] Installer will be downloaded to: $INSTALLER"
echo ""

# Download Miniforge if not exists
if [ ! -f "$INSTALLER" ]; then
  echo "[INFO] Downloading Miniforge3..."
  curl -L -o "$INSTALLER" \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  chmod +x "$INSTALLER"
else
  echo "[INFO] Installer already exists: $INSTALLER"
fi

# Install Miniforge
if [ -d "$INSTALL_DIR" ]; then
  echo "[WARNING] $INSTALL_DIR already exists"
  read -p "Remove and reinstall? (y/N): " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm -rf "$INSTALL_DIR"
  else
    echo "[INFO] Skipping installation"
    exit 0
  fi
fi

echo "[INFO] Installing Miniforge to $INSTALL_DIR..."
bash "$INSTALLER" -b -p "$INSTALL_DIR"

# Initialize conda
echo "[INFO] Initializing conda for $(basename $SHELL)..."
"$INSTALL_DIR/bin/conda" init "$(basename $SHELL)"

# Verify installation
echo ""
echo "[INFO] Verifying installation..."
"$INSTALL_DIR/bin/conda" --version
python3 -c "import conda; print('Conda module:', conda.__file__)" || {
  echo "[ERROR] Conda Python module not found after installation"
  exit 1
}

echo ""
echo "=========================================="
echo "Installation complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Restart your shell or run: source ~/.$(basename $SHELL)rc"
echo "2. Verify: conda --version"
echo "3. Create Snakemake env: conda env create -f envs/smk.yaml"
echo "4. Activate: conda activate smk"
echo ""
echo "To use with Snakemake, ensure PATH includes:"
echo "  $INSTALL_DIR/bin"
echo ""

