# Dockerfile for SeqNado ML CPU
FROM python:3.12-bookworm

LABEL author="asmith"
LABEL software="SeqNado"
LABEL version="1.0"

# Copy repositories
COPY ../repos/lanceotron-mcc /opt/repos/lanceotron-mcc

# Upgrade pip, setuptools, and wheel
RUN python -m pip install --no-cache-dir --upgrade pip setuptools wheel

# Install PyTorch and related packages from PyTorch index
RUN python -m pip install --no-cache-dir --upgrade --no-input \
    torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# Install other Python dependencies
RUN python -m pip install --no-cache-dir \
    loguru \
    pandas \
    "polars<1.9" \
    pydantic \
    pyBigWig \
    tqdm \
    typer \
    cyclopts \
    joblib \
    pyranges \
    crested \
    modisco-lite \
    natsort \
    tensorflow-cpu

# Install local packages
RUN python -m pip install --no-cache-dir --no-deps \
    /opt/repos/quantnado \
    /opt/repos/lanceotron-mcc \
    lanceotron

# Clean up to reduce image size
RUN rm -rf /opt/repos && \
    python -m pip cache purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Default runscript
CMD ["/bin/bash"]

# Healthcheck (optional, basic import test)
HEALTHCHECK --interval=1m --timeout=10s --start-period=30s --retries=3 \
  CMD python -c "import lanceotron, quantnado, torch, tensorflow, lanceotron_mcc" || exit 1
