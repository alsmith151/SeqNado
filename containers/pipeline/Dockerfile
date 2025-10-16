FROM mambaorg/micromamba:latest

# Copy files
COPY environment_pipeline.yml /opt/environment_pipeline.yml
COPY multiqc_config.yaml /opt/seqnado/multiqc_config.yaml
COPY seqnado.png /opt/seqnado/seqnado.png

# Set environment variables
ENV PATH="/opt/conda/bin:$PATH"
ENV CONDA_DEFAULT_ENV=base
ENV CONDA_PREFIX="/opt/conda"
ENV LD_LIBRARY_PATH="/usr/lib:/lib"
ENV MAMBA_ROOT_PREFIX="/opt/conda"

# Ensure we have root for package installation and system-level setup
USER root

# Main build process
RUN set -e && \
    # Fix apt permissions and update
    mkdir -p /var/lib/apt/lists/partial && \
    chmod 755 /var/lib/apt/lists/partial && \
    apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      curl ca-certificates gcc libc6-dev pkg-config libssl-dev wget && \
    rm -rf /var/lib/apt/lists/* && \
    \
    # Install latest Quarto (architecture-aware)
    QUARTO_VERSION=$(curl -s https://api.github.com/repos/quarto-dev/quarto-cli/releases/latest | grep -o '"tag_name": "[^"]*' | cut -d'"' -f4 | sed 's/^v//') && \
    # Try to download a matching deb for the current architecture. prefer native arm64 when available.
    ARCH=$(dpkg --print-architecture 2>/dev/null || echo "$(uname -m)") && \
    echo "Detected architecture: ${ARCH}" && \
    QUARTO_TMP=/tmp/quarto.deb && \
    ( \
      if [ "${ARCH}" = "amd64" ] || [ "${ARCH}" = "x86_64" ]; then \
        urls=( "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" ); \
      else \
        # try common arm64/aarch64 filenames, then fall back to amd64
        urls=( \
          "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-aarch64.deb" \
          "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-arm64.deb" \
          "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
        ); \
      fi && \
      for url in "${urls[@]}"; do \
        echo "Trying $url"; \
        if curl -fSL --retry 3 --retry-delay 2 -o "${QUARTO_TMP}" "$url"; then \
          echo "Downloaded quarto from $url"; break; \
        fi; \
      done; \
    ) && \
    dpkg -i "${QUARTO_TMP}" && \
    rm -f "${QUARTO_TMP}" && \
    \
    # Micromamba environment setup (ensure Python <3.12)
    micromamba install -y -n base -f /opt/environment_pipeline.yml && \
    \
    # Install Rust
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y && \
    export PATH="/root/.cargo/bin:$PATH" && \
    cargo install bamnado && \
    \
    # Move bamnado binary to conda bin and cleanup Rust
    mv /root/.cargo/bin/bamnado /opt/conda/bin/ && \
    rm -rf /root/.cargo /root/.rustup && \
    \
    # Ensure MCCNado is installed. If not install it using pip
    if ! command -v mccnado &> /dev/null; then \
      /opt/conda/bin/python -m pip install mccnado; \
    fi && \
    \
    # Remove build tools
    apt-get purge -y wget gcc libc6-dev pkg-config libssl-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    \
    # Remove unnecessary data and dev files
    find /opt/conda/ -type f -name '*.a' -delete && \
    find /opt/conda/ -type f -name '*.pyc' -delete && \
    find /opt/conda/ -name '__pycache__' -type d -exec rm -rf {} + && \
    find /opt/conda/ -name '*.egg-info' -type d -exec rm -rf {} + && \
    find /opt/conda/ -name 'tests' -type d -exec rm -rf {} + && \
    rm -rf /root/.cache && \
    micromamba clean -ay

# Add labels
LABEL Author="asmith"
LABEL Version="0.4.0"
LABEL Description="Minimal SeqNado container with BamNado and Python tools"
LABEL BamNado.Version="dynamic"
LABEL BuildDate="$(date '+%Y-%m-%d')"

# Switch back to micromamba user (set by the base image)
USER $MAMBA_USER

# Health check / test
RUN echo "Testing minimal container..." && \
    echo "$(python --version)" && \
    echo "$(bamnado --help > /dev/null 2>&1 && echo "BamNado test passed" || exit 1)" && \
    echo "All tests passed!"

# Default command
ENTRYPOINT ["bash", "-c"]
CMD ["exec \"$@\""]