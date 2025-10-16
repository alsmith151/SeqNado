# Use official Ubuntu base image, forcing x86_64 (amd64) architecture for consistency across Apple Silicon and other builders
FROM --platform=linux/amd64 ubuntu:24.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies + pipx
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git build-essential python3 python3-pip python3-venv python3-full pipx \
    software-properties-common ca-certificates gnupg lsb-release \
    fuse3 libfuse3-dev fuse2fs squashfuse \
 && rm -rf /var/lib/apt/lists/*

# Install Apptainer (auto-detect architecture)
ARG APPTAINER_VERSION=1.4.1
RUN ARCH="$(dpkg --print-architecture)" && \
    wget -q "https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/apptainer_${APPTAINER_VERSION}_${ARCH}.deb" && \
    apt-get update && apt-get install -y ./apptainer_${APPTAINER_VERSION}_${ARCH}.deb || apt-get -f install -y && \
    rm -f apptainer_${APPTAINER_VERSION}_${ARCH}.deb

# Configure Apptainer for containerized environment
RUN mkdir -p /etc/apptainer && \
    { \
      echo "mount hostfs = no"; \
      echo "user bind control = yes"; \
      echo "allow setuid = no"; \
      echo "image driver = overlay"; \
      echo "enable overlay = try"; \
      echo "enable underlay = no"; \
    } >> /etc/apptainer/apptainer.conf

# Install seqnado with pipx
RUN pipx install seqnado && pipx ensurepath

# Create non-root user for better security
RUN useradd -m -s /bin/bash seqnado_user && \
    mkdir -p /home/seqnado_user/.apptainer/cache /home/seqnado_user/workspace && \
    chown -R seqnado_user:seqnado_user /home/seqnado_user

# Switch to non-root user
USER seqnado_user
WORKDIR /home/seqnado_user/workspace

# Environment variables for Apptainer/Singularity
ENV APPTAINER_TMPDIR=/home/seqnado_user/.apptainer \
    APPTAINER_CACHEDIR=/home/seqnado_user/.apptainer/cache \
    SINGULARITY_TMPDIR=/home/seqnado_user/.apptainer \
    SINGULARITY_CACHEDIR=/home/seqnado_user/.apptainer/cache \
    APPTAINER_DISABLE_CACHE=false \
    SINGULARITY_DISABLE_CACHE=false

# Default command
CMD ["/bin/bash"]
