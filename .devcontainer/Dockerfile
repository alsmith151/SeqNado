# Use official Ubuntu base image, forcing x86_64 (amd64) architecture for consistency across Apple Silicon and other builders
FROM --platform=linux/amd64 ubuntu:latest

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (compiler basics, networking, Apptainer runtime deps)
# Notes:
#  - squashfs-tools / squashfuse: needed for sandboxing and mounting images
#  - fuse3 & libfuse3-dev: runtime FUSE support
#  - uidmap, slirp4netns: enable unprivileged (user namespace) execution
#  - wget, curl, git: common utilities
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git ca-certificates gnupg lsb-release software-properties-common \
    build-essential pkg-config libseccomp-dev squashfs-tools squashfuse \
    fuse3 libfuse3-dev uidmap slirp4netns fakeroot fuse-overlayfs \
 && rm -rf /var/lib/apt/lists/*


# Install Apptainer
RUN add-apt-repository -y ppa:apptainer/ppa \
    && apt-get update && apt-get install -y apptainer-suid \
    && rm -rf /var/lib/apt/lists/*

# Install uv (copy static binaries from upstream image)
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

# Create non-root user for better security
RUN useradd -m -s /bin/bash seqnado_user && \
    mkdir -p /home/seqnado_user/.apptainer/cache /home/seqnado_user/workspace && \
    chown -R seqnado_user:seqnado_user /home/seqnado_user

# Switch to non-root user
USER seqnado_user
WORKDIR /home/seqnado_user/workspace

# Create virtual environment with uv and clone seqnado for editable install
RUN uv venv .venv

# Update PATH to include virtual environment explicitly (ensure /usr/bin present first for apptainer)
ENV PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/home/seqnado_user/workspace/.venv/bin:$PATH"

# Environment variables for Apptainer/Singularity
ENV APPTAINER_TMPDIR=/home/seqnado_user/.apptainer \
    APPTAINER_CACHEDIR=/home/seqnado_user/.apptainer/cache \
    SINGULARITY_TMPDIR=/home/seqnado_user/.apptainer \
    SINGULARITY_CACHEDIR=/home/seqnado_user/.apptainer/cache \
    APPTAINER_DISABLE_CACHE=false \
    SINGULARITY_DISABLE_CACHE=false \
    APPTAINER_WRITABLE_TMPFS=1

# Default command with helpful message about running Apptainer
CMD ["/bin/bash", "-c", "echo 'Apptainer Docker container ready. Use --privileged flag when running this container.'; echo 'Example: docker run --privileged -it <image_name>'; echo 'For unprivileged mode, use: apptainer --userns <command>'; /bin/bash"]