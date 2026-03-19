# ──────────────────────────────────────────────────────────────────
# ProteinScope — CPU build (default)
#
# For GPU/CUDA acceleration use the CUDA build below or pass:
#   docker build --build-arg BASE=nvidia/cuda:12.4.1-cudnn9-runtime-ubuntu22.04 .
# ──────────────────────────────────────────────────────────────────
ARG BASE=python:3.13-slim
FROM ${BASE}

# Install system deps
RUN apt-get update && apt-get install -y \
    libxml2-dev libxslt-dev \
    libjpeg-dev libpng-dev \
    curl wget \
    && rm -rf /var/lib/apt/lists/*

# Python may not be present in CUDA base image — install if needed
RUN python3 --version 2>/dev/null || (apt-get update && apt-get install -y python3 python3-pip && rm -rf /var/lib/apt/lists/*)

WORKDIR /app

# Copy requirements first for layer caching
COPY proteinscope/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY proteinscope/ ./proteinscope/
WORKDIR /app/proteinscope

# Create data directories
RUN mkdir -p .data/pfam .data/esm2 diagrams output

EXPOSE 8000

CMD ["uvicorn", "web.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]


# ──────────────────────────────────────────────────────────────────
# CUDA build example (comment out the CPU ARG above and use this):
#
#   docker build \
#     --build-arg BASE=nvidia/cuda:12.4.1-cudnn9-runtime-ubuntu22.04 \
#     -t proteinscope:cuda .
#
# Then run with GPU access:
#   docker run --gpus all -p 8000:8000 proteinscope:cuda
#
# Or use docker-compose.gpu.yml which sets:
#   deploy.resources.reservations.devices:
#     - driver: nvidia
#       count: 1
#       capabilities: [gpu]
# ──────────────────────────────────────────────────────────────────
