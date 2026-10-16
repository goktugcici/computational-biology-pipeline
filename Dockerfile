FROM python:3.11-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy source code and data
COPY src/ ./src/
COPY data/ ./data/
COPY reports/ ./reports/
COPY dashboard/ ./dashboard/

# Default entrypoint
CMD ["python", "src/main.py", "--generate_report", "--pan_cancer", "EGFR"]
