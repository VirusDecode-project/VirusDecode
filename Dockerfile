# Use a multi-stage build to separate frontend, backend, and bioinformatics tools
# Stage 1: Base image with essential tools
FROM ubuntu:22.04 AS base

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    MUSCLE_VERSION=3.8.1551 \
    LINEARDESIGN_REPO=https://github.com/LinearDesignSoftware/LinearDesign.git

# Install general dependencies
RUN apt-get update && \
    apt-get install -y \
    openjdk-21-jdk \
    curl \
    wget \
    build-essential \
    python2 \
    python3.11 \
    python3-pip \
    git && \
    curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get install -y nodejs && \
    # npm install -g npm@latest && \
    npm install -g serve && \
    apt-get clean

# Install Python packages for bioinformatics
RUN pip3 install biopython==1.83 requests

# Install MUSCLE tool
RUN mkdir /muscle && cd /muscle && \
    wget https://www.drive5.com/muscle/muscle_src_${MUSCLE_VERSION}.tar.gz && \
    tar -xvzf muscle_src_${MUSCLE_VERSION}.tar.gz && \
    make && \
    cp muscle /usr/local/bin/ && \
    cd .. && rm -rf /muscle

# Install LinearDesign for bioinformatics analysis
RUN git clone ${LINEARDESIGN_REPO} /LinearDesign && \
    cd /LinearDesign && \
    make && \
    cd ..

# Stage 2: Build the backend (Spring Boot)
FROM base AS backend

# Set the working directory to the backend folder
WORKDIR /VirusDecode/backend

# Copy the backend files from the host machine
COPY ./backend /VirusDecode/backend

# Grant executable permissions for Gradle wrapper
RUN chmod +x ./gradlew

# Build the Spring Boot application
RUN ./gradlew build

# Stage 3: Build the frontend (React)
FROM base AS frontend

# Set the working directory to the frontend folder
WORKDIR /VirusDecode/frontend

# Copy the frontend files from the host machine
COPY ./frontend /VirusDecode/frontend

# Install Node.js dependencies
RUN npm install

# Build the React frontend
RUN npm run build

# Stage 4: Final stage with runtime environment
FROM base AS final

# Create desired directory structure
WORKDIR /

# Copy backend and frontend builds
COPY --from=backend /VirusDecode/backend/build/libs/virusdecode.jar /VirusDecode/backend/
COPY --from=frontend /VirusDecode/frontend/build /VirusDecode/frontend/build

# Copy Python scripts and LinearDesign into final image
COPY ./virusdecode.py /VirusDecode/
COPY --from=base /LinearDesign /LinearDesign

# Expose necessary ports
EXPOSE 8080 3000
