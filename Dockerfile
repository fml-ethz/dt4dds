FROM python:3.10-slim

# Install OpenJDK-11 and wget
RUN apt-get update && \
    apt-get install -y default-jre wget && \
    apt-get clean;

# Create a directory for the application
WORKDIR /app

# Copy the requirements file into the container
COPY requirements.txt .

# Install the dependencies
RUN pip install dt4dds

# Install BBMap
RUN mkdir bbmap && wget https://sourceforge.net/projects/bbmap/files/BBMap_38.99.tar.gz/download && tar -zxvf download --strip-components=1 -C bbmap

# Symlink BBMap to ~/.local/bin/bbmap/
RUN mkdir -p ~/.local/bin/ && ln -s /app/bbmap ~/.local/bin/bbmap

# Create a directory for the data and go there
RUN mkdir -p /data
WORKDIR /data