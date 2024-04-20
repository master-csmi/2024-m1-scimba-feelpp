# Use Feel++ as the base image
FROM ghcr.io/feelpp/feelpp:jammy

LABEL maintainer="Helya Amiri <helya.amiri@etu.unistra.fr> , Rayen Tlili <rayen.tlili@etu.unistra.fr>"
LABEL description="Docker image with Feel++ and Scimba."

# Update the system and install additional dependencies for Scimba and Git
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    git

    # Install PyTorch
RUN pip3 install torch

# Check git version
RUN git --version

# Install Scimba directly
# The --no-cache-dir option is used to disable the cache during the installation, which can reduce the size of your Docker image.
RUN pip3 install --no-cache-dir scimba

# Set the default command to run when the container starts
CMD ["python3", "-c", "import feelpp; import scimba"]
