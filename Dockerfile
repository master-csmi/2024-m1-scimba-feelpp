# Use Feel++ as the base image
FROM ghcr.io/feelpp/feelpp:jammy

LABEL maintainer="Helya Amiri <helya.amiri@etu.unistra.fr> , Rayen Tlili <rayen.tlili@etu.unistra.fr>"
LABEL description="Docker image with Feel++ and Scimba."

# Install additional dependencies for Scimba
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip



# Install Scimba directly
RUN pip3 install --no-cache-dir scimba

# Set the default command to run when the container starts
CMD ["python3", "-c", "import feelpp; import scimba"]
