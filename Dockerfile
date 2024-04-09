# Start with the Feel++ base image
FROM ghcr.io/feelpp/feelpp:jammy

# Set labels for metadata
LABEL maintainer="Helya Amiri <helya.amiri@etu.unistra.fr>, Rayen Tlili <rayen.tlili@etu.unistra.fr>"
LABEL description="Docker image with Feel++, Scimba, and PyTorch."

# Install PyTorch.
RUN pip3 install torch

# Copy the local Scimba directory to the container.
COPY scimba/ /scimba/

# Install Scimba and its dependencies
WORKDIR /scimba
USER root

RUN pip3 install .

# Set the default command to launch Python.
CMD ["python3"]
