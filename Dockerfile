# Use Feel++ as the base image
FROM ghcr.io/feelpp/feelpp:jammy

LABEL maintainer="Helya Amiri <helya.amiri@etu.unistra.fr> , Rayen Tlili <rayen.tlili@etu.unistra.fr>"
LABEL description="Docker image with Feel++ and Scimba."

ADD scimba/ /scimba/

WORKDIR /scimba

RUN pip3 install torch 
#RUN pip install .
# Set the default command to run when the container starts
# CMD ["python3", "-c", "import feelpp; import scimba"]
