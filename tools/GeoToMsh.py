import gmsh

def geo_to_msh(geo_file, msh_file, mesh_size=0.1):
    gmsh.initialize()

    # Load the .geo file
    gmsh.open(geo_file)

    # Set the mesh size
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)

    # Generate the mesh
    gmsh.model.mesh.generate(2)  # For 2D mesh; use 3 for 3D mesh

    # Save the mesh to a .msh file
    gmsh.write(msh_file)

    gmsh.finalize()
