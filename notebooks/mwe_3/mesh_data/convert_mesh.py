import meshio


def convert_med_to_xdmf(
    med_file: str,
    cell_file: str | None = "mesh_domains.xdmf",
    facet_file: str | None = "mesh_boundaries.xdmf",
    cell_type: str | None = "tetra",
    facet_type: str | None = "triangle",
):
    """Converts a .med mesh to .xdmf

    Args:
        med_file: the name of the MED file
        cell_file: the name of the file containing the volume markers. Defaults to "mesh_domains.xdmf".
        facet_file: the name of the file containing the surface markers.. Defaults to "mesh_boundaries.xdmf".
        cell_type: The topology of the cells. Defaults to "tetra".
        facet_type: The topology of the facets. Defaults to "triangle".

    Returns:
        dict, dict: the correspondence dict, the cell types
    """

    msh = meshio.read(med_file)

    correspondance_dict = {-k: v for k, v in msh.cell_tags.items()}

    cell_data_types = msh.cell_data_dict["cell_tags"].keys()

    for mesh_block in msh.cells:
        if mesh_block.type == cell_type:
            meshio.write_points_cells(
                cell_file,
                msh.points,
                [mesh_block],
                cell_data={"f": [-1 * msh.cell_data_dict["cell_tags"][cell_type]]},
            )
        elif mesh_block.type == facet_type:
            meshio.write_points_cells(
                facet_file,
                msh.points,
                [mesh_block],
                cell_data={"f": [-1 * msh.cell_data_dict["cell_tags"][facet_type]]},
            )

    return correspondance_dict, cell_data_types


med_file = "box_mesh.med"

corr_dict, cell_types = convert_med_to_xdmf(
    med_file=med_file, cell_type="triangle", facet_type="line"
)

print(corr_dict)
print(cell_types)
