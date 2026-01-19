function create_cuboid_grid_with_quantum_well(
        cell_up,
        cell_qw,
        cell_low,
        boundary_bottom;
        L_up = 6.0,
        L_mid = 4.0,
        L_low = 10.0,
        L_width = 50.0,
        h = 1.0,
        H = 2.0
    )

    builder = SimplexGridBuilder(; Generator = TetGen)

    # bottom
    p01 = point!(builder, 0, 0, 0)
    p02 = point!(builder, L_width, 0, 0)
    p03 = point!(builder, L_width, L_width, 0)
    p04 = point!(builder, 0, L_width, 0)
    facetregion!(builder, boundary_bottom)
    facet!(builder, p01, p02, p03, p04)

    # lower intersection
    p11 = point!(builder, 0, 0, L_low)
    p12 = point!(builder, L_width, 0, L_low)
    p13 = point!(builder, L_width, L_width, L_low)
    p14 = point!(builder, 0, L_width, L_low)
    facetregion!(builder, 2)
    facet!(builder, p11, p12, p13, p14)

    # upper intersection
    p21 = point!(builder, 0, 0, L_low + L_mid)
    p22 = point!(builder, L_width, 0, L_low + L_mid)
    p23 = point!(builder, L_width, L_width, L_low + L_mid)
    p24 = point!(builder, 0, L_width, L_low + L_mid)
    facetregion!(builder, 3)
    facet!(builder, p21, p22, p23, p24)

    # top
    p31 = point!(builder, 0, 0, L_low + L_mid + L_up)
    p32 = point!(builder, L_width, 0, L_low + L_mid + L_up)
    p33 = point!(builder, L_width, L_width, L_low + L_mid + L_up)
    p34 = point!(builder, 0, L_width, L_low + L_mid + L_up)
    facetregion!(builder, 4)
    facet!(builder, p31, p32, p33, p34)

    # front face
    facetregion!(builder, 5)
    facet!(builder, p01, p02, p12, p11)
    facet!(builder, p11, p12, p22, p21)
    facet!(builder, p21, p22, p32, p31)

    # right face
    facetregion!(builder, 6)
    facet!(builder, p02, p03, p13, p12)
    facet!(builder, p12, p13, p23, p22)
    facet!(builder, p22, p23, p33, p32)

    # back face
    facetregion!(builder, 7)
    facet!(builder, p03, p04, p14, p13)
    facet!(builder, p13, p14, p24, p23)
    facet!(builder, p23, p24, p34, p33)

    # left face
    facetregion!(builder, 8)
    facet!(builder, p04, p01, p11, p14)
    facet!(builder, p14, p11, p21, p24)
    facet!(builder, p24, p21, p31, p34)

    # lower cell region
    cellregion!(builder, cell_low)
    maxvolume!(builder, H)
    regionpoint!(builder, L_width / 2, L_width / 2, L_low / 2)

    # quantum well cell region
    cellregion!(builder, cell_qw)
    maxvolume!(builder, h)
    regionpoint!(builder, L_width / 2, L_width / 2, L_low + L_mid / 2)

    # upper cell region
    cellregion!(builder, cell_up)
    maxvolume!(builder, H)
    regionpoint!(builder, L_width / 2, L_width / 2, L_low + L_mid + L_up / 2)

    return simplexgrid(builder)
end
