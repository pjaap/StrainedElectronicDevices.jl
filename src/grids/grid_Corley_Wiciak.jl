"""
    from 10.1103/PhysRevApplied.20.024056

    lifted to 3D

    all quantities in [nm] !!
"""
function grid_Corley_Wiciak_3D(;
        region_TiN = 1,
        region_VS_fine = 2,
        region_VS_coarse = 3,
        bregion_left = 1,
        bregion_right = 2,
        bregion_bottom = 3,
        bregion_front = 4,
        bregion_back = 5,
        bregion_default = 6,
        electrode_width = 160.0,
        pitch_width = 140.0,
        electrode_height = 30.0,
        VS_fine_height = 200.0,
        VS_coarse_height = 400.0, # 4000 in the paper, reduced here
        depth = 50,
        h = 1,
        kwargs...
    )
    builder = SimplexGridBuilder(; Generator = TetGen)

    # accumulated thickness of all layers
    level = zeros(4)
    level[2] = VS_coarse_height
    level[3] = level[2] + VS_fine_height
    level[4] = level[3] + electrode_height

    # total width (pitch_width at both ends)
    width = electrode_width + pitch_width

    #  z      9--10      <- level 4
    #  ^      |   |
    #  |  5---6---7---8  <- level 3
    #  |  |           |
    #  |  3-----------4  <- level 2
    #  |  |           |
    #  |  1-----------2  <- level 1
    #  |
    #  +----------------> x

    # create point array as above with given y-value
    make_points(y) = [
        point!(builder, -width / 2, y, level[1]),
        point!(builder, width / 2, y, level[1]),

        point!(builder, -width / 2, y, level[2]),
        point!(builder, width / 2, y, level[2]),

        point!(builder, -width / 2, y, level[3]),
        point!(builder, -electrode_width / 2, y, level[3]),
        point!(builder, electrode_width / 2, y, level[3]),
        point!(builder, width / 2, y, level[3]),

        point!(builder, -electrode_width / 2, y, level[4]),
        point!(builder, electrode_width / 2, y, level[4]),
    ]

    points_front = make_points(-depth / 2)
    points_back = make_points(depth / 2)

    # we split the lower block in a half for at least two grid cells along y
    points_center = [
        point!(builder, -width / 2, 0.0, level[1]),
        point!(builder, width / 2, 0.0, level[1]),
        point!(builder, -width / 2, 0.0, level[2]),
        point!(builder, width / 2, 0.0, level[2]),
    ]

    facetregion!(builder, bregion_left)
    begin
        facet!(builder, [points_front[1], points_front[3], points_center[3], points_center[1]])
        facet!(builder, [points_center[1], points_center[3], points_back[3], points_back[1]])
        facet!(builder, [points_front[3], points_front[5], points_back[5], points_back[3], points_center[3]])
    end

    facetregion!(builder, bregion_right)
    begin
        facet!(builder, [points_front[2], points_front[4], points_center[4], points_center[2]])
        facet!(builder, [points_center[2], points_center[4], points_back[4], points_back[2]])
        facet!(builder, [points_front[4], points_front[8], points_back[8], points_back[4], points_center[4]])
    end

    facetregion!(builder, bregion_bottom)
    facet!(builder, [points_front[1], points_front[2], points_center[2], points_center[1]])
    facet!(builder, [points_center[1], points_center[2], points_back[2], points_back[1]])

    facetregion!(builder, bregion_front)
    begin
        facet!(builder, [points_front[1], points_front[2], points_front[4], points_front[3]])
        facet!(builder, [points_front[3], points_front[4], points_front[8], points_front[7], points_front[6], points_front[5]])
        facet!(builder, [points_front[6], points_front[7], points_front[10], points_front[9]])
    end

    facetregion!(builder, bregion_back)
    begin
        facet!(builder, [points_back[1], points_back[2], points_back[4], points_back[3]])
        facet!(builder, [points_back[3], points_back[4], points_back[8], points_back[7], points_back[6], points_back[5]])
        facet!(builder, [points_back[6], points_back[7], points_back[10], points_back[9]])
    end
    facetregion!(builder, bregion_default)
    begin
        # horizontal facets
        facet!(builder, [points_front[3], points_front[4], points_center[4], points_center[3]])
        facet!(builder, [points_center[3], points_center[4], points_back[4], points_back[3]])
        facet!(builder, [points_front[5], points_front[6], points_back[6], points_back[5]])
        facet!(builder, [points_front[6], points_front[7], points_back[7], points_back[6]])
        facet!(builder, [points_front[7], points_front[8], points_back[8], points_back[7]])
        facet!(builder, [points_front[9], points_front[10], points_back[10], points_back[9]])
        # vertical facets
        facet!(builder, [points_center[1], points_center[2], points_center[4], points_center[3]])
        facet!(builder, [points_front[6], points_front[9], points_back[9], points_back[6]])
        facet!(builder, [points_front[7], points_front[10], points_back[10], points_back[7]])
    end

    # mark cell regions
    cellregion!(builder, region_TiN)
    maxvolume!(builder, h)
    regionpoint!(builder, 0, 0, (level[3] + level[4]) / 2)

    cellregion!(builder, region_VS_fine)
    maxvolume!(builder, h)
    regionpoint!(builder, 0, 0, (level[2] + level[3]) / 2)

    cellregion!(builder, region_VS_coarse)
    maxvolume!(builder, Inf32)
    regionpoint!(builder, 0, depth / 4, (level[1] + level[2]) / 2)
    regionpoint!(builder, 0, -depth / 4, (level[1] + level[2]) / 2)

    return simplexgrid(builder)
end
