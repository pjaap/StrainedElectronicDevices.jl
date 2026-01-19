# from 10.48550/arXiv.2312.09235
function create_cuboid_grid_with_wire(;
        wire_width = 80,
        x_length = 700,
        y_length = 700,
        virtual_substrate = 2000,
        quantum_well = 9,
        spacer = 30,
        insulator = 10,
        wire_height = 60,
        cover = 2,
        front = 1,
        back = 2,
        left = 3,
        right = 4,
        bottom = 5,
        top = 6,
        internal = 7,
        H = 1.0e9,
        h = 1.0e9
    )

    builder = SimplexGridBuilder(; Generator = TetGen)

    # virtual virtual_substrate
    low = 1.5 - quantum_well - virtual_substrate
    hi = low + virtual_substrate

    p11 = point!(builder, -x_length / 2, -y_length / 2, low)
    p12 = point!(builder, +x_length / 2, -y_length / 2, low)
    p13 = point!(builder, -x_length / 2, +y_length / 2, low)
    p14 = point!(builder, +x_length / 2, +y_length / 2, low)

    facetregion!(builder, bottom)
    facet!(builder, p11, p12, p14, p13)

    p21 = point!(builder, -x_length / 2, -y_length / 2, hi)
    p22 = point!(builder, +x_length / 2, -y_length / 2, hi)
    p23 = point!(builder, -x_length / 2, +y_length / 2, hi)
    p24 = point!(builder, +x_length / 2, +y_length / 2, hi)

    facetregion!(builder, front)
    facet!(builder, p11, p12, p22, p21)
    facetregion!(builder, back)
    facet!(builder, p13, p14, p24, p23)
    facetregion!(builder, left)
    facet!(builder, p11, p21, p23, p13)
    facetregion!(builder, right)
    facet!(builder, p12, p22, p24, p14)

    # QW
    low = 1.5 - quantum_well
    hi = low + quantum_well

    p11 = point!(builder, -x_length / 2, -y_length / 2, low)
    p12 = point!(builder, +x_length / 2, -y_length / 2, low)
    p13 = point!(builder, -x_length / 2, +y_length / 2, low)
    p14 = point!(builder, +x_length / 2, +y_length / 2, low)

    facetregion!(builder, internal)
    facet!(builder, p11, p12, p14, p13)

    p21 = point!(builder, -x_length / 2, -y_length / 2, hi)
    p22 = point!(builder, +x_length / 2, -y_length / 2, hi)
    p23 = point!(builder, -x_length / 2, +y_length / 2, hi)
    p24 = point!(builder, +x_length / 2, +y_length / 2, hi)

    facetregion!(builder, front)
    facet!(builder, p11, p12, p22, p21)
    facetregion!(builder, back)
    facet!(builder, p13, p14, p24, p23)
    facetregion!(builder, left)
    facet!(builder, p11, p21, p23, p13)
    facetregion!(builder, right)
    facet!(builder, p12, p22, p24, p14)

    # spacer
    low = 1.5
    hi = low + spacer

    p11 = point!(builder, -x_length / 2, -y_length / 2, low)
    p12 = point!(builder, +x_length / 2, -y_length / 2, low)
    p13 = point!(builder, -x_length / 2, +y_length / 2, low)
    p14 = point!(builder, +x_length / 2, +y_length / 2, low)

    facetregion!(builder, internal)
    facet!(builder, p11, p12, p14, p13)

    p21 = point!(builder, -x_length / 2, -y_length / 2, hi)
    p22 = point!(builder, +x_length / 2, -y_length / 2, hi)
    p23 = point!(builder, -x_length / 2, +y_length / 2, hi)
    p24 = point!(builder, +x_length / 2, +y_length / 2, hi)

    facetregion!(builder, front)
    facet!(builder, p11, p12, p22, p21)
    facetregion!(builder, back)
    facet!(builder, p13, p14, p24, p23)
    facetregion!(builder, left)
    facet!(builder, p11, p21, p23, p13)
    facetregion!(builder, right)
    facet!(builder, p12, p22, p24, p14)


    # insulator
    low = 1.5 + spacer
    hi = low + insulator
    mid = hi - cover

    p11 = point!(builder, -x_length / 2, -y_length / 2, low)
    p12 = point!(builder, +x_length / 2, -y_length / 2, low)
    p13 = point!(builder, -x_length / 2, +y_length / 2, low)
    p14 = point!(builder, +x_length / 2, +y_length / 2, low)

    facetregion!(builder, internal)
    facet!(builder, p11, p12, p14, p13)

    p21 = point!(builder, -x_length / 2, -y_length / 2, hi)
    p22 = point!(builder, -wire_width / 2 - cover, -y_length / 2, hi)
    p23 = point!(builder, -wire_width / 2, -y_length / 2, mid)
    p24 = point!(builder, +wire_width / 2, -y_length / 2, mid)
    p25 = point!(builder, +wire_width / 2 + cover, -y_length / 2, hi)
    p26 = point!(builder, +x_length / 2, -y_length / 2, hi)

    p31 = point!(builder, -x_length / 2, +y_length / 2, hi)
    p32 = point!(builder, -wire_width / 2 - cover, +y_length / 2, hi)
    p33 = point!(builder, -wire_width / 2, +y_length / 2, mid)
    p34 = point!(builder, +wire_width / 2, +y_length / 2, mid)
    p35 = point!(builder, +wire_width / 2 + cover, +y_length / 2, hi)
    p36 = point!(builder, +x_length / 2, +y_length / 2, hi)

    facetregion!(builder, front)
    facet!(builder, [p11, p12, p26, p25, p24, p23, p22, p21])
    facetregion!(builder, back)
    facet!(builder, [p13, p14, p36, p35, p34, p33, p32, p31])
    facetregion!(builder, left)
    facet!(builder, p11, p21, p31, p13)
    facetregion!(builder, right)
    facet!(builder, p12, p26, p36, p14)

    facetregion!(builder, top)
    facet!(builder, p21, p22, p32, p31)
    facet!(builder, p25, p26, p36, p35)

    facetregion!(builder, internal)
    facet!(builder, p22, p23, p33, p32)
    facetregion!(builder, internal)
    facet!(builder, p23, p24, p34, p33)
    facetregion!(builder, internal)
    facet!(builder, p24, p25, p35, p34)

    p41 = point!(builder, -wire_width / 2 - cover, -y_length / 2, hi + wire_height)
    p42 = point!(builder, +wire_width / 2 + cover, -y_length / 2, hi + wire_height)
    p43 = point!(builder, -wire_width / 2 - cover, +y_length / 2, hi + wire_height)
    p44 = point!(builder, +wire_width / 2 + cover, +y_length / 2, hi + wire_height)

    p51 = point!(builder, -wire_width / 2, -y_length / 2, hi + wire_height - cover)
    p52 = point!(builder, +wire_width / 2, -y_length / 2, hi + wire_height - cover)
    p53 = point!(builder, -wire_width / 2, +y_length / 2, hi + wire_height - cover)
    p54 = point!(builder, +wire_width / 2, +y_length / 2, hi + wire_height - cover)

    facetregion!(builder, front)
    facet!(builder, [p22, p23, p51, p52, p24, p25, p42, p41])
    facet!(builder, [p23, p24, p52, p51])

    facetregion!(builder, back)
    facet!(builder, [p32, p33, p53, p54, p34, p35, p44, p43])
    facet!(builder, [p33, p34, p54, p53])


    facetregion!(builder, top)
    facet!(builder, p22, p41, p43, p32)
    facet!(builder, p41, p42, p44, p43)
    facet!(builder, p42, p25, p35, p44)

    facetregion!(builder, internal)
    facet!(builder, p23, p51, p53, p33)
    facet!(builder, p51, p52, p54, p53)
    facet!(builder, p24, p52, p54, p34)


    # virtual_substrate cell region
    cellregion!(builder, 1)
    maxvolume!(builder, Inf64)
    regionpoint!(builder, 0, 0, 1.5 - quantum_well - virtual_substrate / 2)

    # quantum well cell region
    cellregion!(builder, 2)
    maxvolume!(builder, h)
    regionpoint!(builder, 0, 0, 0)

    # spacer cell region
    cellregion!(builder, 3)
    maxvolume!(builder, H)
    regionpoint!(builder, 0, 0, 1.5 + spacer / 2)


    # insulator cell region
    cellregion!(builder, 4)
    maxvolume!(builder, H)
    regionpoint!(builder, 0, 0, 1.5 + spacer + (insulator - cover) / 2)
    regionpoint!(builder, 0, 0, 1.5 + spacer + insulator + wire_height - cover / 2)


    # # wire cell region
    cellregion!(builder, 5)
    maxvolume!(builder, H)
    regionpoint!(builder, 0, 0, 1.5 + spacer + insulator - cover + wire_height / 2)

    return simplexgrid(builder)
end
