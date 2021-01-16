import numpy
import MDAnalysis as mda
import matplotlib.pylab as plt

EPS0 = 8.85418781762e-12
ELC = 1.6021766208e-19
GMX_DOUBLE_MIN = 2.22507386e-308

ce = 0
cb = 0


# Should take in a index file, with number of groups you want
#   to look at, and return array of size NxM,
#       - where N is number of groups and M is group values
def parse_ndx(filename, num_groups):
    # Non dictionary method
    grp_names = []
    grp_values = []

    data = []

    with open(filename, "r") as f:
        for line in f:
            if line[0] == "[":
                grp_names.append(line.split(' ')[1])
                continue

            if len(line) > 1:
                data.extend(list(map(int, line.split())))

            elif line == '\n':
                if len(data) != 0:
                    grp_values.append(data)
                    data = []

        # This catches the last group at bottom of files, cheesy
        grp_values.append(data)

    i = 0
    print("\tSelect a group:")
    for key in grp_names:
        print("\t", i, ":", key)
        i += 1

    # loop over number of groups user wants (-nr command)
    group_selection = []
    for k in range(0, num_groups):

        while True:
            # usr_input = int(input("Choose group {} -> ".format(k)))
            usr_input = 0
            if usr_input < i:
                group_selection.append(grp_values[usr_input])
                break
            else:
                print("Error choosing group number")

        print("Group '{}' was selected".format(grp_names[usr_input]))

    return group_selection, grp_names[usr_input]


# this function integrates the array data and returns the resulting array */
# routine uses simple trapezoid rule
def p_integrate(data, ndata, slWidth):
    global ce, cb

    result = []

    if ndata <= 2:
        print("Warning: nr of slices very small. This will result in nonsense.")

    print("Integrating from slice {} to slice {}".format(cb, ndata - ce))

    for mySlice in range(cb, (ndata - ce)):
        total = 0
        for i in range(cb, mySlice):
            total += (slWidth * (data[i] + 0.5 * (data[i + 1] - data[i])))
        result.append(total)

    return result


# output_data = calc_potential(universe, gnx_in, fudge_z_in, ngrps,
#                             nslices, axis_in, bSpherical_in, bCorrect_in)
def calc_potential(u, gnx, fudge_z, num_groups, num_slices,
                   axis, bSpherical, bCorrect):
    ax1 = 0
    ax2 = 0
    num_frames = 0
    box = u.trajectory.ts.triclinic_dimensions / 10

    index = gnx

    if axis == 0:
        ax1 = 1
        ax2 = 2
    elif axis == 1:
        ax1 = 0
        ax2 = 2
    elif axis == 2:
        ax1 = 0
        ax2 = 1
    else:
        print("Invalid axes. Terminating")
        exit()

    natoms = len(u.atoms)

    if natoms == 0:
        print("Count not read coordinates from statusfile.")

    if not num_slices:
        num_slices = box[axis, axis] * 10

    print("Dividing the box in {} slices".format(num_slices))

    slField = numpy.zeros((num_groups, num_slices))
    slCharge = numpy.zeros((num_groups, num_slices))
    slPotential = numpy.zeros((num_groups, num_slices))

    '''
    SECTION ON TRAJECTORY ANALYSIS
    '''
    # ------ START PROCESSING TRAJECTORY ------
    print("Number of frames :", len(u.trajectory))

    for ts in range(0, 1):
        # for ts in range(0, len(u.trajectory) - 1):

        box = u.trajectory.ts.triclinic_dimensions / 10
        slWidth = box[axis, axis] / num_slices

        if ts % 500 == 0:
            print("Frame =", ts)

        # Position of all atoms
        x0 = u.trajectory.ts.positions / 10

        select_string = "bynum 0"

        for i in gnx[0]:
            select_string += " " + str(i)

        # Center of mass, only used in spherical mode
        xcm = -1 * u.select_atoms(select_string).center_of_mass() / 10

        for n in range(0, num_groups):
            '''
            Check whether we actually have all positions of the requested index group
            in the trajectory file
            '''
            if len(gnx[n]) > natoms:
                print("You selected a group with {} atoms, but only {} atoms were "
                      "found in the trajectory.".format(len(gnx[n]), natoms))
                exit(0)

            for i in range(0, len(gnx[n])):

                if bSpherical:

                    x0[index[n][i], :] += xcm

                    # only distance from origin counts, not sign
                    slice_val = int(numpy.linalg.norm(x0[index[n][i], :]) / slWidth)

                    '''
                    This is a nice check for spherical groups but not for all water
                    in a cubic box since a lot will fall outside the sphere
                        if slice > num_slices:
                            print("WARNING: slice = {}".format(slice_val))
                    '''
                    slCharge[n, slice_val] += u.atoms.charges[index[n][i]]

                else:
                    z = x0[index[n][i]][axis]
                    z += fudge_z

                    if z < 0:
                        z += box[axis][axis]
                    if z > box[axis][axis]:
                        z -= box[axis][axis]

                    slice_val = int(z / slWidth)

                    slCharge[n, slice_val] += u.atoms.charges[index[n][i]]

        num_frames += 1  # Counter to display number of frames processed
        
        if len(u.trajectory) > 1:
        	u.trajectory.next()  # Gets next frame of trajectories

    # INTEGRATION STAGE --------------------------------------------------------

    # slCharge now contains the total charge per slice, summed over all frames.
    #   Now divide by num_frames and integrate twice

    if bSpherical:
        print("Read {} frames from trajectory. Calculating potential in "
              "spherical coordinates".format(num_frames))
    else:
        print("Read {} frames from trajectory. Calculating "
              "potential.".format(num_frames))

    for n in range(0, num_groups):
        for i in range(0, num_slices):

            # Spherical mode
            if bSpherical:
                '''
                charge per volume id now the summed charge, divided by
                the nr of frames and by the volume of the slice its in.
                4pi r^2 dr
                '''
                slVolume = 4 * numpy.pi * (i ** 2) * (slWidth ** 2) * slWidth

                if slVolume == 0:
                    slCharge[n, i] = 0
                else:
                    slCharge[n, i] /= (num_frames * slVolume)

            # Non-spherical mode
            else:
                # Get charge per volume
                slCharge[n, i] = slCharge[n, i] * num_slices / (
                    num_frames * box[axis][axis] * box[ax1][ax1] * box[ax2][ax2])

    # Now we have charge densities

    if bCorrect and not bSpherical:
        for n in range(0, num_groups):
            nn = 0
            qsum = 0

            for i in range(0, num_slices):
                if numpy.abs(slCharge[n, i]) >= GMX_DOUBLE_MIN:
                    nn += 1
                    qsum += slCharge[n, i]

            qsum /= nn

            for i in range(0, num_slices):
                if numpy.abs(slCharge[n, i]) >= GMX_DOUBLE_MIN:
                    slCharge[n, i] -= qsum

    for n in range(0, num_groups):
        slField[n, :] = p_integrate(slCharge[n, :], num_slices, slWidth)

    if bCorrect and not bSpherical:
        for n in range(0, num_groups):
            nn = 0
            qsum = 0

            for i in range(0, num_slices):
                if numpy.abs(slCharge[n, i]) >= GMX_DOUBLE_MIN:
                    nn += 1
                    qsum += slField[n, i]

            qsum /= nn

            for i in range(0, num_slices):
                if numpy.abs(slCharge[n, i]) >= GMX_DOUBLE_MIN:
                    slField[n, i] -= qsum

    for n in range(0, num_groups):
        slPotential[n, :] = p_integrate(slField[n, :], num_slices, slWidth)

    # Now correct for EPS0, and in spherical case for r
    for n in range(0, num_groups):
        for i in range(0, num_slices):

            if bSpherical:
                slPotential[n, i] = ELC * slPotential[n, i] * -1.0E9 / (EPS0 * i * slWidth)
                slField[n, i] = ELC * slField[n, i] * 1E18 / (EPS0 * i * slWidth)
            else:
                slPotential[n, i] = ELC * slPotential[n, i] * -1.0E9 / EPS0
                slField[n, i] = ELC * slField[n, i] * 1E18 / EPS0

    return [numpy.array([slPotential, slCharge, slField]), slWidth, num_slices]


def plot_potential(data, group_name, slWidth, nslices):
    global cb, ce

    x_axis = [(s * slWidth) for s in range(ce, (nslices - ce))]
    # for s in range(ce, (nslices - ce)):
    #     x_axis.append((s * slWidth))

    titles = ['Electrostatic Potential', 'Charge Distribution', 'Electric Field']
    ylabels = ['Potential (V)', 'Charge density (q/nm^3', 'Field (V/nm)']
    xlabels = ['Box (nm)', 'Box (nm)', 'Box (nm)']

    # data[0, :] *= 10
    data[2, :] /= 1E9

    with open("potential.dat", "w") as f:
        f.write("x\t\ty\n")
        for i in range(0, (len(x_axis))):
            f.write("{}\t\t{}\n".format(x_axis[i], data[0, 0, i]))

    with open("charge.dat", "w") as f:
        f.write("x\t\ty\n")
        for i in range(0, (len(x_axis))):
            f.write("{}\t\t{}\n".format(x_axis[i], data[1, 0, i]))

    with open("field.dat", "w") as f:
        f.write("x\t\ty\n")
        for i in range(0, (len(x_axis))):
            f.write("{}\t\t{}\n".format(x_axis[i], data[2, 0, i]))

    for i in range(0, len(titles)):
        plt.figure(titles[i])
        plt.plot(x_axis, data[i, 0], label=group_name)
        plt.title(titles[i])
        plt.xlabel(xlabels[i])
        plt.ylabel(ylabels[i])
        plt.legend(loc="best")


# ---------------------------------------------------------------------------------
def main(file_names, other_options):
    global cb, ce

    axtitle = other_options[1]
    axis = ord(axtitle) - ord('X')

    fudge_z_in = other_options[2]

    nslices = other_options[5]

    cb = other_options[6]
    ce = other_options[7]

    ngrps = other_options[8]

    bSpherical_in = other_options[9]
    bCorrect_in = other_options[10]

    xtc_filename = file_names[0]
    tpr_filename = file_names[1]
    ndx_filename = file_names[2]

    universe = mda.Universe(tpr_filename, xtc_filename)

    gnx_in, grp_name = parse_ndx(ndx_filename, ngrps)

    output_data = calc_potential(universe, gnx_in, fudge_z_in, ngrps,
                                 nslices, axis, bSpherical_in, bCorrect_in)

    plot_potential(output_data[0], grp_name, output_data[1], output_data[2])

    plt.show()
