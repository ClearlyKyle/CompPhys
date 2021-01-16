import optparse
import calc_potential
import make_index
import sys

def parse_args():
    # Purely for testing, makes it faster to input args on PyCharm
    # ---------------------------------------------------------

    # NO INDEX FILE TEST                        - PASSED!
    # args = "-f freevolume.xtc -s freevolume.tpr"

    # STANDARD TEST                             - PASSED!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx"

    # Correct TEST                              - PASSED!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --correct"

    # Spherical TEST                            - FAILED!?
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --spherical"

    # Membrane direction TEST                   - PASSED!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx -d X"

    # Number of slices test TEST                - PASSED!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --sl 100"

    # Discard number of first slices TEST       - ERROR!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --cb 2"

    # Discard number of last slices TEST        - ERROR!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --ce 2"

    # Translate all coordinates TEST            - PASSED!
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --tz 1"

    # Translate all coordinates TEST
    # args = "-f freevolume.xtc -s freevolume.tpr -n index.ndx --ng 2"

    # TESTING A REAL MEMBRANE
    # args = "-f step7_11_0-50ns.xtc -s step7_11.tpr -n index.ndx --sl 100"

    # ---------------------------------------------------------

    usage = "Usage: %prog [options] arg1 arg2"
    version = "Poisson Solver v1.0"
    description = "[THISMODULE] computes the electrostatical potential " \
                  "across the box. The potential is calculated by first " \
                  "summing the charges per slice and then integrating twice " \
                  "of this charge distribution. Periodic boundaries are " \
                  "not taken into account. Reference of potential is " \
                  "taken to be the left side of the box. It is also possible " \
                  "to calculate the potential in spherical coordinates as " \
                  "function of r by calculating a charge distribution in " \
                  "spherical slices and twice integrating them. epsilon_r " \
                  "is taken as 1, but 2 is more appropriate in many cases."

    parser = optparse.OptionParser(usage=usage, version=version, description=description)

    '''
    parser.add_option("-f", "--file",
                      action="store",  # default is store, don't really need this
                      type="string",  # default is string, so not needed when needed input is string
                      dest="filename",  # variable name, to use it, use -> options.filename
                      default="FILE.dat",  # setting default value if no arg value is seen
                      help="Enter output filename",  # set help message
                      metavar="FILENAME")  # used in help messages in a place of an expected argument
    '''

    # ---------------------------------------------------------d------
    # File options
    group = optparse.OptionGroup(parser, "Options to specify input and output files")

    # Input Files
    group.add_option("-f", "--trajfile",  # Input File
                     action="store",
                     type="string",
                     dest="input_traj_filename",
                     metavar="INPUT_FILENAME",
                     help="Trajectory: xtc gro (Currently only .xtc tested)")

    group.add_option("-s", "--topfile",  # Input File
                     action="store",
                     type="string",
                     dest="input_top_filename",
                     metavar="INPUT_FILENAME",
                     help="Topology: tpr (Currently only .xtc tested)")

    group.add_option("-n", "--indexfile",  # Input File
                     action="store",
                     type="string",
                     dest="input_index_filename",
                     metavar="INPUT_FILENAME",
                     default=False,
                     help="Index: ndx file")

    # Output Files
    group.add_option("--op",  # Output Potential File
                     action="store",
                     type="string",
                     dest="output_potential_filename",
                     default="POTENTIAL_DATA.dat",
                     metavar="OUTPUT_POTENTIAL_FILENAME",
                     help="Output the potential (.xvg)")

    group.add_option("--oc",  # Output Charge File
                     action="store",
                     type="string",
                     dest="output_charge_filename",
                     default="CHARGE_DATA.dat",
                     metavar="OUTPUT_CHARGE_FILENAME",
                     help="Output the charge (.xvg)")

    group.add_option("--of",  # Output Field File
                     action="store",
                     type="string",
                     dest="output_field_filename",
                     default="FIELD_DATA.dat",
                     metavar="OUTPUT_FIELD_FILENAME",
                     help="Output the field (.xvg)")

    parser.add_option_group(group)

    # ---------------------------------------------------------------
    # Other options
    group = optparse.OptionGroup(parser, "Other Options")
    group.add_option("-t",  # Time calculation
                     action="store_true",
                     dest="time_calc",
                     default=False,
                     help="Calculate time for the calculation (s)")

    group.add_option("-d",  # Take the normal on the membrane in direction X, Y or Z.
                     action="store",
                     type="string",
                     dest="axtitle",
                     default='Z',
                     metavar="<string> (Z)",
                     help="Take the normal on the membrane in direction X, Y or Z.")

    group.add_option("--tz",  # Translate all coordinates by this distance in the direction of the box
                     action="store",
                     type="int",
                     dest="fudge_z",
                     default=0,
                     metavar="<real> (0)",
                     help="Translate all coordinates by this distance in the direction of the box")

    group.add_option("--b",  # First frame (ps) to read from trajectory
                     action="store",
                     type="int",
                     dest="first_frame",
                     default=0,
                     metavar="<time> (0)",
                     help="First frame (ps) to read from trajectory")

    group.add_option("--e",  # Last frame (ps) to read from trajectory
                     action="store",
                     type="int",
                     dest="last_frame",
                     default=0,
                     metavar="<time> (0)",
                     help="Last frame (ps) to read from trajectory")

    group.add_option("--sl",  # Dividing the box in this number of slices.
                     action="store",
                     type="int",
                     dest="slice_num",
                     default=10,
                     metavar="<int> (10)",
                     help="Calculate potential as function of boxlength, dividing the box in this number of slices.")

    group.add_option("--cb",  # Skip number of slices
                     action="store",
                     type="int",
                     dest="skip_init_val",
                     default=0,
                     metavar="<int> (0)",
                     help="Discard this number of first slices of box.")

    group.add_option("--ce",  # Skip end of box slices
                     action="store",
                     type="int",
                     dest="skip_end_val",
                     default=0,
                     metavar="<int> (0)",
                     help="Discard this number of last slices of box")

    group.add_option("--ng",  # Number of groups
                     action="store",
                     type="int",
                     dest="num_groups",
                     default=1,
                     metavar="<int> (1)",
                     help="Number of groups to consider")

    group.add_option("--spherical",  # Calculate spherical thingy
                     action="store_true",
                     dest="bSpherical",
                     default=False,
                     help="Calculate spherical thingy")

    group.add_option("--correct",  # Assume net zero charge of groups to improve accuracy
                     action="store_true",
                     dest="bCorrect",
                     default=False,
                     help="Assume net zero charge of groups to improve accuracy")

    parser.add_option_group(group)

    # ---------------------------------------------------------------
    # Debug options
    group = optparse.OptionGroup(parser, "Debug Options")
    group.add_option("--debug",
                     action="store_true",
                     help="Print debug information, currently not implemented")

    parser.add_option_group(group)

    # (options, args) = parser.parse_args(args.split())
    (options, args) = parser.parse_args()

    # ---------------------------------------------------------------
    if len(args) == 1:
        parser.error("incorrect number of arguments")

    # Make an index file if one is not given in command line
    if not options.input_index_filename:
        options.input_index_filename = make_index.make_index(options.input_traj_filename, options.input_top_filename)

    filenames = [options.input_traj_filename, options.input_top_filename, options.input_index_filename,
                 options.output_potential_filename, options.output_charge_filename, options.output_field_filename]

    other_options = [options.time_calc, options.axtitle, options.fudge_z, options.first_frame, options.last_frame,
                     options.slice_num, options.skip_init_val, options.skip_end_val, options.num_groups,
                     options.bSpherical, options.bCorrect]

    # Call to the main function for calculating the potential
    calc_potential.main(filenames, other_options)


# --------------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------------
if __name__ == "__main__":
    parse_args()
