'''
Brag Reflector

Kyle McIntosh
'''
import pylab
import numpy
import sys


# ----------------------------------------------------------
# Calculate the reflectivity from given index, values
def calc_reflection(index, N, wave_info, wave_change):
    # Index variable holds all the refractive index values we need
    #   I change them to better variable names to help with readability and
    #   helps when implementing the function
    n_1 = index[0]
    n_H = index[1]
    n_L = index[2]
    n_S = index[3]

    # Get the angle from the normal and wavelength of the light
    in_angle = wave_info[0]
    wavelength = wave_info[1]

    theta_1 = in_angle * numpy.pi / 180
    theta_2 = numpy.arcsin((n_1 / n_H) * numpy.sin(theta_1))
    theta_3 = numpy.arcsin((n_H / n_L) * numpy.sin(theta_2))

    h_2 = wavelength / (4 * n_H)
    h_3 = wavelength / (4 * n_L)

    p_2 = n_H * numpy.cos(theta_2)
    p_3 = n_L * numpy.cos(theta_3)

    # create two empty arrays to gold matrix data
    #   also note have to specify its complex data type
    M1 = numpy.zeros((2, 2), dtype=numpy.complex)
    M2 = numpy.zeros((2, 2), dtype=numpy.complex)

    # find angle of light into the substrate material
    theta_S = numpy.arcsin((n_L / n_S) * numpy.sin(theta_1))

    p_S = n_S * numpy.cos(theta_S)
    p_1 = n_1 * numpy.cos(theta_1)

    R = []

    for i in wave_change:
        # re-calculate beta2 and beta3 for each wavelength in the range
        beta_2 = ((2 * numpy.pi) / i) * n_H * h_2 * numpy.cos(theta_2)
        beta_3 = ((2 * numpy.pi) / i) * n_L * h_3 * numpy.cos(theta_3)

        # values for Matrix 1
        M1[0, 0] = numpy.cos(beta_2)
        M1[0, 1] = (-1j / p_2) * numpy.sin(beta_2)
        M1[1, 0] = (-1j * p_2) * numpy.sin(beta_2)
        M1[1, 1] = numpy.cos(beta_2)

        # values for Matrix 1
        M2[0, 0] = numpy.cos(beta_3)
        M2[0, 1] = (-1j / p_3) * numpy.sin(beta_3)
        M2[1, 0] = (-1j * p_3) * numpy.sin(beta_3)
        M2[1, 1] = numpy.cos(beta_3)

        # The way numpy handles matrix multiplication you can just use M1 * M2 for complex
        #   numbers, using dot(M1, M2) is equal to M1 * M2
        #   linalg.matrix.power(M, N) multiplies a matrix by its self N times, so this N is
        #   the number of layers of our system
        M = numpy.linalg.matrix_power((numpy.dot(M1, M2)), N)

        little_r = (((M[0, 0] + (M[0, 1] * p_S)) * p_1) - (M[1, 0] + (M[1, 1] * p_S))) / (((M[0, 0] + (M[0, 1] * p_S)) * p_1) + (M[1, 0] + (M[1, 1] * p_S)))
        R.append(numpy.abs(little_r ** 2))

    return R


# ----------------------------------------------------------
# Plot the Reflectivity versus the wavelength of visible light
def plot_graph(x, y):
    # Takes the x axis as N and y axis as reflectivity
    pylab.plot(x, y, label='N = {}'.format(N))
    pylab.title("Reflectivity vs Wavelength for N layered Bragg Reflector")
    pylab.xlabel("Wavelength (m)")
    pylab.ylabel("Reflectivity")
    pylab.legend()
    pylab.show()


# ----------------------------------------------------------
# Get user input value for the index of materials
def get_index_value(name):
    while True:
        try:
            # Ask for user to input a name for the parameter required
            user_input = input("Enter value for {} -> ".format(name))

            # Return the user input as a float number
            return numpy.abs(float(user_input))

        # Should catch any common errors inputting
        except TypeError:
            print("Invalid Input, check the value inputted.")


# ----------------------------------------------------------
# Change one of the parameters that have previously been set
def change_value():
    global N  # breaks if this isn't here, not sure why "index" or "wave_data" don't need this

    print("Example of changing value 'n_1 = 1'")

    while True:
        while True:
            user_input = input("Enter what you would like to change ->")

            # Remove the spaces and split the word on the '=' sign
            #   makes manipulation later easier
            splitted = user_input.replace(" ", "").split("=")

            if splitted[0] == "n_1":
                index[0] = float(splitted[-1])
            elif splitted[0] in ["n_H", "n_h"]:
                index[1] = float(splitted[-1])
            elif splitted[0] in ["n_L", "n_l"]:
                index[2] = float(splitted[-1])
            elif splitted[0] in ["n_S", "n_s"]:
                index[3] = float(splitted[-1])
            elif splitted[0] in ["N", "n"]:
                N = int(splitted[-1])
            elif splitted[0].lower() in ["wavelength", "wave"]:
                wave_data[1] = float(splitted[-1])
            elif splitted[0].lower() == "angle":
                wave_data[0] = float(splitted[-1])
            else:
                print("Invalid input, can only use key words:")
                print("\t\tn_1, n_H, n_L, n_S, N, wavelength, angle")
                continue

            break

        exit_loop = input("Change another value y/n? ->")
        if exit_loop.lower() == "n":
            break


# ----------------------------------------------------------
# Checks if data is already inputted, so user cant plot no data
def check_data_loaded(index):
    if index == [0, 0, 0, 0]:
        print("Error: No data values are currently inputted.")
        return 1
    else:
        return 0


# ----------------------------------------------------------
# Opening Instruction screen only seen on start-up
def display_start():
    print("\nBRAGG REFLECTOR 	by Kyle McIntosh"
          "\n-----------------------------------------------------"
          "\nA program to calculate the reflectivity of alternating"
          "\nhigh and low refractive index materials."
          "\nUsing this arrangement it can be used to produce"
          "\nvery high reflective mirrors."
          "\n"
          "\nn_1 = Refractive index above the layers"
          "\nn_H = High refractive index material"
          "\nn_L = Lower refractive index"
          "\nn_S = Subtrate refractive index"
          "\n"
          "\nTo input powers 'e' must be used in place of 'x10'"
          # "\n   - or: m = mili, u = micro, n = nano, p = pico"
          "\n")


# ----------------------------------------------------------
# Shows the aviable inputs on screen
def display_help():
    print("\nOptions ---------------------------------------------"
          "\n\tplot\t- Plot the reflectivity vs N"
          "\n\tinput\t- Input data values"
          "\n\tsample\t- Imports data from a sample file 'sample_data.txt'"
          "\n\tchange\t- Change a specific value for a parameter"
          "\n\tdisplay\t- Show a little diagram of the current setup"
          "\n\tprange\t- Enter multiple N values to plot"
          "\n\t"
          "\n\tquit - exits the program"
          "\n=====================================================")


# ----------------------------------------------------------
# Displays the current parameters inputted into the system
#   as a little diagram
def display_layers(refractive_indexs, N, wave_info):
    # Displays the layers with refractive index's
    print("                         |   {}m @ {}degrees".format(wave_info[1], wave_info[0]))
    print("\t n_1 = {0:.2f}         \|/".format(refractive_indexs[0]))
    print("\t---------------------------------")
    print("\t n_H = {0:.2f}".format(refractive_indexs[1]))
    print("\t---------------------------------")
    print("\t n_L = {0:.2f}".format(refractive_indexs[2]))
    print("\t---------------------------------")
    print("\t n_S = {0:.2f}".format(refractive_indexs[3]))
    print("")


# ----------------------------------------------------------
# array to store the index values of the layers
#   [n_1, n_H, n_L, n_S]
index = [0, 0, 0, 0]

# Number of layers to the system
N = 0

# Data for the incoming light
#   [angle, wavelength]
wave_data = [0, 0]

# ----------------------------------------------------------
# ----------------------------------------------------------
# Load initial function on program start
display_start()  # Display the start screen on startup
display_help()  # Displays the help screen on startup

# Main program loop
while True:
    try:
        # ----------------------------------------------------------
        # Input for the users console command
        menu = input("--->").lower()

        # ----------------------------------------------------------
        # Plot current data inputted
        if menu in ["plot", "p"]:
            if check_data_loaded(index) == 0:
                print("Plotting, Reflectivity vs Wavelength")

                # Range which to plot over
                temp = numpy.arange(300, 700) * 10 ** -9
                R = calc_reflection(index, N, wave_data, temp)

                plot_graph(temp, R)

        # ----------------------------------------------------------
        # Plot a range of N values
        elif menu in ["prange", "pr"]:
            if check_data_loaded(index) == 0:
                temp = numpy.arange(300, 700) * 10 ** -9

                N_values = [int(x) for x in input("Enter range of N values, separated by ',' ->").split(',')]

                for i in N_values:
                    R = calc_reflection(index, i, wave_data, temp)
                    pylab.plot(temp, R, label='N = {}'.format(i))

                pylab.title("Reflectivity vs Wavelength for N layered Bragg Reflector")
                pylab.xlabel("Wavelength (m)")
                pylab.ylabel("Reflectivity")
                pylab.legend()
                pylab.show()

        # ----------------------------------------------------------
        # Get user input for the values needed
        elif menu in ["input", "i"]:
            index[0] = get_index_value("n_1")
            index[1] = get_index_value("n_H")
            index[2] = get_index_value("n_L")
            index[3] = get_index_value("n_S")

            # Validate that N layers has to be an integer
            while True:
                try:
                    N = int(input("Enter value for N \t->"))
                except ValueError:
                    print("Has to be Integer value.")
                    continue
                break

            # Allows the use to use p for pico meters.. ext
            wave_input = input("Enter value for wavelength ->")
            if wave_input[-1].lower() == "m":
                wave_data[1] = float(wave_input[:-1]) * 10 ** -3
            elif wave_input[-1].lower() == "u":
                wave_data[1] = float(wave_input[:-1]) * 10 ** -6
            elif wave_input[-1].lower() == "n":
                wave_data[1] = float(wave_input[:-1]) * 10 ** -9
            elif wave_input[-1].lower() == "p":
                wave_data[1] = float(wave_input[:-1]) * 10 ** -12
            else:
                wave_data[1] = float(wave_input)

            # Get the angle from normal
            #   Can input the angle in degree "deg" or radians"
            wave_input = input("Enter angle from normal \t->")

            if wave_input[-3:].lower() == "deg":  # Handle degrees
                wave_data[0] = float(wave_input[:-3])

            elif wave_input[-3:].lower() == "rad":  # Handle radians
                wave_data[0] = float(wave_input[:-3]) * 180 / numpy.pi

            # Can alternatively input the angle like 'pi/2' or '3pi'
            #   this may be prone to bugs, tested as much as I could
            elif "pi" in wave_input:
                if "/" in wave_input:
                    temp = wave_input.split('/')
                    if temp[0].split('pi')[0] == '':
                        wave_data[0] = numpy.pi / float(temp[1]) * 180 / numpy.pi
                    else:
                        wave_data[0] = float(temp[0].split('pi')[0]) * numpy.pi / float(temp[1]) * 180 / numpy.pi
                else:
                    if wave_input.split('pi')[0] == '':
                        wave_data[0] = 180
                    else:
                        wave_data[0] = float(wave_input[0].split('pi')[0]) * numpy.pi * 180 / numpy.pi
            else:  # just assume degrees
                wave_data[0] = float(wave_input)

        # ----------------------------------------------------------
        # Run sample data from a file
        elif menu in ["sample", "s"]:
            file_data = numpy.loadtxt('sample_data.txt', comments='#')

            # Import first 4 values of sample file as the index values
            for i in range(0, 4):
                index[i] = file_data[i]

            wave_data[0] = file_data[4]  # angle from normal
            wave_data[1] = file_data[5]  # wavelength of light in

            N = int(file_data[6])  # number of layers

            # display the diagram of the system
            display_layers(index, N, wave_data)

        # ----------------------------------------------------------
        # Change one of the values stored
        elif menu in ["change"]:
            change_value()

        # ----------------------------------------------------------
        # Show current data inputted in a diagram
        elif menu in ["display", "d", "disp"]:
            # If there is no data currently in the system then it will
            #   not show the diagram
            if index == [0, 0, 0, 0]:
                print("No index data inputted currently.")
            else:
                display_layers(index, N, wave_data)

        # ----------------------------------------------------------
        # Exit the program
        elif menu in ["exit", "end", "stop", "quit", "q"]:
            print("Exiting program")
            sys.exit()

        # ----------------------------------------------------------
        # No input or invalid input used
        else:
            print("No option {}.\nType 'help' for list of commands.".format(menu))

    # ----------------------------------------------------------
    # Fail safe for NameError
    except NameError:
        print("No data currently inputted.")
