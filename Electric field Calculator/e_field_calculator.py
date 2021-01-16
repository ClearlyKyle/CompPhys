import numpy, tkinter, sys
from tkinter import filedialog

# Hides the default window created by tk
root = tkinter.Tk().withdraw()


def open_file():
    # Allows the option of csv or txt file input
    file_path = filedialog.askopenfilename(title="Select file to open.", initialdir=dir, filetypes=[("txt", ".txt"), ("csv", ".csv"), ("All files", "*.*")], defaultextension=".txt")

    # Changes how numpy handles the input based on whether csv or txt file
    if '.csv' in file_path:
        return numpy.loadtxt('{}'.format(file_path), delimiter=',')
    elif '.txt' in file_path:
        return numpy.loadtxt('{}'.format(file_path))
    else:
        print("No file selected")


def save_file(datatosave):
    # Gets the user to specify a location to save a file to
    file_path = filedialog.asksaveasfilename(title="Save the data as...", initialdir=dir, initialfile='New Text Document.txt', filetypes=[("txt", ".txt"), ("csv", ".csv"), ("All files", "*.*")], defaultextension=".txt")

    # Check for no location selected
    if file_path is None:
        return "File not saved, no save destination selected."

    numpy.savetxt('{}'.format(file_path), datatosave)


# Displays the data currently loaded, formatted on screen
def displaydata(d):
    print("=" * 30, "\n {:6}\t{:6}\t{:4}\t{}".format('x', 'y', 'z', 'Charge q'))
    print('=' * 30)
    print('\n'.join('     '.join(str(cell) for cell in row) for row in d))


# Electric field calculation
def electricfield(data):
    # Asks for points for the electric field
    point = input("Points for field:").split(',')

    # Calculate the distance between the points
    x_d = float(point[0]) - data[:, 0]
    y_d = float(point[1]) - data[:, 1]
    z_d = float(point[2]) - data[:, 2]

    distance = numpy.sqrt(x_d ** 2 + y_d ** 2 + z_d ** 2)

    # Unit vector calculation
    unit_vector = numpy.matrix([x_d / distance, y_d / distance, z_d / distance])

    constant_k = 8.9875517873681764 * 10 ** 9

    # Calculates the field in each direction
    field = numpy.array(((constant_k * data[:, 3]) / distance ** 2) * unit_vector.T)

    # Prints the electric field
    print("Field value: {:.4} N/C".format(numpy.sum([field])))

    # Print the direction of the field
    print("Direction: ", end='')
    try:
        if numpy.sum([field[:, 0]]) != 0: print("({:.4g}N/C)i ".format(numpy.sum([field[:, 0]])), end='')
        if numpy.sum([field[:, 1]]) != 0: print("({:.4g}N/C)j ".format(numpy.sum([field[:, 1]])), end='')
        if numpy.sum([field[:, 2]]) != 0: print("({:.4g}N/C)k ".format(numpy.sum([field[:, 2]])), end='')
        print()

    except IndexError:
        print("No direction")

    return point


# User can manually enter the data
def manualinput():
    # Input validation for number of charges
    while True:
        try:
            num_charges = int(input("How many Charges are needed?"))
            break
        except ValueError:
            print("Invalid input!")

    data = numpy.zeros([num_charges, 4])  # Empty matrix

    # Input validation for coordinates of each charge
    while True:
        try:
            for i in range(num_charges):
                data[i:, ] = input("Enter the X, Y and Z coordinates for charge q{},\nfollowed by the charge value --:".format(i + 1)).split(',')
            break
        except ValueError:
            print("Invalid input!")

    return data


# Shows the aviable inputs on screem
def display_help():
    print("\nOptions ---------------------------------------------"
          "\nsave\t- Save the file to a specific location"
          "\nopen\t- Open a file with data values"
          "\n\t\t\t(.csv or .txt file only)"
          "\ninput\t- Manualy input coordinates and charge values"
          "\nfield\t- Enter coordinates to find the electric field"
          "\n\t\t\tat a specific point."
          "\nadd\t\t- Add new coordinates and charge values to"
          "\n\t\t\tdata already inputted"
          "\ndisp\t- Display current data"
          "\nexit\t- Ends the program"
          "\n=====================================================")


# Opening Instruction screen only seen on start-up
print("\nELECTRIC FIELD CALCULATOR 	by Kyle McIntosh"
      "\n-----------------------------------------------------"
      "\nA simple Python program to calculate the electric"
      "\nfield at a point in 3D space due to one or more"
      "\ncharges."
      "\n"
      "\nCoordinates for X Y and Z should be inputted separated"
      "\nby a ','. To use powers, 'e' must be used in place of"
      "\n'x10'")

data_values = numpy.array([])
display_help()  # Displays the commands on startup
while True:
    try:
        # Input for the users console command
        menu = input("--->").lower()

        if menu == "save" or menu == 's':
            print("Saving data...")
            save_file(data_values)

        elif menu == "open" or menu == 'o':
            print("Select a file to open...")
            data_values = open_file()
            displaydata(data_values)

        elif menu == "input" or menu == 'i':
            data_values = manualinput()
            displaydata(data_values)

        elif menu == "field" or menu == 'f':
            point = electricfield(data_values)

        elif menu == "add":
            for i in range(int(input("How many charges would you like to add:"))):
                try:
                    data = input("Enter the X, Y and Z coordinates for charge q{},\nfollowed by the charge value --:".format(i + 1)).split(',')
                    data_values = numpy.vstack([data_values, data])
                except ValueError:
                    print("Incorrect Data inputted."); i -= 1

            print(data_values)
            displaydata(data_values)

        elif menu == "disp":
            displaydata(data_values)

        elif menu == "help":
            display_help()

        elif menu == "exit":
            sys.exit()

        else:
            print("No option {}.\nType 'help' for list of commands.".format(menu))

    except NameError:
        print("No data currently inputted.")
