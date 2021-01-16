import MDAnalysis as mda


def make_index(traj_filename, top_filename):

    # Initialise with the tpr and xtc files
    u = mda.Universe(top_filename, traj_filename)

    print("You are missing an index file, would you like"
          "to create one now? "
          "\n(This will be one group with"
          "all atoms in the system")

    index_filename = input("Enter a name for index file >")
    # index_filename = "test.ndx"

    if index_filename == "":
        print("No name entered, using default name 'index.ndx'")
        index_filename = "index.ndx"

    if index_filename[-4:] != ".ndx":
        index_filename += ".ndx"

    print(index_filename)

    # Write all atoms to index file
    with mda.selections.gromacs.SelectionWriter(index_filename, mode="w") as ndx:
        ndx.write(u.atoms)

    return index_filename

