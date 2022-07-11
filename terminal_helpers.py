"""
This file defines terminal printout functions used throughout the package.
"""


def terminal_table(mpc_code, band, good_epochs):
    """
    Generates an output table of source ids, frames indicies, and flux values,
    for a given set of epochs and a given band set.
    
    Arguments: mpc_code (str) -- the mpc designation for the asteroid
               band (int) -- the w band being targeted
               good_epochs (list) -- a list of valid source ids
    Returns: None (terminal printout)
    """

    band_file_map = {2: ("input_data/" + mpc_code + ".txt", 
                    "input_data/" + mpc_code + ".tbl", 0),
                    3: ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_3band.tbl", 22),
                    4: ("input_data/" + mpc_code + ".txt",
                    "input_data/" + mpc_code + "_cryo.tbl", 44)}

    # Table formatting below
    dash_string = '-' * 65 + '-' * band_file_map[band][2]
    test_string = "| Frame | Source Id "
    for i in range(1, band + 1):
        test_string += f"|  W{i} Flux | W{i} Sigma "
    test_string += '|'
    print(dash_string)
    print(test_string)
    print(dash_string)
    for idx, epoch in enumerate(good_epochs):
        epoch_string = '|' + f"{idx * 2 + 1}-{idx * 2 + 2}".rjust(6)
        epoch_string += f" | {epoch}"
        for i in range(int(band) * 2):
            epoch_string += ' | ' + f"{good_epochs[epoch][4 + i]}".rjust(8) 
        epoch_string += " | "
        print(epoch_string)
    print(dash_string)