"""Here we will try to make all of the available plots... 
Let's see if it's possible!
"""


"""
- first get the arguments.
- then read in data specified
    - this will be either gitm, sami, or both
- then get the data we want 
    - e.g. diff from filter, diff btwn runs, etc.
- then launch plotting routines
    - keos, maps, polar dials, etc.
    - save plots as we go, print out dirs for ffmpeg
        - and/or a ffmpeg command.
- exit gracefully?

OR

- get the args
- launch plotting routines for fig types
- get the data for each plotting routine
    - loop thru making plots
    - make ffmpeg commands and/or print them out
    - return data if it's needed elsewhere
- exit :)

"""


def main(args):

    # print out help info
    if args.var_help:
        # go thru each of the paths provided and print out data vars & shapes
        for help_path in args.var_help:
            print(f'help for {help_path}')
            print('----------------')
            print('')

            # get the data
            data = get_data(help_path)

            for var in data['vars']:
                print(var, data[var].shape)

    # set up a keyword to give data back if multiple plots requested.
    if np.sum([args.keogram, args.map, args.polardials]) > 1:
        return_data = True
    else:
        return_data = False

    if args.keogram:
        pass
    if args.map:
        pass
    if args.polardials:
        pass
