import argparse


def read_in_data(file_list,
                 cols_to_drop:list=None,
                 selection_dict:dict=None):
    ds = None
        
    return ds


def make_keogram():
    pass

def make_map():
    pass


def main():
    print('hi')
    pass



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     """Make plots from netCDF files.
        
        
        """)
    parser.add_argument('file_list', metavar='file_list', type=str, nargs='+',
                        help='List of files to plot')
    parser.add_argument('--cols_to_drop', metavar='cols_to_drop', type=str, nargs='+',
                        help='List of columns to drop from the data')
    parser.add_argument('--selection_dict', metavar='selection_dict', type=str, nargs='+',
                        help='Dictionary of selection criteria for the data')
    args = parser.parse_args()
    main(args.file_list, args.cols_to_drop, args.selection_dict)

