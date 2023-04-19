# Misc Utility Functions


def str_to_ut(in_str):
    """Convert a string to a datetime object.

    Parameters
    ----------
    in_str : str
        String to convert to datetime object.

    Returns
    -------
    datetime
        Datetime object.

    """

    import datetime

    return datetime.datetime.strptime(
        in_str.ljust(14, '0'), '%Y%m%d%H%M%S')
