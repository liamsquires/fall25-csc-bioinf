from os.path import dirname, join, realpath


def data_dir(subdir):
    return join(dirname(realpath(__file__)), subdir, "data")