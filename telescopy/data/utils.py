"""
Put utility functions for accessing the data folder here
Possibly also add the skycalc_cli interface here

"""

def data_database():
    pass
    
def set_data_directory():
    pass
    
def set_cwd():
    pass
    
def get_cwd():
    pass
    
def pwd():
    return get_cwd()

def get_local_packages(path=None):

    if path is None:
        path = os.path.join(dkeys.PKG_DIR,
                            dkeys.INST_PKG_LOCAL_PATH,
                            dkeys.INST_PKG_LOCAL_DB_NAME)

    if not os.path.exists(path):
        raise ValueError(path + " doesn't exist")

    list_of_packages = ioascii.read(path)

    return list_of_packages



