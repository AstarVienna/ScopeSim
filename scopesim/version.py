from importlib import metadata
try:
    version = metadata.version(__package__)
except:
    version = 0.8
