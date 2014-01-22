import os
import yaml
from pkg_resources import Requirement, resource_filename

def config():
    pkg = __name__.split(".")[0]
    cfgfile = "{0}.yaml".format(pkg)
    
    global_config_dir = resource_filename(Requirement.parse(pkg), "config")

    paths = [os.curdir, os.path.expanduser("~"), global_config_dir]

    for path in paths:
        for fname in cfgfile, "." + cfgfile:
            try:
                with open(os.path.join(path, fname)) as fo:
                    return yaml.load(fo.read())
            except IOError:
              pass
        
VERSION = "1.0.1"
