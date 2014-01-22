import os
import yaml
from pkg_resources import Requirement, resource_filename

# XDG compliant configuration directories
HOME = os.path.expanduser("~")
XDG_CONFIG_HOME = os.environ.get("XDG_CONFIG_HOME", os.path.join(HOME, ".config"))
XDG_CONFIG_DIRS = [XDG_CONFIG_HOME] + os.environ.get("XDG_CONFIG_DIRS", "/etc/xdg").split(":")

def config():
    pkg = __name__.split(".")[0]
    cfgfile = "{0}.yaml".format(pkg)
    
    global_config_dir = resource_filename(Requirement.parse(pkg), "config")

    paths = [os.curdir] + [os.path.join(d, pkg) for d in XDG_CONFIG_DIRS] + [global_config_dir]
    print paths
    for path in paths:
        try:
            with open(os.path.join(path, cfgfile)) as fo:
                return yaml.load(fo.read())
        except IOError:
            pass
    raise Exception, "Configuration file {0} not found".format(cfgfile)

VERSION = "1.0.2"
