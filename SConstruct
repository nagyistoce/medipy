import os.path

from deployement.scons import Environment

vars = Variables()
vars.Add(BoolVariable("optimized", "Set to 1 to build for release", 0))
vars.Add(BoolVariable("verbose", "Set to 1 to print build commands", 0))
vars.Add(PathVariable("build_dir", "Location of the build directory", 
                      os.path.abspath("build/medipy"),
                      PathVariable.PathIsDirCreate))

env = Environment(variables=vars, ENV = {"PATH" : os.environ["PATH"]})
Help(vars.GenerateHelpText(env))
Export("env")

# Add lib_dir to targets if not specified on the command line. Necessary, as it
# is outside the build tree.
#if env["lib_dir"] not in BUILD_TARGETS :
#    BUILD_TARGETS.append(env["lib_dir"])

dirs = ["lib", "plugins", "apps"]

if env["build_dir"] != "." :
    for dir in dirs :
        SConscript(os.path.join(dir, "SConstruct"), 
                   variant_dir=os.path.join(env["build_dir"], dir))
else :
    for dir in dirs :
        SConscript(os.path.join(dir, "SConstruct"))
