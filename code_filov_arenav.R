library(pomp)

# Compiling C code and loading the dll
dyn.unload("filov_bat_model.dll")

# compile code
system("R CMD SHLIB filov_bat_model.c")

# load
dyn.load("filov_bat_model.dll")
