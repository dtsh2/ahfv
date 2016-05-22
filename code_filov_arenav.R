library(pomp)

# filovirus - bat model
# Compiling C code and loading the dll
dyn.unload("filov_bat_model.dll")

# compile code
system("R CMD SHLIB filov_bat_model.c")

# load
dyn.load("filov_bat_model.dll")

# arenavirus - rodent model
# Compiling C code and loading the dll
dyn.unload("arenav_rodent_model.dll")

# compile code
system("R CMD SHLIB arenav_rodent_model.c")

# load
dyn.load("arenav_rodent_model.dll")

