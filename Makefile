
FC=R
LAMBDAR="~/Research/UWA/LAMBDAR/Current_Dev/"
opts="-l"
LIB="/gamah/awright/src/R/lib/"

sudoinstall:
	sudo $(FC) CMD INSTALL $(LAMBDAR)
install:
	$(FC) CMD INSTALL $(LAMBDAR) $(opts) $(LIB)

build:
	$(FC) CMD build $(LAMBDAR) 

check:
	$(FC) CMD check $(LAMBDAR) $(opts) $(LIB)
