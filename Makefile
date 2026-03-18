include config.mk

.DEFAULT_GOAL := all

PROGRAMS := \
	$(BIN)/jamming_by_growth \
	$(BIN)/shear_yeast_linearshear

$(BIN)/jamming_by_growth: $(SRC)/jamming_by_growth.f
$(BIN)/shear_yeast_linearshear: $(SRC)/shear_yeast_linearshear.f

$(PROGRAMS):
	@echo "BUILDING $@"
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $< -o $@
	@echo "BUILDING IS DONE"

.PHONY: all clean

all: $(PROGRAMS)

clean:
	@echo "Cleaning..."
	$(RM) $(PROGRAMS)
