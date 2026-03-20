include config.mk

.DEFAULT_GOAL := all

PROGRAMS := \
	$(BIN)/jamming_by_growth \
	$(BIN)/jamming_by_growth_lineage \
	$(BIN)/box_compress_bext \
	$(BIN)/shear_yeast_linearshear

$(BIN)/jamming_by_growth: $(SRC)/jamming_by_growth.f
$(BIN)/jamming_by_growth_lineage: $(SRC)/jamming_by_growth_lineage.f
$(BIN)/box_compress_bext: $(SRC)/box_compress_bext.f
$(BIN)/shear_yeast_linearshear: $(SRC)/shear_yeast_linearshear.f

jamming_by_growth: $(BIN)/jamming_by_growth
jamming_by_growth_lineage: $(BIN)/jamming_by_growth_lineage
box_compress_bext: $(BIN)/box_compress_bext
shear_yeast_linearshear: $(BIN)/shear_yeast_linearshear

$(PROGRAMS):
	@echo "BUILDING $@"
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $< -o $@
	@echo "BUILDING IS DONE"

.PHONY: all clean jamming_by_growth jamming_by_growth_lineage box_compress_bext shear_yeast_linearshear

all: $(PROGRAMS)

clean:
	@echo "Cleaning..."
	$(RM) $(PROGRAMS)
