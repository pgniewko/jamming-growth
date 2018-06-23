include config.mk

TARGET      := $(BIN)/jamming_by_growth
SHEAR       := $(BIN)/shear_yeast_linearshear

JAMMING_SOURCES     := $(SRC)/jamming_by_growth.f
SHEAR_SOURCES       := $(SRC)/shear_yeast_linearshear.f

$(TARGET):$(JAMMING_SOURCES) 
	@echo BUILDING MAIN JAMMING CODE
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@ 
	@echo BUILDING IS DONE

$(SHEAR):$(SHEAR_SOURCES)
	@echo BUILDING SHEAR CALCULATIONS CODE
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

# Tell make that these are phony targets
.PHONY: all clean

all: $(TARGET) $(SHEAR)
clean:
	@echo Cleaning...
	rm -f $(TARGET) $(SHEAR)


