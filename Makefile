include config.mk

TARGET1      := $(BIN)/cg_cj_finder
TARGET2      := $(BIN)/growth_2d
TARGET3      := $(BIN)/dumbbell
TARGET4      := $(BIN)/growth_extra
TARGET5      := $(BIN)/growth_2d_rndremover
TARGET6      := $(BIN)/jamming_by_growth
TARGET7      := $(BIN)/shear_yeast_linearshear
TARGET8      := $(BIN)/jamming_by_growth_adaptivedt

SOURCES1	 := $(SRC)/CG_CJ_finder_LS.f
SOURCES2     := $(SRC)/growth_jamming.f
SOURCES3     := $(SRC)/dumbbell_jam.f
SOURCES4     := $(SRC)/growth_extra.f
SOURCES5     := $(SRC)/growth_jamming_ext.f
SOURCES6     := $(SRC)/jamming_by_growth.f
SOURCES7     := $(SRC)/shear_yeast_linearshear.f
SOURCES8     := $(SRC)/jamming_by_growth_adaptivedt.f

$(TARGET1):$(SOURCES1) 
	@echo BUILDING TARGET 1
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@ 
	@echo BUILDING IS DONE

$(TARGET2):$(SOURCES2)
	@echo BUILDING TARGET 2
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

$(TARGET3):$(SOURCES3)
	@echo BUILDING TARGET 3
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE
	
$(TARGET4):$(SOURCES4)
	@echo BUILDING TARGET 4
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

$(TARGET5):$(SOURCES5)
	@echo BUILDING TARGET 5
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

$(TARGET6):$(SOURCES6)
	@echo BUILDING TARGET 6
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

$(TARGET7):$(SOURCES7)
	@echo BUILDING TARGET 7
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

$(TARGET8):$(SOURCES8)
	@echo BUILDING TARGET 8
	@mkdir -p $(@D)
	$(FC) $(FTNFLAGS) $^ -o $@
	@echo BUILDING IS DONE

# Tell make that these are phony targets
.PHONY: all build clean

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8)

clean:
	@echo Cleaning...
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGER7) $(TARGET8)



