PACKAGE     := CAFAna
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
SUBDIRS     := Core Vars Experiment
LIBCXXFILES := $(wildcard *.cxx)

########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_ifdhc.mk

LIBLINK    := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -l$(PACKAGE) -lCAFAnaCore -lCAFAnaVars -lCAFAnaExperiment

override LIBLIBS += `root-config --libs` -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) -L$(IFDHC_FQ_DIR)/lib/ -lifdh
