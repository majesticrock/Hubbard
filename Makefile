BUILD_DIR = build
MKL_BUILD_DIR = build_MKL
CLUSTER_BUILD_DIR = build_cluster

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCLUSTER_BUILD=OFF ..

MKL: $(MKL_BUILD_DIR)/Makefile
	@$(MAKE) -C $(MKL_BUILD_DIR)

$(MKL_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(MKL_BUILD_DIR)
	@cd $(MKL_BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_MKL=ON ..

cluster: $(CLUSTER_BUILD_DIR)/Makefile
	@$(MAKE) -C $(CLUSTER_BUILD_DIR)

$(CLUSTER_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(CLUSTER_BUILD_DIR)
	@cd $(CLUSTER_BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCLUSTER_BUILD=ON ..

clean:
	@rm -rf $(BUILD_DIR) $(CLUSTER_BUILD_DIR) $(MKL_BUILD_DIR)

.PHONY: all clean cluster MKL
