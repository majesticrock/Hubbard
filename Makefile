BUILD_DIR = build

all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

# Run cmake to configure the project with custom install prefix
$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

# Clean up the build directory
clean:
	@rm -rf $(BUILD_DIR)

.PHONY: all clean
