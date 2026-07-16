all: PRESET = default
cascadelake: PRESET = cascadelake
icelake: PRESET = icelake
debug: PRESET = debug

all cascadelake icelake debug:
	@cmake --preset $(PRESET)
	@cmake --build --preset $(PRESET)

clean:
	@rm -rf build
	@rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug