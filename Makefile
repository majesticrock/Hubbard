all: PRESET = default
cascadelake: PRESET = cascadelake
icelake: PRESET = icelake
debug: PRESET = debug
test: PRESET = test

all cascadelake icelake debug:
	cmake --preset $(PRESET)
	+cmake --build --preset $(PRESET)

test:
	cmake --preset $(PRESET)
	+cmake --build --preset $(PRESET)
	ctest --preset $(PRESET)

clean:
	rm -rf build
	rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug test