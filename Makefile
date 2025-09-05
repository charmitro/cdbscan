CC = cc
AR = ar
CFLAGS = -Wall -O2 -fPIC -Iinclude
PREFIX = /usr/local

all: libcdbscan.a libcdbscan.so

libcdbscan.a: src/cdbscan.o
	$(AR) rcs $@ $^

libcdbscan.so: src/cdbscan.o
	$(CC) -shared -o $@ $^ -lm $(LDFLAGS)

src/cdbscan.o: src/cdbscan.c include/cdbscan.h
	$(CC) $(CFLAGS) -c -o $@ $<

examples: examples/example examples/example_distances examples/example_normalize examples/example_estimate_eps examples/example_kdtree

examples/example: examples/example.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

examples/example_distances: examples/example_distances.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

examples/example_normalize: examples/example_normalize.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

examples/example_estimate_eps: examples/example_estimate_eps.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

examples/example_kdtree: examples/example_kdtree.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

install: libcdbscan.a libcdbscan.so
	install -d $(DESTDIR)$(PREFIX)/lib
	install -d $(DESTDIR)$(PREFIX)/include
	install -m 644 libcdbscan.a $(DESTDIR)$(PREFIX)/lib/
	install -m 755 libcdbscan.so $(DESTDIR)$(PREFIX)/lib/
	install -m 644 include/cdbscan.h $(DESTDIR)$(PREFIX)/include/

tests: tests/test_core_points tests/test_density_reachability tests/test_border_noise tests/test_cluster_properties tests/test_kdtree

tests/test_core_points: tests/test_core_points.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

tests/test_density_reachability: tests/test_density_reachability.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

tests/test_border_noise: tests/test_border_noise.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

tests/test_cluster_properties: tests/test_cluster_properties.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

tests/test_kdtree: tests/test_kdtree.c libcdbscan.a
	$(CC) $(CFLAGS) -o $@ $< -L. -lcdbscan -lm $(LDFLAGS)

test: tests
	@echo "Running specification tests..."
	@echo "=============================="
	@LD_LIBRARY_PATH=.:$$LD_LIBRARY_PATH ./tests/test_core_points
	@echo
	@LD_LIBRARY_PATH=.:$$LD_LIBRARY_PATH ./tests/test_density_reachability
	@echo
	@LD_LIBRARY_PATH=.:$$LD_LIBRARY_PATH ./tests/test_border_noise
	@echo
	@LD_LIBRARY_PATH=.:$$LD_LIBRARY_PATH ./tests/test_cluster_properties
	@echo
	@echo "[SUCCESS] All specification tests passed!"

format:
	@echo "Formatting C source files..."
	@clang-format -i src/*.c include/*.h examples/*.c tests/*.c
	@echo "Formatting complete."

clean:
	rm -f libcdbscan.a libcdbscan.so src/*.o
	rm -f examples/example examples/example_distances examples/example_normalize examples/example_estimate_eps examples/example_kdtree
	rm -f tests/test_core_points tests/test_density_reachability tests/test_border_noise tests/test_cluster_properties tests/test_kdtree

.PHONY: all install clean examples tests test format
