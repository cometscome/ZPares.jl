MAKE        = make --no-print-directory
MAKE_DIRS   = src examples
default:
	@for subdir in $(MAKE_DIRS) ; do \
		(cd $$subdir && $(MAKE)) ;\
	done

clean:
	@for subdir in $(MAKE_DIRS) ; do \
		(cd $$subdir && $(MAKE) clean) ;\
	done
