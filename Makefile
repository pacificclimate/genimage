COREDIR = core
TOOLSDIR = tools

.PHONY: core
core:
	$(MAKE) -C $(COREDIR);

.PHONY: tools
tools:
	$(MAKE) -C $(TOOLSDIR);

clean:
	$(MAKE) -C $(COREDIR) clean
	$(MAKE) -C $(TOOLSDIR) clean
