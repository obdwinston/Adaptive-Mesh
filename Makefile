FC = gfortran
FFLAGS = -O2 -march=native -std=f2008
PYTHON = .venv/bin/python3

SRCDIR = src
OBJDIR = build

OBJS = $(OBJDIR)/config.o $(OBJDIR)/grid.o $(OBJDIR)/flow.o $(OBJDIR)/body.o $(OBJDIR)/io.o $(OBJDIR)/main.o

TARGET = main

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(OBJDIR)/config.o: $(SRCDIR)/config.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -J$(OBJDIR)

$(OBJDIR)/grid.o: $(SRCDIR)/grid.f90 $(OBJDIR)/config.o | $(OBJDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -I$(OBJDIR) -J$(OBJDIR)

$(OBJDIR)/flow.o: $(SRCDIR)/flow.f90 $(OBJDIR)/config.o $(OBJDIR)/grid.o | $(OBJDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -I$(OBJDIR) -J$(OBJDIR)

$(OBJDIR)/body.o: $(SRCDIR)/body.f90 $(OBJDIR)/config.o $(OBJDIR)/grid.o | $(OBJDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -I$(OBJDIR) -J$(OBJDIR)

$(OBJDIR)/io.o: $(SRCDIR)/io.f90 $(OBJDIR)/grid.o | $(OBJDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -I$(OBJDIR) -J$(OBJDIR)

$(OBJDIR)/main.o: $(SRCDIR)/main.f90 $(OBJDIR)/config.o $(OBJDIR)/grid.o $(OBJDIR)/flow.o $(OBJDIR)/body.o $(OBJDIR)/io.o | $(OBJDIR)
	$(FC) $(FFLAGS) -c -o $@ $< -I$(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

run: $(TARGET)
	mkdir -p output
	./$(TARGET)

animate:
	$(PYTHON) animate.py

setup:
	python3 -m venv .venv
	.venv/bin/pip install -r requirements.txt

clean:
	rm -rf $(OBJDIR) $(TARGET) output/

.PHONY: clean run animate setup
