# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Isrc

# --- Application ---
SRCDIR := src
BUILDDIR := build
TARGET := neuron_model

SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SRCS))

# --- Tests ---
TESTSRCDIR := tests
TESTBUILDDIR := build/tests
TEST_TARGET := test_runner

# Test sources are all .cpp files in the tests directory
TEST_SRCS := $(wildcard $(TESTSRCDIR)/*.cpp)
TEST_OBJS := $(patsubst $(TESTSRCDIR)/%.cpp,$(TESTBUILDDIR)/%.o,$(TEST_SRCS))

# Test dependencies are all source files from the app, EXCEPT for main.cpp
APP_MAIN_SRC := $(SRCDIR)/main.cpp
LIB_SRCS := $(filter-out $(APP_MAIN_SRC), $(wildcard $(SRCDIR)/*.cpp))
LIB_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(TESTBUILDDIR)/%.o,$(LIB_SRCS))

# --- Rules ---

# Default rule
all: $(TARGET)

# Main application linking
$(TARGET): $(OBJS)
	@echo "Linking main application..."
	$(CXX) $(CXXFLAGS) -o $@ $^

# Test runner linking
$(TEST_TARGET): $(TEST_OBJS) $(LIB_OBJS)
	@echo "Linking test runner..."
	$(CXX) $(CXXFLAGS) -o $@ $^

# App object files compilation
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Test source object files compilation
$(TESTBUILDDIR)/%.o: $(TESTSRCDIR)/%.cpp | $(TESTBUILDDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Test dependency object files compilation
$(TESTBUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(TESTBUILDDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Create directories
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(TESTBUILDDIR):
	mkdir -p $(TESTBUILDDIR)

# --- Phony Targets ---
.PHONY: all test clean

# Test execution target
test: $(TEST_TARGET)
	@echo "Running tests..."
	./$(TEST_TARGET)

# Clean up build artifacts
clean:
	@echo "Cleaning up..."
	rm -rf $(BUILDDIR) $(TARGET) $(TEST_TARGET)
