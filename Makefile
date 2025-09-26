# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -O3 -Wall -Isrc

# CUDA compiler and flags
# Note: You may need to change sm_70 to match your GPU's compute capability.
NVCC := nvcc
NVCCFLAGS := -std=c++17 -O3 -arch=sm_70 --ptxas-options=-v -Isrc
LDFLAGS := -lcudart # Link against the CUDA runtime library

# --- Application ---
SRCDIR := src
BUILDDIR := build
TARGET := neuron_model

# Differentiate source file types
CXX_SRCS := $(wildcard $(SRCDIR)/*.cpp)
CUDA_SRCS := $(wildcard $(SRCDIR)/*.cu)

# Generate object file paths
CXX_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(CXX_SRCS))
CUDA_OBJS := $(patsubst $(SRCDIR)/%.cu,$(BUILDDIR)/%.o,$(CUDA_SRCS))
OBJS := $(CXX_OBJS) $(CUDA_OBJS)

# --- Tests ---
TESTSRCDIR := tests
TESTBUILDDIR := build/tests
TEST_TARGET := test_runner

# Test sources are all .cpp files in the tests directory
TEST_SRCS := $(wildcard $(TESTSRCDIR)/*.cpp)
TEST_OBJS := $(patsubst $(TESTSRCDIR)/%.cpp,$(TESTBUILDDIR)/%.o,$(TEST_SRCS))

# Test dependencies are all C++ source files from the app, EXCLUDING main.cpp
APP_MAIN_SRC := $(SRCDIR)/main.cpp
# Note: Test build intentionally excludes .cu files
LIB_CXX_SRCS := $(filter-out $(APP_MAIN_SRC), $(CXX_SRCS))
LIB_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(TESTBUILDDIR)/%.o,$(LIB_CXX_SRCS))

# --- Rules ---

.PHONY: all test clean

# Default rule
all: $(TARGET)

# Main application linking
$(TARGET): $(OBJS)
	@echo "Linking main application: $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Test runner linking (no CUDA dependencies)
test: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS) $(LIB_OBJS)
	@echo "Linking test runner: $(TEST_TARGET)..."
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compilation rule for C++ files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	@echo "Compiling C++: $< -> $@"
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compilation rule for CUDA files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cu | $(BUILDDIR)
	@echo "Compiling CUDA: $< -> $@"
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

# Compilation rule for test source files
$(TESTBUILDDIR)/%.o: $(TESTSRCDIR)/%.cpp | $(TESTBUILDDIR)
	@echo "Compiling Test C++: $< -> $@"
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Compilation rule for library files used by tests
$(TESTBUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(TESTBUILDDIR)
	@echo "Compiling Lib for Test: $< -> $@"
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Create build directories
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(TESTBUILDDIR):
	mkdir -p $(TESTBUILDDIR)

# Clean up build artifacts
clean:
	@echo "Cleaning up build artifacts..."
	rm -rf $(BUILDDIR) $(TARGET) $(TEST_TARGET)