# Kompilator i standard
CXX = g++
CXXFLAGS = -std=c++20 -Wall -Wextra -g -O0

# Nazwa pliku wynikowego
TARGET = Opt

# Pliki źródłowe (wyciągnięte z Twojego zrzutu ekranu i CMake)
SRCS = main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp user_funs.cpp

# Pliki obiektowe (zamienia rozszerzenie .cpp na .o)
OBJS = $(SRCS:.cpp=.o)

# Domyślny cel kompilacji
all: $(TARGET)

# Linkowanie plików obiektowych do pliku wykonywalnego
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Reguła kompilacji poszczególnych plików .cpp do .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Czyszczenie plików tymczasowych
clean:
	rm -f $(OBJS) $(TARGET)

# Oznaczenie celów, które nie są plikami
.PHONY: all clean

