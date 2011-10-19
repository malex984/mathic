#ifndef COLUMN_PRINTER_GUARD
#define COLUMN_PRINTER_GUARD

#include <sstream>
#include <vector>
#include <cstdio>

class ColumnPrinter {
 public:
  ColumnPrinter(size_t columnCount = 0);
  ~ColumnPrinter();

  void setPrefix(const std::string& prefix);
  void addColumn(bool flushLeft = true,
				 const std::string& prefix = "  ",
                 const std::string& suffix = "");
  size_t getColumnCount() const;

  std::ostream& operator[](size_t col);

  void print(std::ostream& out) const;

 private:
  struct Col {
	std::string prefix;
	std::stringstream text;
	std::string suffix;
	bool flushLeft;
  };
  std::vector<Col*> _cols;
  std::string _prefix;
};

std::ostream& operator<<(std::ostream& out, const ColumnPrinter& printer);
void print(FILE* out, const ColumnPrinter& pr);

#endif
