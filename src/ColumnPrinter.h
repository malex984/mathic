#ifndef MATHIC_COLUMN_PRINTER_GUARD
#define MATHIC_COLUMN_PRINTER_GUARD

#include "stdinc.h"
#include <sstream>
#include <vector>
#include <cstdio>
#include <memory>
#include <string>

namespace mathic {
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

    static std::string commafy(unsigned long long l);
    static std::string percent(
      unsigned long long numerator,
      unsigned long long denominator);

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
}

#endif
