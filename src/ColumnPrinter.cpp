#include "ColumnPrinter.h"

namespace mathic {
  namespace {
	size_t getLineWidth(const std::string& str, size_t pos) {
	  size_t startPos = pos;
	  while (pos < str.size() && str[pos] != '\n')
		++pos;
	  return pos - startPos;
	}

	void printSpaces(std::ostream& out, size_t howMany) {
	  while (howMany > 0) {
		out << ' ';
		--howMany;
	  }
	}
  }

  ColumnPrinter::ColumnPrinter(size_t columnCount):
	_cols() {
	while (columnCount > 0) {
	  addColumn();
	  --columnCount;
	}
  }

  ColumnPrinter::~ColumnPrinter() {
	for (std::vector<Col*>::iterator it = _cols.begin();
		 it != _cols.end(); ++it)
	  delete *it;
  }

  void ColumnPrinter::setPrefix(const std::string& prefix) {
	_prefix = prefix;
  }

  void ColumnPrinter::addColumn(bool flushLeft,
								const std::string& prefix,
								const std::string& suffix) {
	std::auto_ptr<Col> col(new Col());
	col->prefix = prefix;
	col->suffix = suffix;
	col->flushLeft = flushLeft;

	_cols.push_back(0);
	_cols.back() = col.release(); // push_back didn't throw, so safe to release
  }

  size_t ColumnPrinter::getColumnCount() const {
	return _cols.size();
  }

  std::ostream& ColumnPrinter::operator[](size_t col) {
	MATHIC_ASSERT(col < getColumnCount());
	return _cols[col]->text;
  }

  void ColumnPrinter::print(std::ostream& out) const {
	// Calculate the width of each column.
	std::vector<size_t> widths(getColumnCount());
	for (size_t col = 0; col < getColumnCount(); ++col) {
	  const std::string& text = _cols[col]->text.str();
	  size_t maxWidth = 0;

	  size_t pos = 0;
	  while (pos < text.size()) {
		size_t width = getLineWidth(text, pos);
		if (width > maxWidth)
		  maxWidth = width;

		// We can't just increment pos unconditionally by width + 1, as
		// that could result in an overflow.
		pos += width;
		if (text[pos] == '\n')
		  ++pos;
	  }
	  widths[col] = maxWidth;
	}

	// Print each row
	std::vector<size_t> poses(getColumnCount());
	while (true) {
	  bool done = true;
	  for (size_t col = 0; col < getColumnCount(); ++col) {
		if (poses[col] < _cols[col]->text.str().size()) {
		  done = false;
		  break;
		}
	  }
	  if (done)
		break;

	  out << _prefix;
	  for (size_t col = 0; col < getColumnCount(); ++col) {
		out << _cols[col]->prefix;

		const std::string& text = _cols[col]->text.str();
		size_t& pos = poses[col];
		size_t width = getLineWidth(text, pos);

		if (!_cols[col]->flushLeft)
		  printSpaces(out, widths[col] - width);

		while (pos < text.size()) {
		  if (text[pos] == '\n') {
			++pos;
			break;
		  }
		  out << text[pos];
		  ++pos;
		}

		if (_cols[col]->flushLeft)
		  printSpaces(out, widths[col] - width);

		out << _cols[col]->suffix;
	  }
	  out << '\n';
	}
  }

  std::ostream& operator<<(std::ostream& out, const ColumnPrinter& printer) {
	printer.print(out);
	return out;
  }

  void print(FILE* out, const ColumnPrinter& pr) {
	std::ostringstream str;
	str << pr;
	fputs(str.str().c_str(), out);
  }
}