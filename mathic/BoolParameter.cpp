#include <mathic/BoolParameter.h>

#include <mathic/error.h>

namespace mathic {
  BoolParameter::BoolParameter(const std::string& name,
                               const std::string& description,
                               bool defaultValue):
    CliParameter(name, description),
    _value(defaultValue) {
  }

  std::string BoolParameter::argumentType() const {
    return "[BOOL]";
  }

  std::string BoolParameter::valueAsString() const {
    return _value ? "on" : "off";
  }

  void BoolParameter::processArgument(const std::string& argument) {
    if (argument.empty() || argument == "on" || argument == "1")
      _value = true;
    else if (argument == "off" || argument == "0")
      _value = false;
    else {
      reportError("Option -" + name() + " was given the argument \"" +
        argument +
        "\". The only valid arguments are \"on\", \"0\", \"off\" and \"1\".");
    }
  }
}
