#ifndef MATHIC_CLI_PARSER_GUARD
#define MATHIC_CLI_PARSER_GUARD

#include "NameFactory.h"
#include "Action.h"

namespace mathic {
  class CliParser {
  public:
    CliParser();

    template<class ConcreteAction>
    void registerAction(const std::string& name);

    std::auto_ptr<Action> parse(int argc, char** argv);
    std::auto_ptr<Action> parse(const std::vector<std::string>& commandLine);

  private:
    NameFactory<Action> _actions;
  };

  template<class ConcreteAction>
  void registerAction(const std::string& name) {
    nameFactoryRegister<ConcreteAction>(_actions, name);
  };
}

#endif
