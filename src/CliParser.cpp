#include "CliParser.h"

#include "error.h"
#include "CliParameter.h"

namespace mathic {
  namespace {
    // We are using a NameFactory just for its prefix-finding functionality, so
    // we set the thing being created to be just void pointers that are null.
    typedef void* Dummy;
    typedef NameFactory<Dummy> ParamNames;
    ParamNames makeParamNames(std::vector<CliParameter*> params) {
      struct HoldsFunction { // work around for no local functions in old C++
        static std::auto_ptr<Dummy> dummyCreate() {
          return std::auto_ptr<Dummy>(0);
        }
      };

      ParamNames names("option");
      for (size_t i = 0; i < params.size(); ++i)
        names.registerProduct(params[i]->name(), HoldsFunction::dummyCreate);
      return names;
    }
  }

  CliParser::CliParser(): _actions("action") {}

  std::auto_ptr<Action> CliParser::parse(int argc, char** argv) {
    std::vector<std::string> commandLine(argv, argv + argc);
    return parse(commandLine);
  }

  std::auto_ptr<Action> CliParser::parse
    (const std::vector<std::string>& commandLine) {
    if (commandLine.empty())
      throwError<UnknownNameException>("No action specified.");
    std::auto_ptr<Action> action =
      createWithPrefix(_actions, commandLine[0]);

    std::vector<CliParameter*> params;
    action->pushBackParameters(params);
    ParamNames paramNames = makeParamNames(params);

    size_t i = 1;
    while (i < commandLine.size()) {
      std::string const& token = commandLine[i];
      for (++i; i < commandLine.size(); ++i)
        if (!commandLine[i].empty())
          break;
      if (token.empty())
        continue;
      if (token[0] != '-')
        reportError("Expected an option when reading \"" +
                    token + "\", but options start with a dash (-).\n");
      std::string noDash(token.begin() + 1, token.end());
      std::string name = uniqueNameWithPrefix(paramNames, noDash);

      std::string optionArgument;
      if (i < commandLine.size() && commandLine[i][0] != '-') {
        optionArgument = commandLine[i];
        for (++i; i < commandLine.size(); ++i)
          if (!commandLine[i].empty())
            break;
      }

      for (std::vector<CliParameter*>::iterator it = params.begin();; ++it) {
        if (it == params.end()) {
          // shouldn't get here as name was recognized above
          reportInternalError("Processing non-existent option \""
            + name + "\".");
        }
        if ((*it)->name() == name) {
          (*it)->processArgument(optionArgument);
          break;
        }
      }
    }
    return action;
  }
}
