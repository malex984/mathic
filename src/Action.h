#ifndef MATHIC_ACTION_GUARD
#define MATHIC_ACTION_GUARD

#include "stdinc.h"
#include <vector>

namespace mathic {
  class CliParameter;

  class Action {
  public:
    virtual ~Action();

    // Do what it is this action does.
    virtual void performAction() = 0;

    // ***************************************
    // **** Information provided by each Action

    // The name of the action.
    virtual const char* name() const = 0;

    // More detailed explanation of what the action does.
    virtual const char* description() const = 0;

    // One-line summary of the description.
    virtual const char* shortDescription() const = 0;

    // Append the parameters for this action to the passed-in container.
    // Do not clear the passed-in container.
    virtual void pushBackParameters(std::vector<CliParameter*>& parameters);
  };
}

#endif
