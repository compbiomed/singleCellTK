#source: https://github.com/klmr/decorator

`%@%` = function (decorator, f) UseMethod('%@%')

`%@%.default` = function (decorator, f)
  stop(deparse(substitute(decorator)), ' is not a decorator')

`%@%.decorator` = function (decorator, f) {
  pretty_decorators = as.list(match.call())[-1]
  # Patch delayed decorator.
  if (! is.null({pretty_patched = attr(decorator, 'calls')}))
    pretty_decorators = c(pretty_patched, pretty_decorators[-1])

  # Handle operator precedence so that we can chain decorators.
  if (inherits(f, 'decorator'))
    .delayed_decorate(decorator, f, pretty_decorators)
  else
    prettify(decorator(f), f, pretty_decorators[-length(pretty_decorators)])
}

decorator = function (f)
  structure(f, class = 'decorator')

decorator = decorator(decorator)

print.decorated = function (x, useSource = TRUE, ...) {
  bare = function (f) {
    bare = unclass(f)
    attr(bare, 'decorators') = NULL
    bare
  }

  fun_def = utils::capture.output(print.function(bare(x), useSource = useSource, ...))
  for (decorator in attr(x, 'decorators'))
    cat(deparse(decorator), '%@%\n')
  cat(fun_def, sep = '\n')
  invisible(x)
}

# modules::register_S3_method('print', 'decorated', print.decorated)

prettify = function (f, original, decorator_calls) {
  attr(f, 'srcref') = pretty_code(original)
  attr(f, 'decorators') = decorator_calls
  class(f) = c(class(f), 'decorated')
  f
}

pretty_code = function (f) {
  srcref = attr(f, 'srcref')
  if (is.null(srcref)) body(f) else srcref
}

.delayed_decorate = function (d1, d2, decorator_calls)
  structure(decorator(function (f) d1(d2(f))), calls = decorator_calls)

#'
#' @name simpleLog
#' @title A decorator that prints the arguments to the decorated function
#' @param f A function to decorate
#' @return Prints message
#'
simpleLog <- decorator %@% function(f) {
  function(...) {
    message(deparse(match.call()))
    f(...)
  }
}