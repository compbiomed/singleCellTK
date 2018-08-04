if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package Style", {
    lintr::expect_lint_free(linters = c(lintr::trailing_whitespace_linter,
                                        lintr::no_tab_linter,
                                        lintr::T_and_F_symbol_linter,
                                        lintr::semicolon_terminator_linter,
                                        lintr::infix_spaces_linter,
                                        lintr::closed_curly_linter,
                                        lintr::assignment_linter,
                                        lintr::commas_linter,
                                        lintr::spaces_left_parentheses_linter))
  })
}
