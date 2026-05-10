#ifndef _instantiate_matrix_evaluator_hpp_
#define _instantiate_matrix_evaluator_hpp_

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cctype>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <traits.hpp>

#include "nrgchain.hpp"

namespace NRG::Spawn {

struct MatrixEvalContext {
  std::map<std::string, double> variables;
  NRG::Tools::NrgChain::WilsonData wilson;
};

class MatrixEvaluator {
 private:
  struct Value {
    bool scalar = true;
    double number = 0.0;
    std::vector<Value> items;

    static Value scalar_value(const double value) { return Value{true, value, {}}; }
    static Value list_value(std::vector<Value> values) { return Value{false, 0.0, std::move(values)}; }
  };

  MatrixEvalContext context;
  std::string input;
  std::string source;
  size_t pos = 0;
  size_t line = 1;
  size_t column = 1;

  [[nodiscard]] auto error(const std::string &message) const {
    return std::runtime_error(fmt::format("{}:{}:{}: {}", source, line, column, message));
  }

  [[nodiscard]] bool eof() const noexcept { return pos >= input.size(); }
  [[nodiscard]] char peek() const noexcept { return eof() ? '\0' : input[pos]; }

  char get() {
    const auto ch = peek();
    if (eof()) return ch;
    ++pos;
    if (ch == '\n') {
      ++line;
      column = 1;
    } else {
      ++column;
    }
    return ch;
  }

  void skip_ignored() {
    while (!eof()) {
      while (!eof() && std::isspace(static_cast<unsigned char>(peek()))) get();
      if (peek() != '#') return;
      while (!eof() && peek() != '\n') get();
    }
  }

  bool consume(const char expected) {
    skip_ignored();
    if (peek() != expected) return false;
    get();
    return true;
  }

  void expect(const char expected) {
    if (!consume(expected)) throw error(fmt::format("Expected '{}'.", expected));
  }

  [[nodiscard]] bool starts_identifier() const {
    const auto ch = peek();
    return std::isalpha(static_cast<unsigned char>(ch)) || ch == '_';
  }

  std::string parse_identifier() {
    skip_ignored();
    if (!starts_identifier()) throw error("Expected identifier.");
    std::string name;
    while (!eof()) {
      const auto ch = peek();
      if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '_') break;
      name.push_back(get());
    }
    return name;
  }

  double parse_number() {
    skip_ignored();
    const auto start = input.c_str() + pos;
    char *end = nullptr;
    const auto value = std::strtod(start, &end);
    if (end == start) throw error("Expected number.");
    while (input.c_str() + pos != end) get();
    if (!std::isfinite(value)) throw error("Non-finite number.");
    return value;
  }

  double lookup_variable(const std::string &name) const {
    if (name == "Pi" || name == "pi") return M_PI;
    const auto it = context.variables.find(name);
    if (it == context.variables.end()) throw error("Unknown variable: " + name);
    return it->second;
  }

  const NRG::Tools::NrgChain::WilsonChannel &channel(const int one_based_channel) const {
    if (one_based_channel < 1 || static_cast<size_t>(one_based_channel) > context.wilson.channels.size())
      throw error(fmt::format("Channel {} is out of range.", one_based_channel));
    return context.wilson.channels[static_cast<size_t>(one_based_channel - 1)];
  }

  double coefficient(const std::vector<double> &values, const int index) const {
    if (index < 0 || static_cast<size_t>(index) >= values.size())
      throw error(fmt::format("Coefficient index {} is out of range.", index));
    return values[static_cast<size_t>(index)];
  }

  static int integral_argument(const double value, const std::string &name) {
    const auto rounded = std::round(value);
    if (std::abs(value - rounded) > 1e-12)
      throw std::runtime_error(fmt::format("{} argument must be integral.", name));
    return static_cast<int>(rounded);
  }

  double call_function(const std::string &name, const std::vector<double> &args) const {
    if ((name == "Sqrt" || name == "sqrt") && args.size() == 1) return std::sqrt(args[0]);
    if (name == "exp" && args.size() == 1) return std::exp(args[0]);
    if (name == "log" && args.size() == 1) return std::log(args[0]);
    if (name == "Conjugate" && args.size() == 1) return args[0];
    if (name == "gammaPolCh" && args.size() == 1) {
      const auto ch = integral_argument(args[0], name);
      return std::sqrt(channel(ch).theta / M_PI);
    }
    if (name == "coefxi" && args.size() == 2) {
      const auto ch = integral_argument(args[0], name);
      const auto idx = integral_argument(args[1], name);
      return coefficient(channel(ch).xi, idx);
    }
    if (name == "coefzeta" && args.size() == 2) {
      const auto ch = integral_argument(args[0], name);
      const auto idx = integral_argument(args[1], name);
      return coefficient(channel(ch).zeta, idx);
    }
    throw error(fmt::format("Unsupported function {} with {} argument(s).", name, args.size()));
  }

  std::vector<double> parse_arguments(const char closing) {
    std::vector<double> args;
    skip_ignored();
    if (consume(closing)) return args;
    while (true) {
      args.push_back(parse_expression());
      if (consume(closing)) return args;
      expect(',');
    }
  }

  double parse_primary() {
    skip_ignored();
    if (consume('(')) {
      const auto value = parse_expression();
      expect(')');
      return value;
    }
    if (starts_identifier()) {
      const auto name = parse_identifier();
      if (consume('(')) return call_function(name, parse_arguments(')'));
      if (consume('[')) return call_function(name, parse_arguments(']'));
      if (name == "parse") throw error("parse directives are not supported by the in-process matrix evaluator yet.");
      return lookup_variable(name);
    }
    return parse_number();
  }

  double parse_unary() {
    skip_ignored();
    if (consume('+')) return parse_unary();
    if (consume('-')) return -parse_unary();
    return parse_primary();
  }

  double parse_product() {
    auto value = parse_unary();
    while (true) {
      if (consume('*')) {
        value *= parse_unary();
      } else if (consume('/')) {
        const auto denom = parse_unary();
        if (denom == 0.0) throw error("Division by zero.");
        value /= denom;
      } else {
        return value;
      }
    }
  }

  double parse_expression() {
    auto value = parse_product();
    while (true) {
      if (consume('+')) {
        value += parse_product();
      } else if (consume('-')) {
        value -= parse_product();
      } else {
        return value;
      }
    }
  }

  Value parse_value() {
    skip_ignored();
    if (!consume('{')) return Value::scalar_value(parse_expression());

    std::vector<Value> items;
    skip_ignored();
    if (consume('}')) return Value::list_value(std::move(items));
    while (true) {
      items.push_back(parse_value());
      if (consume('}')) return Value::list_value(std::move(items));
      expect(',');
    }
  }

  static double scalar(const Value &value) {
    if (!value.scalar) throw std::runtime_error("Scalar value expected inside matrix.");
    return value.number;
  }

  static Matrix_traits<double> to_matrix(const Value &value) {
    if (value.scalar) {
      Matrix_traits<double> matrix(1, 1);
      matrix(0, 0) = value.number;
      return matrix;
    }
    if (value.items.empty()) return Matrix_traits<double>(0, 0);

    const bool matrix_form = !value.items.front().scalar;
    if (!matrix_form) {
      Matrix_traits<double> matrix(1, static_cast<::Eigen::Index>(value.items.size()));
      for (size_t col = 0; col < value.items.size(); ++col) matrix(0, static_cast<::Eigen::Index>(col)) = scalar(value.items[col]);
      return matrix;
    }

    const auto cols = value.items.front().items.size();
    Matrix_traits<double> matrix(static_cast<::Eigen::Index>(value.items.size()), static_cast<::Eigen::Index>(cols));
    for (size_t row = 0; row < value.items.size(); ++row) {
      if (value.items[row].scalar) throw std::runtime_error("Cannot mix scalars and rows in a matrix.");
      if (value.items[row].items.size() != cols) throw std::runtime_error("Matrix rows must have equal length.");
      for (size_t col = 0; col < cols; ++col) matrix(static_cast<::Eigen::Index>(row), static_cast<::Eigen::Index>(col)) = scalar(value.items[row].items[col]);
    }
    return matrix;
  }

  static void reject_shell_escapes(const std::string &text, const std::string &source_name) {
    size_t line_no = 1;
    size_t start = 0;
    while (start <= text.size()) {
      const auto end = text.find('\n', start);
      const auto line_text = text.substr(start, end == std::string::npos ? std::string::npos : end - start);
      const auto first = line_text.find_first_not_of(" \t\r");
      if (first != std::string::npos && line_text[first] == '!')
        throw std::runtime_error(fmt::format("{}:{}: shell escapes are not supported in matrix templates.", source_name, line_no));
      if (end == std::string::npos) break;
      start = end + 1;
      ++line_no;
    }
  }

 public:
  explicit MatrixEvaluator(MatrixEvalContext context_) : context(std::move(context_)) {}

  Matrix_traits<double> evaluate_matrix(std::string text, std::string source_name = "<matrix>") {
    reject_shell_escapes(text, source_name);
    input = std::move(text);
    source = std::move(source_name);
    pos = 0;
    line = 1;
    column = 1;
    auto value = parse_value();
    skip_ignored();
    if (!eof()) throw error("Unexpected trailing input.");
    return to_matrix(value);
  }
};

} // namespace NRG::Spawn

#endif
