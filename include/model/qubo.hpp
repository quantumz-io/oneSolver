/*! file qubo.hpp
 * 
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "helpers/hash.hpp"
#include "helpers/insert.hpp"
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>

#include <boost/bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/qi.hpp>
#include <iostream>

namespace qi = boost::spirit::qi;

#ifndef QUBO_MODEL_HPP__
#define QUBO_MODEL_HPP__

namespace qubo {

using boost::fusion::at_c;
using qi::double_;
using qi::int_;
using qi::phrase_parse;
/**
 * @brief Storage class for linear coefficients of the QUBO model.
 * 
 * @tparam NodeType integer number type used to index nodes of the QUBO model.
 * @tparam CoefType real type used to store values of the QUBO model.
 */
template <class NodeType, class CoefType>
using LinearCoef = std::unordered_map<NodeType, CoefType>;

/**
 * @brief Storage class for quadratic coefficients of the QUBO model.
 * 
 * @tparam NodeType integer number type used to index nodes of the QUBO model.
 * @tparam CoefType real type used to store values of the QUBO model.
 */
template <class NodeType, class CoefType>
using QuadraticCoef = std::unordered_map<std::pair<NodeType, NodeType>,
                                         CoefType, helpers::hash_pair>;
/**
 * @brief Container class for storage of the QUBO model.
 * 
 * @tparam NodeType integer number type used to index nodes of the QUBO model.
 * @tparam CoefType real type used to store values of the QUBO model.
 */
template <class NodeType, class CoefType> class QUBOModel {
protected:
 /**
  * @brief QUBO model linear coefficients store.
  * 
  */
  LinearCoef<NodeType, CoefType> linear;
  /**
   * @brief QUBO model quadratic coefficients store.
   * 
   */
  QuadraticCoef<NodeType, CoefType> quadratic;
  /**
   * @brief Number of nodes (variables) in the QUBO model.
   * 
   */
  ulong num_nodes = 0;

public:
  /**
   * @brief Construct a new QUBOModel object.
   * 
   * @param c_linear Linear coefficients of the QUBO model.
   * @param c_quadratic Quadratic coefficients of the QUBO model.
   */
  QUBOModel(LinearCoef<NodeType, CoefType> &c_linear,
            QuadraticCoef<NodeType, CoefType> &c_quadratic);
  
  /**
   * @brief Returns a textual representation of the QUBO model.
   * 
   * @return std::string Textual representation of the stored QUBO model.
   */
  std::string str() const;

  /**
   * @brief Add a binary variable (node) to the QUBO model.
   * 
   * @param vi Variable addres.
   * @param hi Value of the variable coefficient.
   */
  void add_variable(const NodeType &vi, const CoefType &hi);

  /**
   * @brief Add a connection (coupling) between binaray variables.
   * 
   * @param connection Addresses of the variables to connect (couple).
   * @param Ji Strength of the connection (coupling).
   */
  void add_connection(const std::pair<NodeType, NodeType> &connection,
                      const CoefType &Ji);

  /**
   * @brief Get the variable object.
   * 
   * @param vi Address of the variable.
   * @return const CoefType Variable coefficient.
   */
  const CoefType get_variable(const NodeType &vi) const;

  /**
   * @brief Get the connection object.
   * 
   * @param connection Addresses of the variables.
   * @return const CoefType Connection (coupling) strength between variables.
   */
  const CoefType
  get_connection(const std::pair<NodeType, NodeType> &connection) const;
  
  /**
   * @brief Set the number of nodes (variables) in the QUBO model.
   * 
   * @param num_nodes Number of nodes to set.
   */
  void set_nodes(int num_nodes);

  /**
   * @brief Get number of nodes (variables) in the QUBO model.
   * 
   * @return ulong Number of nodes.
   */
  ulong get_nodes() const;

  /**
   * @brief Print textual representation of the QUBO model to a stream.
   * 
   * @param os Ouput stream reference.
   * @param qubos QUBO model.
   * @return std::ostream& Modified ouput stream reference.
   */
  friend std::ostream &operator<<(std::ostream &os,
                                  qubo::QUBOModel<NodeType, CoefType> &qubos);
  
  /**
   * @brief Create a QUBO model from the input stream.
   * 
   * @param is Input stream reference.
   * @param qubos QUBO model.
   * @return std::istream& Modified input stream reference.
   */
  friend std::istream &operator>>(std::istream &is,
                                  qubo::QUBOModel<NodeType, CoefType> &qubos);
  
  /**
   * @brief Load a QUBO model for an input stream.
   * 
   * @param stream Input stream refernce.
   * @return QUBOModel<int, double> Constructed QUBO model object.
   */
  static QUBOModel<int, double> load(std::istream &stream);
};

template <class NodeType, class CoefType>
QUBOModel<NodeType, CoefType>::QUBOModel(
    qubo::LinearCoef<NodeType, CoefType> &c_linear,
    qubo::QuadraticCoef<NodeType, CoefType> &c_quadratic) {
  linear = c_linear;
  quadratic = c_quadratic;
}

template <class NodeType, class CoefType>
void QUBOModel<NodeType, CoefType>::set_nodes(int number) {
  QUBOModel::num_nodes = number;
}
template <class NodeType, class CoefType>
ulong QUBOModel<NodeType, CoefType>::get_nodes() const {
  return num_nodes;
}

template <class NodeType, class CoefType>
std::string QUBOModel<NodeType, CoefType>::str() const {

  std::string s;

  s.append("QUBO model");

  for (auto &qii : linear) {
    s.append(" " + std::to_string(qii.first) + "--" +
             std::to_string(qii.first) + ":" + std::to_string(qii.second) +
             " ");
  }

  for (auto &qij : quadratic) {
    s.append(std::to_string(qij.first.first) + "--" +
             std::to_string(qij.first.second) + ":" +
             std::to_string(qij.second));
  }

  return s;
}

template <class NodeType, class CoefType>
void QUBOModel<NodeType, CoefType>::add_variable(const NodeType &vi,
                                                 const CoefType &hi) {
  helpers::insert_model(linear, vi, hi);
}

template <class NodeType, class CoefType>
void QUBOModel<NodeType, CoefType>::add_connection(
    const std::pair<NodeType, NodeType> &connection, const CoefType &Ji) {
  helpers::insert_model(quadratic, connection, Ji);
}

template <class NodeType, class CoefType>
const CoefType QUBOModel<NodeType, CoefType>::get_connection(
    const std::pair<NodeType, NodeType> &connection) const {
  if (quadratic.count(connection) == 0) {

    return 0;

  } else {

    return quadratic.at(connection);
  }
}

template <class NodeType, class CoefType>
const CoefType
QUBOModel<NodeType, CoefType>::get_variable(const NodeType &vi) const {
  if (linear.count(vi) == 0) {

    return 0;

  } else {

    return linear.at(vi);
  }
}

/**
 * @brief Print textual representation of the QUBO model to a stream.
 * 
 * @param os Ouput stream reference.
 * @param qubos QUBO model.
 * @return std::ostream& Modified ouput stream reference.
 */
template <class NodeType, class CoefType>
std::ostream &operator<<(std::ostream &os,
                         const qubo::QUBOModel<NodeType, CoefType> &qubos) {
  os << qubos.str();
  return os;
}

/**
 * @brief Create a QUBO model from the input stream.
 * 
 * @param is Input stream reference.
 * @param qubos QUBO model.
 * @return std::istream& Modified input stream reference.
 */
template <class NodeType, class CoefType>
std::istream &operator>>(std::istream &is,
                         const qubo::QUBOModel<NodeType, CoefType> &qubos) {
  NodeType J1;
  NodeType J2;
  CoefType value;
  is >> J1 >> J2 >> value;
  if (J1 == J2) {
    qubos.add_variable(J1, value);
  } else {
    qubos.add_connection(std::make_pair(J1, J2), value);
  }

  return is;
}

  /*!
  * @brief Builder for QUBOModel.
  */
  struct QUBOBuilder {
  private:
    qubo::LinearCoef<int, double> linear_c{};
    qubo::QuadraticCoef<int, double> quadratic_c{};
    int num_quadratic = -1;
    int num_linear = -1;
    int max_nodes = -1;
    std::set<int> used_nodes;

  public:
    /*!
    * @brief Add new QUBO Model coefficient.
    *
    * Methods of this class are meant to be used in conjunction with
    * appropriate boost spirit based parser which should call them
    * as actions after parsing respective portions of the input file.
    *
    * @param element a triple (i, j, q_ij) defining coefficient to
    *                be added. The indices pair has to be ordered such
    *                that i <= j.
    * @throw std::invalid_argument if i > j.
    */
    void add_element(boost::fusion::vector<int, int, double> element) {
      auto i = at_c<0>(element);
      auto j = at_c<1>(element);
      auto coef = at_c<2>(element);
      if (i == j) {
        linear_c.insert({i, coef});
      } else if (i < j) {
        quadratic_c.insert({{i, j}, coef});
      } else {
        throw std::invalid_argument("Incorrect file, encountered coefficient "
                                    "from lower triangle of QUBO matrix.");
      }
      used_nodes.insert(i);
      used_nodes.insert(j);
    }

    /*!
    * @brief Set number of quadratic coefficients in the model.
    *
    * @param value number of quadratic coefficients.
    */
    void set_num_quadratic(int value) { num_quadratic = value; }

    /*!
    * @brief Set number of linear coefficients in the model.
    *
    * @param value number of linear coefficients.
    */
    void set_num_linear(int value) { num_linear = value; }

    /*!
    * @brief Set maximum number of nodes in the model.
    *
    * @param value maximm number of nodes in the model.
    */
    void set_max_nodes(int value) { max_nodes = value; }

    /*!
    * @brief Build QUBO Model based on information provided this far
    *        to the builder.
    *
    * @throw std::invalid_argument if one of the following conditions
    *        is met:
    *        - number of quadratic and/or linear coefficients is different
    *          from the declared one.
    *        - number of quadratic/linear coefficients or maximum number
    *          of nodes was not provided.
    *        - there are more variables than the declard maximum number.
    * @returns Constructed QUBO model.
    */
    qubo::QUBOModel<int, double> build_qubo() {
      if (linear_c.size() == 0 && quadratic_c.size() == 0) {
        throw std::invalid_argument("An empty input, no coefficients defined.");
      }
      if (num_quadratic == -1 || num_linear == -1 || max_nodes == -1) {
        throw std::invalid_argument(
            "No header line or header line misformatted.");
      }
      if (num_linear > max_nodes) {
        throw std::invalid_argument(
            "Number of linear terms is greater than num nodes.");
      }
      if (quadratic_c.size() != num_quadratic) {
        throw std::invalid_argument(
            "Number of quadratic terms is not equal to the declared one.");
      }
      if (linear_c.size() != num_linear) {
        throw std::invalid_argument(
            "Number of linear terms is not equal to the declared one.");
      }

      qubo::QUBOModel q(linear_c, quadratic_c);
      q.set_nodes(*used_nodes.rbegin() + 1);
      return q;
    }
  };

  /*!
  * @brief Parse QUBO using boost spirit based parser.
  *
  * @param first Iterator to the start of the input stream containing
  *              definition of QUBO.
  * @param last Iterator to the end of the input stream containing
  *              definition of QUBO.
  * @throw std::invalid_argument if parsing failed, either from syntactic
  *        or semantic reasons.
  * @returns Constructed QUBO.
  */
  template <typename Iterator>
  qubo::QUBOModel<int, double> parse_qubo(Iterator first, Iterator last) {
    QUBOBuilder builder;
    auto real = qi::lexeme[double_];
    auto index = qi::lexeme[qi::uint_ >> " "];
    auto coefficients =
        (index >> index >>
        real)[boost::bind(&QUBOBuilder::add_element, &builder, _1)] >>
        (qi::eol | qi::eoi);
    auto comment = qi::lit('c') >> *(qi::print) >> (qi::eol | qi::eoi);
    auto header =
        (qi::lit("p qubo") >> qi::uint_ >>
        qi::uint_[boost::bind(&QUBOBuilder::set_max_nodes, &builder, _1)] >>
        qi::uint_[boost::bind(&QUBOBuilder::set_num_linear, &builder, _1)] >>
        qi::uint_[boost::bind(&QUBOBuilder::set_num_quadratic, &builder, _1)] >>
        (qi::eol | qi::eoi));

    auto result = phrase_parse(first, last,
                              (*comment >> -header >> *(comment | coefficients)),
                              qi::blank);

    if (!result || first != last) {
      throw std::invalid_argument("Parsing failed. Incorrect file format.");
    }
    return builder.build_qubo();
  }

  /*!
  * @brief Load QUBO from an input stream.
  *
  * @param stream an input stream containing definition of QUBO.
  * @returns A loaded instance.
  */
  template <class NodeType, class CoefType>
  qubo::QUBOModel<int, double>
  qubo::QUBOModel<NodeType, CoefType>::load(std::istream &stream) {
    boost::spirit::istream_iterator first(stream), last;
    return parse_qubo(first, last);
  }

}; // namespace qubo
#endif