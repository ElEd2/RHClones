// Generated by rstantools.  Do not edit by hand.

/*
    RHClones is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RHClones is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RHClones.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_pop_model_normal_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_pop_model_normal");
    reader.add_event(60, 58, "end", "model_pop_model_normal");
    return reader;
}
#include <stan_meta_header.hpp>
class model_pop_model_normal
  : public stan::model::model_base_crtp<model_pop_model_normal> {
private:
        int N;
        vector_d age;
        std::vector<int> counts;
        std::vector<int> num_crypts;
public:
    model_pop_model_normal(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_pop_model_normal(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_pop_model_normal_namespace::model_pop_model_normal";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 3;
            validate_non_negative_index("age", "N", N);
            context__.validate_dims("data initialization", "age", "vector_d", context__.to_vec(N));
            age = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("age");
            pos__ = 0;
            size_t age_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < age_j_1_max__; ++j_1__) {
                age(j_1__) = vals_r__[pos__++];
            }
            check_greater_or_equal(function__, "age", age, 0);
            current_statement_begin__ = 4;
            validate_non_negative_index("counts", "N", N);
            context__.validate_dims("data initialization", "counts", "int", context__.to_vec(N));
            counts = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("counts");
            pos__ = 0;
            size_t counts_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < counts_k_0_max__; ++k_0__) {
                counts[k_0__] = vals_i__[pos__++];
            }
            size_t counts_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < counts_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "counts[i_0__]", counts[i_0__], 0);
            }
            current_statement_begin__ = 5;
            validate_non_negative_index("num_crypts", "N", N);
            context__.validate_dims("data initialization", "num_crypts", "int", context__.to_vec(N));
            num_crypts = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("num_crypts");
            pos__ = 0;
            size_t num_crypts_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < num_crypts_k_0_max__; ++k_0__) {
                num_crypts[k_0__] = vals_i__[pos__++];
            }
            size_t num_crypts_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < num_crypts_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "num_crypts[i_0__]", num_crypts[i_0__], 0);
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 9;
            validate_non_negative_index("a_i", "N", N);
            num_params_r__ += N;
            current_statement_begin__ = 12;
            num_params_r__ += 1;
            current_statement_begin__ = 13;
            num_params_r__ += 1;
            current_statement_begin__ = 14;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_pop_model_normal() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 9;
        if (!(context__.contains_r("a_i")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable a_i missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("a_i");
        pos__ = 0U;
        validate_non_negative_index("a_i", "N", N);
        context__.validate_dims("parameter initialization", "a_i", "vector_d", context__.to_vec(N));
        Eigen::Matrix<double, Eigen::Dynamic, 1> a_i(N);
        size_t a_i_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < a_i_j_1_max__; ++j_1__) {
            a_i(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lb_unconstrain(0, a_i);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable a_i: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 12;
        if (!(context__.contains_r("x_intercept")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable x_intercept missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("x_intercept");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "x_intercept", "double", context__.to_vec());
        double x_intercept(0);
        x_intercept = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(x_intercept);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable x_intercept: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 13;
        if (!(context__.contains_r("pop_mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable pop_mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("pop_mu");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "pop_mu", "double", context__.to_vec());
        double pop_mu(0);
        pop_mu = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, pop_mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable pop_mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 14;
        if (!(context__.contains_r("pop_sd")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable pop_sd missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("pop_sd");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "pop_sd", "double", context__.to_vec());
        double pop_sd(0);
        pop_sd = vals_r__[pos__++];
        try {
            writer__.scalar_lb_unconstrain(0, pop_sd);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable pop_sd: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 9;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> a_i;
            (void) a_i;  // dummy to suppress unused var warning
            if (jacobian__)
                a_i = in__.vector_lb_constrain(0, N, lp__);
            else
                a_i = in__.vector_lb_constrain(0, N);
            current_statement_begin__ = 12;
            local_scalar_t__ x_intercept;
            (void) x_intercept;  // dummy to suppress unused var warning
            if (jacobian__)
                x_intercept = in__.scalar_constrain(lp__);
            else
                x_intercept = in__.scalar_constrain();
            current_statement_begin__ = 13;
            local_scalar_t__ pop_mu;
            (void) pop_mu;  // dummy to suppress unused var warning
            if (jacobian__)
                pop_mu = in__.scalar_lb_constrain(0, lp__);
            else
                pop_mu = in__.scalar_lb_constrain(0);
            current_statement_begin__ = 14;
            local_scalar_t__ pop_sd;
            (void) pop_sd;  // dummy to suppress unused var warning
            if (jacobian__)
                pop_sd = in__.scalar_lb_constrain(0, lp__);
            else
                pop_sd = in__.scalar_lb_constrain(0);
            // transformed parameters
            current_statement_begin__ = 19;
            validate_non_negative_index("p_i", "N", N);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> p_i(N);
            stan::math::initialize(p_i, DUMMY_VAR__);
            stan::math::fill(p_i, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 20;
            stan::math::assign(p_i, elt_multiply(a_i, subtract(age, x_intercept)));
            current_statement_begin__ = 21;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 22;
                if (as_bool(logical_gt(get_base1(p_i, i, "p_i", 1), 1))) {
                    current_statement_begin__ = 22;
                    stan::model::assign(p_i, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                (1 - 1e-16), 
                                "assigning variable p_i");
                }
                current_statement_begin__ = 23;
                if (as_bool(logical_lt(get_base1(p_i, i, "p_i", 1), 0))) {
                    current_statement_begin__ = 23;
                    stan::model::assign(p_i, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                1e-16, 
                                "assigning variable p_i");
                }
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 19;
            size_t p_i_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < p_i_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(p_i(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: p_i" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable p_i: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            // model body
            current_statement_begin__ = 29;
            lp_accum__.add(normal_log<propto__>(x_intercept, 0, 10));
            current_statement_begin__ = 30;
            lp_accum__.add(gamma_log<propto__>(pop_mu, 1e-2, 1e-2));
            current_statement_begin__ = 31;
            lp_accum__.add(gamma_log<propto__>(pop_sd, 1e-2, 1e-2));
            current_statement_begin__ = 34;
            lp_accum__.add(normal_log<propto__>(a_i, pop_mu, pop_sd));
            current_statement_begin__ = 35;
            lp_accum__.add(binomial_log<propto__>(counts, num_crypts, p_i));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("a_i");
        names__.push_back("x_intercept");
        names__.push_back("pop_mu");
        names__.push_back("pop_sd");
        names__.push_back("p_i");
        names__.push_back("y_pred");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_pop_model_normal_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> a_i = in__.vector_lb_constrain(0, N);
        size_t a_i_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < a_i_j_1_max__; ++j_1__) {
            vars__.push_back(a_i(j_1__));
        }
        double x_intercept = in__.scalar_constrain();
        vars__.push_back(x_intercept);
        double pop_mu = in__.scalar_lb_constrain(0);
        vars__.push_back(pop_mu);
        double pop_sd = in__.scalar_lb_constrain(0);
        vars__.push_back(pop_sd);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 19;
            validate_non_negative_index("p_i", "N", N);
            Eigen::Matrix<double, Eigen::Dynamic, 1> p_i(N);
            stan::math::initialize(p_i, DUMMY_VAR__);
            stan::math::fill(p_i, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 20;
            stan::math::assign(p_i, elt_multiply(a_i, subtract(age, x_intercept)));
            current_statement_begin__ = 21;
            for (int i = 1; i <= N; ++i) {
                current_statement_begin__ = 22;
                if (as_bool(logical_gt(get_base1(p_i, i, "p_i", 1), 1))) {
                    current_statement_begin__ = 22;
                    stan::model::assign(p_i, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                (1 - 1e-16), 
                                "assigning variable p_i");
                }
                current_statement_begin__ = 23;
                if (as_bool(logical_lt(get_base1(p_i, i, "p_i", 1), 0))) {
                    current_statement_begin__ = 23;
                    stan::model::assign(p_i, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                1e-16, 
                                "assigning variable p_i");
                }
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            // write transformed parameters
            if (include_tparams__) {
                size_t p_i_j_1_max__ = N;
                for (size_t j_1__ = 0; j_1__ < p_i_j_1_max__; ++j_1__) {
                    vars__.push_back(p_i(j_1__));
                }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 39;
            validate_non_negative_index("y_pred", "N", N);
            std::vector<int> y_pred(N, int(0));
            stan::math::fill(y_pred, std::numeric_limits<int>::min());
            current_statement_begin__ = 40;
            validate_non_negative_index("log_lik", "N", N);
            Eigen::Matrix<double, Eigen::Dynamic, 1> log_lik(N);
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 42;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 44;
                stan::model::assign(log_lik, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            binomial_log(get_base1(counts, n, "counts", 1), get_base1(num_crypts, n, "num_crypts", 1), get_base1(p_i, n, "p_i", 1)), 
                            "assigning variable log_lik");
            }
            current_statement_begin__ = 47;
            for (int n = 1; n <= N; ++n) {
                current_statement_begin__ = 49;
                stan::model::assign(y_pred, 
                            stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                            binomial_rng(get_base1(num_crypts, n, "num_crypts", 1), get_base1(p_i, n, "p_i", 1), base_rng__), 
                            "assigning variable y_pred");
            }
            // validate, write generated quantities
            current_statement_begin__ = 39;
            size_t y_pred_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < y_pred_k_0_max__; ++k_0__) {
                vars__.push_back(y_pred[k_0__]);
            }
            current_statement_begin__ = 40;
            size_t log_lik_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
                vars__.push_back(log_lik(j_1__));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_pop_model_normal";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t a_i_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < a_i_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "a_i" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "x_intercept";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pop_mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pop_sd";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t p_i_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < p_i_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "p_i" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        size_t y_pred_k_0_max__ = N;
        for (size_t k_0__ = 0; k_0__ < y_pred_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_pred" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t a_i_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < a_i_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "a_i" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "x_intercept";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pop_mu";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "pop_sd";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t p_i_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < p_i_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "p_i" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
        size_t y_pred_k_0_max__ = N;
        for (size_t k_0__ = 0; k_0__ < y_pred_k_0_max__; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "y_pred" << '.' << k_0__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t log_lik_j_1_max__ = N;
        for (size_t j_1__ = 0; j_1__ < log_lik_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_pop_model_normal_namespace::model_pop_model_normal stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif