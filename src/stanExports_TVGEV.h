// Generated by rstantools.  Do not edit by hand.

/*
    beverstan is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    beverstan is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beverstan.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_TVGEV_namespace {
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
    reader.add_event(0, 0, "start", "model_TVGEV");
    reader.add_event(117, 115, "end", "model_TVGEV");
    return reader;
}
template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& mu,
             const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& sigma,
             const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& xi, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 12;
        validate_non_negative_index("t", "rows(y)", rows(y));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> t(rows(y));
        stan::math::initialize(t, DUMMY_VAR__);
        stan::math::fill(t, DUMMY_VAR__);
        current_statement_begin__ = 13;
        validate_non_negative_index("lp", "rows(y)", rows(y));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lp(rows(y));
        stan::math::initialize(lp, DUMMY_VAR__);
        stan::math::fill(lp, DUMMY_VAR__);
        current_statement_begin__ = 14;
        int n(0);
        (void) n;  // dummy to suppress unused var warning
        stan::math::fill(n, std::numeric_limits<int>::min());
        current_statement_begin__ = 15;
        stan::math::assign(n, rows(y));
        current_statement_begin__ = 16;
        for (int i = 1; i <= n; ++i) {
            current_statement_begin__ = 17;
            stan::model::assign(t, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        (logical_eq(get_base1(xi, i, "xi", 1), 0) ? stan::math::exp(((get_base1(mu, i, "mu", 1) - get_base1(y, i, "y", 1)) / get_base1(sigma, i, "sigma", 1))) : pow((1.0 + (get_base1(xi, i, "xi", 1) * ((get_base1(y, i, "y", 1) - get_base1(mu, i, "mu", 1)) / get_base1(sigma, i, "sigma", 1)))), (-(1.0) / get_base1(xi, i, "xi", 1))) ), 
                        "assigning variable t");
            current_statement_begin__ = 19;
            stan::model::assign(lp, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        ((-(stan::math::log(get_base1(sigma, i, "sigma", 1))) + ((get_base1(xi, i, "xi", 1) + 1.0) * stan::math::log(get_base1(t, i, "t", 1)))) - get_base1(t, i, "t", 1)), 
                        "assigning variable lp");
        }
        current_statement_begin__ = 21;
        return stan::math::promote_scalar<fun_return_scalar_t__>(sum(lp));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& mu,
             const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& sigma,
             const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& xi, std::ostream* pstream__) {
    return gev_lpdf<false>(y,mu,sigma,xi, pstream__);
}
struct gev_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& mu,
             const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& sigma,
             const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& xi, std::ostream* pstream__) const {
        return gev_lpdf(y, mu, sigma, xi, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
gev_lcdf(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& mu,
             const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& sigma,
             const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& xi, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 25;
        validate_non_negative_index("t", "rows(y)", rows(y));
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> t(rows(y));
        stan::math::initialize(t, DUMMY_VAR__);
        stan::math::fill(t, DUMMY_VAR__);
        current_statement_begin__ = 26;
        int n(0);
        (void) n;  // dummy to suppress unused var warning
        stan::math::fill(n, std::numeric_limits<int>::min());
        current_statement_begin__ = 27;
        stan::math::assign(n, rows(y));
        current_statement_begin__ = 28;
        for (int i = 1; i <= n; ++i) {
            current_statement_begin__ = 29;
            stan::model::assign(t, 
                        stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                        (logical_eq(get_base1(xi, i, "xi", 1), 0) ? stan::math::exp(((get_base1(mu, i, "mu", 1) - get_base1(y, i, "y", 1)) / get_base1(sigma, i, "sigma", 1))) : pow((1.0 + (get_base1(xi, i, "xi", 1) * ((get_base1(y, i, "y", 1) - get_base1(mu, i, "mu", 1)) / get_base1(sigma, i, "sigma", 1)))), (-(1.0) / get_base1(xi, i, "xi", 1))) ), 
                        "assigning variable t");
        }
        current_statement_begin__ = 32;
        return stan::math::promote_scalar<fun_return_scalar_t__>(-(sum(t)));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct gev_lcdf_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& y,
             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& mu,
             const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& sigma,
             const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& xi, std::ostream* pstream__) const {
        return gev_lcdf(y, mu, sigma, xi, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_TVGEV
  : public stan::model::model_base_crtp<model_TVGEV> {
private:
        int n;
        vector_d y;
        int p_mu;
        int p_sigma;
        int p_xi;
        int cst_mu;
        int cst_sigma;
        int cst_xi;
        matrix_d X_mu;
        matrix_d X_sigma;
        matrix_d X_xi;
        matrix_d cov_psi_mu0;
        matrix_d cov_psi_sigma0;
        matrix_d cov_psi_xi0;
        vector_d mean_psi_mu0;
        vector_d mean_psi_sigma0;
        vector_d mean_psi_xi0;
public:
    model_TVGEV(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_TVGEV(stan::io::var_context& context__,
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
        static const char* function__ = "model_TVGEV_namespace::model_TVGEV";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 39;
            context__.validate_dims("data initialization", "n", "int", context__.to_vec());
            n = int(0);
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            n = vals_i__[pos__++];
            check_greater_or_equal(function__, "n", n, 0);
            current_statement_begin__ = 40;
            validate_non_negative_index("y", "n", n);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(n));
            y = Eigen::Matrix<double, Eigen::Dynamic, 1>(n);
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = n;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                y(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 43;
            context__.validate_dims("data initialization", "p_mu", "int", context__.to_vec());
            p_mu = int(0);
            vals_i__ = context__.vals_i("p_mu");
            pos__ = 0;
            p_mu = vals_i__[pos__++];
            check_greater_or_equal(function__, "p_mu", p_mu, 1);
            current_statement_begin__ = 44;
            context__.validate_dims("data initialization", "p_sigma", "int", context__.to_vec());
            p_sigma = int(0);
            vals_i__ = context__.vals_i("p_sigma");
            pos__ = 0;
            p_sigma = vals_i__[pos__++];
            check_greater_or_equal(function__, "p_sigma", p_sigma, 1);
            current_statement_begin__ = 45;
            context__.validate_dims("data initialization", "p_xi", "int", context__.to_vec());
            p_xi = int(0);
            vals_i__ = context__.vals_i("p_xi");
            pos__ = 0;
            p_xi = vals_i__[pos__++];
            check_greater_or_equal(function__, "p_xi", p_xi, 1);
            current_statement_begin__ = 48;
            context__.validate_dims("data initialization", "cst_mu", "int", context__.to_vec());
            cst_mu = int(0);
            vals_i__ = context__.vals_i("cst_mu");
            pos__ = 0;
            cst_mu = vals_i__[pos__++];
            check_greater_or_equal(function__, "cst_mu", cst_mu, 0);
            current_statement_begin__ = 49;
            context__.validate_dims("data initialization", "cst_sigma", "int", context__.to_vec());
            cst_sigma = int(0);
            vals_i__ = context__.vals_i("cst_sigma");
            pos__ = 0;
            cst_sigma = vals_i__[pos__++];
            check_greater_or_equal(function__, "cst_sigma", cst_sigma, 0);
            current_statement_begin__ = 50;
            context__.validate_dims("data initialization", "cst_xi", "int", context__.to_vec());
            cst_xi = int(0);
            vals_i__ = context__.vals_i("cst_xi");
            pos__ = 0;
            cst_xi = vals_i__[pos__++];
            check_greater_or_equal(function__, "cst_xi", cst_xi, 0);
            current_statement_begin__ = 53;
            validate_non_negative_index("X_mu", "n", n);
            validate_non_negative_index("X_mu", "p_mu", p_mu);
            context__.validate_dims("data initialization", "X_mu", "matrix_d", context__.to_vec(n,p_mu));
            X_mu = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, p_mu);
            vals_r__ = context__.vals_r("X_mu");
            pos__ = 0;
            size_t X_mu_j_2_max__ = p_mu;
            size_t X_mu_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_mu_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_mu_j_1_max__; ++j_1__) {
                    X_mu(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 54;
            validate_non_negative_index("X_sigma", "n", n);
            validate_non_negative_index("X_sigma", "p_sigma", p_sigma);
            context__.validate_dims("data initialization", "X_sigma", "matrix_d", context__.to_vec(n,p_sigma));
            X_sigma = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, p_sigma);
            vals_r__ = context__.vals_r("X_sigma");
            pos__ = 0;
            size_t X_sigma_j_2_max__ = p_sigma;
            size_t X_sigma_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_sigma_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_sigma_j_1_max__; ++j_1__) {
                    X_sigma(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 55;
            validate_non_negative_index("X_xi", "n", n);
            validate_non_negative_index("X_xi", "p_xi", p_xi);
            context__.validate_dims("data initialization", "X_xi", "matrix_d", context__.to_vec(n,p_xi));
            X_xi = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(n, p_xi);
            vals_r__ = context__.vals_r("X_xi");
            pos__ = 0;
            size_t X_xi_j_2_max__ = p_xi;
            size_t X_xi_j_1_max__ = n;
            for (size_t j_2__ = 0; j_2__ < X_xi_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_xi_j_1_max__; ++j_1__) {
                    X_xi(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 58;
            validate_non_negative_index("cov_psi_mu0", "p_mu", p_mu);
            validate_non_negative_index("cov_psi_mu0", "p_mu", p_mu);
            context__.validate_dims("data initialization", "cov_psi_mu0", "matrix_d", context__.to_vec(p_mu,p_mu));
            cov_psi_mu0 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(p_mu, p_mu);
            vals_r__ = context__.vals_r("cov_psi_mu0");
            pos__ = 0;
            size_t cov_psi_mu0_j_2_max__ = p_mu;
            size_t cov_psi_mu0_j_1_max__ = p_mu;
            for (size_t j_2__ = 0; j_2__ < cov_psi_mu0_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < cov_psi_mu0_j_1_max__; ++j_1__) {
                    cov_psi_mu0(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 59;
            validate_non_negative_index("cov_psi_sigma0", "p_sigma", p_sigma);
            validate_non_negative_index("cov_psi_sigma0", "p_sigma", p_sigma);
            context__.validate_dims("data initialization", "cov_psi_sigma0", "matrix_d", context__.to_vec(p_sigma,p_sigma));
            cov_psi_sigma0 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(p_sigma, p_sigma);
            vals_r__ = context__.vals_r("cov_psi_sigma0");
            pos__ = 0;
            size_t cov_psi_sigma0_j_2_max__ = p_sigma;
            size_t cov_psi_sigma0_j_1_max__ = p_sigma;
            for (size_t j_2__ = 0; j_2__ < cov_psi_sigma0_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < cov_psi_sigma0_j_1_max__; ++j_1__) {
                    cov_psi_sigma0(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 60;
            validate_non_negative_index("cov_psi_xi0", "p_xi", p_xi);
            validate_non_negative_index("cov_psi_xi0", "p_xi", p_xi);
            context__.validate_dims("data initialization", "cov_psi_xi0", "matrix_d", context__.to_vec(p_xi,p_xi));
            cov_psi_xi0 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(p_xi, p_xi);
            vals_r__ = context__.vals_r("cov_psi_xi0");
            pos__ = 0;
            size_t cov_psi_xi0_j_2_max__ = p_xi;
            size_t cov_psi_xi0_j_1_max__ = p_xi;
            for (size_t j_2__ = 0; j_2__ < cov_psi_xi0_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < cov_psi_xi0_j_1_max__; ++j_1__) {
                    cov_psi_xi0(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 63;
            validate_non_negative_index("mean_psi_mu0", "p_mu", p_mu);
            context__.validate_dims("data initialization", "mean_psi_mu0", "vector_d", context__.to_vec(p_mu));
            mean_psi_mu0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(p_mu);
            vals_r__ = context__.vals_r("mean_psi_mu0");
            pos__ = 0;
            size_t mean_psi_mu0_j_1_max__ = p_mu;
            for (size_t j_1__ = 0; j_1__ < mean_psi_mu0_j_1_max__; ++j_1__) {
                mean_psi_mu0(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 64;
            validate_non_negative_index("mean_psi_sigma0", "p_sigma", p_sigma);
            context__.validate_dims("data initialization", "mean_psi_sigma0", "vector_d", context__.to_vec(p_sigma));
            mean_psi_sigma0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(p_sigma);
            vals_r__ = context__.vals_r("mean_psi_sigma0");
            pos__ = 0;
            size_t mean_psi_sigma0_j_1_max__ = p_sigma;
            for (size_t j_1__ = 0; j_1__ < mean_psi_sigma0_j_1_max__; ++j_1__) {
                mean_psi_sigma0(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 65;
            validate_non_negative_index("mean_psi_xi0", "p_xi", p_xi);
            context__.validate_dims("data initialization", "mean_psi_xi0", "vector_d", context__.to_vec(p_xi));
            mean_psi_xi0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(p_xi);
            vals_r__ = context__.vals_r("mean_psi_xi0");
            pos__ = 0;
            size_t mean_psi_xi0_j_1_max__ = p_xi;
            for (size_t j_1__ = 0; j_1__ < mean_psi_xi0_j_1_max__; ++j_1__) {
                mean_psi_xi0(j_1__) = vals_r__[pos__++];
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 71;
            validate_non_negative_index("psi_mu", "p_mu", p_mu);
            num_params_r__ += p_mu;
            current_statement_begin__ = 72;
            validate_non_negative_index("psi_sigma", "p_sigma", p_sigma);
            num_params_r__ += p_sigma;
            current_statement_begin__ = 73;
            validate_non_negative_index("psi_xi", "p_xi", p_xi);
            num_params_r__ += p_xi;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_TVGEV() { }
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
        current_statement_begin__ = 71;
        if (!(context__.contains_r("psi_mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable psi_mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("psi_mu");
        pos__ = 0U;
        validate_non_negative_index("psi_mu", "p_mu", p_mu);
        context__.validate_dims("parameter initialization", "psi_mu", "vector_d", context__.to_vec(p_mu));
        Eigen::Matrix<double, Eigen::Dynamic, 1> psi_mu(p_mu);
        size_t psi_mu_j_1_max__ = p_mu;
        for (size_t j_1__ = 0; j_1__ < psi_mu_j_1_max__; ++j_1__) {
            psi_mu(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(psi_mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable psi_mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 72;
        if (!(context__.contains_r("psi_sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable psi_sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("psi_sigma");
        pos__ = 0U;
        validate_non_negative_index("psi_sigma", "p_sigma", p_sigma);
        context__.validate_dims("parameter initialization", "psi_sigma", "vector_d", context__.to_vec(p_sigma));
        Eigen::Matrix<double, Eigen::Dynamic, 1> psi_sigma(p_sigma);
        size_t psi_sigma_j_1_max__ = p_sigma;
        for (size_t j_1__ = 0; j_1__ < psi_sigma_j_1_max__; ++j_1__) {
            psi_sigma(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(psi_sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable psi_sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 73;
        if (!(context__.contains_r("psi_xi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable psi_xi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("psi_xi");
        pos__ = 0U;
        validate_non_negative_index("psi_xi", "p_xi", p_xi);
        context__.validate_dims("parameter initialization", "psi_xi", "vector_d", context__.to_vec(p_xi));
        Eigen::Matrix<double, Eigen::Dynamic, 1> psi_xi(p_xi);
        size_t psi_xi_j_1_max__ = p_xi;
        for (size_t j_1__ = 0; j_1__ < psi_xi_j_1_max__; ++j_1__) {
            psi_xi(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(psi_xi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable psi_xi: ") + e.what()), current_statement_begin__, prog_reader__());
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
            current_statement_begin__ = 71;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> psi_mu;
            (void) psi_mu;  // dummy to suppress unused var warning
            if (jacobian__)
                psi_mu = in__.vector_constrain(p_mu, lp__);
            else
                psi_mu = in__.vector_constrain(p_mu);
            current_statement_begin__ = 72;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> psi_sigma;
            (void) psi_sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                psi_sigma = in__.vector_constrain(p_sigma, lp__);
            else
                psi_sigma = in__.vector_constrain(p_sigma);
            current_statement_begin__ = 73;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> psi_xi;
            (void) psi_xi;  // dummy to suppress unused var warning
            if (jacobian__)
                psi_xi = in__.vector_constrain(p_xi, lp__);
            else
                psi_xi = in__.vector_constrain(p_xi);
            // model body
            {
            current_statement_begin__ = 79;
            validate_non_negative_index("mu", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu(n);
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 80;
            validate_non_negative_index("sigma", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> sigma(n);
            stan::math::initialize(sigma, DUMMY_VAR__);
            stan::math::fill(sigma, DUMMY_VAR__);
            current_statement_begin__ = 81;
            validate_non_negative_index("xi", "n", n);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> xi(n);
            stan::math::initialize(xi, DUMMY_VAR__);
            stan::math::fill(xi, DUMMY_VAR__);
            current_statement_begin__ = 85;
            lp_accum__.add(multi_normal_log<propto__>(psi_mu, mean_psi_mu0, cov_psi_mu0));
            current_statement_begin__ = 86;
            if (as_bool(cst_mu)) {
                current_statement_begin__ = 87;
                for (int i = 1; i <= n; ++i) {
                    current_statement_begin__ = 88;
                    stan::model::assign(mu, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                get_base1(psi_mu, 1, "psi_mu", 1), 
                                "assigning variable mu");
                }
            } else {
                current_statement_begin__ = 91;
                stan::math::assign(mu, multiply(X_mu, psi_mu));
            }
            current_statement_begin__ = 94;
            lp_accum__.add(multi_normal_log<propto__>(psi_sigma, mean_psi_sigma0, cov_psi_sigma0));
            current_statement_begin__ = 95;
            if (as_bool(cst_sigma)) {
                current_statement_begin__ = 96;
                for (int i = 1; i <= n; ++i) {
                    current_statement_begin__ = 97;
                    stan::model::assign(sigma, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                get_base1(psi_sigma, 1, "psi_sigma", 1), 
                                "assigning variable sigma");
                }
            } else {
                current_statement_begin__ = 100;
                stan::math::assign(sigma, multiply(X_sigma, psi_sigma));
            }
            current_statement_begin__ = 103;
            lp_accum__.add(multi_normal_log<propto__>(psi_xi, mean_psi_xi0, cov_psi_xi0));
            current_statement_begin__ = 104;
            if (as_bool(cst_xi)) {
                current_statement_begin__ = 105;
                for (int i = 1; i <= n; ++i) {
                    current_statement_begin__ = 106;
                    stan::model::assign(xi, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list()), 
                                get_base1(psi_xi, 1, "psi_xi", 1), 
                                "assigning variable xi");
                }
            } else {
                current_statement_begin__ = 109;
                stan::math::assign(xi, multiply(X_xi, psi_xi));
            }
            current_statement_begin__ = 113;
            lp_accum__.add(gev_lpdf<propto__>(y, mu, sigma, xi, pstream__));
            }
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
        names__.push_back("psi_mu");
        names__.push_back("psi_sigma");
        names__.push_back("psi_xi");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(p_mu);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p_sigma);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(p_xi);
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
        static const char* function__ = "model_TVGEV_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> psi_mu = in__.vector_constrain(p_mu);
        size_t psi_mu_j_1_max__ = p_mu;
        for (size_t j_1__ = 0; j_1__ < psi_mu_j_1_max__; ++j_1__) {
            vars__.push_back(psi_mu(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> psi_sigma = in__.vector_constrain(p_sigma);
        size_t psi_sigma_j_1_max__ = p_sigma;
        for (size_t j_1__ = 0; j_1__ < psi_sigma_j_1_max__; ++j_1__) {
            vars__.push_back(psi_sigma(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> psi_xi = in__.vector_constrain(p_xi);
        size_t psi_xi_j_1_max__ = p_xi;
        for (size_t j_1__ = 0; j_1__ < psi_xi_j_1_max__; ++j_1__) {
            vars__.push_back(psi_xi(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
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
        return "model_TVGEV";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t psi_mu_j_1_max__ = p_mu;
        for (size_t j_1__ = 0; j_1__ < psi_mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi_mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t psi_sigma_j_1_max__ = p_sigma;
        for (size_t j_1__ = 0; j_1__ < psi_sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi_sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t psi_xi_j_1_max__ = p_xi;
        for (size_t j_1__ = 0; j_1__ < psi_xi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi_xi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t psi_mu_j_1_max__ = p_mu;
        for (size_t j_1__ = 0; j_1__ < psi_mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi_mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t psi_sigma_j_1_max__ = p_sigma;
        for (size_t j_1__ = 0; j_1__ < psi_sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi_sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t psi_xi_j_1_max__ = p_xi;
        for (size_t j_1__ = 0; j_1__ < psi_xi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi_xi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_TVGEV_namespace::model_TVGEV stan_model;
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
