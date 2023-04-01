licenses(["notice"])  # New BSD, portions MIT.

CERES_SRCS = ["internal/ceres/" + filename for filename in [
    "array_utils.cc",
    "block_evaluate_preparer.cc",
    "block_jacobi_preconditioner.cc",
    "block_jacobian_writer.cc",
    "block_random_access_dense_matrix.cc",
    "block_random_access_diagonal_matrix.cc",
    "block_random_access_matrix.cc",
    "block_random_access_sparse_matrix.cc",
    "block_sparse_matrix.cc",
    "block_structure.cc",
    "c_api.cc",
    "callbacks.cc",
    "canonical_views_clustering.cc",
    "cgnr_solver.cc",
    "compressed_col_sparse_matrix_utils.cc",
    "compressed_row_jacobian_writer.cc",
    "compressed_row_sparse_matrix.cc",
    "conditioned_cost_function.cc",
    "conjugate_gradients_solver.cc",
    "context.cc",
    "context_impl.cc",
    "coordinate_descent_minimizer.cc",
    "corrector.cc",
    "cost_function.cc",
    "covariance.cc",
    "covariance_impl.cc",
    "cxsparse.cc",
    "dense_cholesky.cc",
    "dense_normal_cholesky_solver.cc",
    "dense_qr.cc",
    "dense_qr_solver.cc",
    "dense_sparse_matrix.cc",
    "detect_structure.cc",
    "dogleg_strategy.cc",
    "dynamic_compressed_row_jacobian_writer.cc",
    "dynamic_compressed_row_sparse_matrix.cc",
    "dynamic_sparse_normal_cholesky_solver.cc",
    "eigensparse.cc",
    "evaluation_callback.cc",
    "evaluator.cc",
    "file.cc",
    "first_order_function.cc",
    "float_cxsparse.cc",
    "function_sample.cc",
    "gradient_checker.cc",
    "gradient_checking_cost_function.cc",
    "gradient_problem.cc",
    "gradient_problem_solver.cc",
    "implicit_schur_complement.cc",
    "inner_product_computer.cc",
    "is_close.cc",
    "iteration_callback.cc",
    "iterative_refiner.cc",
    "iterative_schur_complement_solver.cc",
    "levenberg_marquardt_strategy.cc",
    "line_search.cc",
    "line_search_direction.cc",
    "line_search_minimizer.cc",
    "line_search_preprocessor.cc",
    "linear_least_squares_problems.cc",
    "linear_operator.cc",
    "linear_solver.cc",
    "local_parameterization.cc",
    "loss_function.cc",
    "low_rank_inverse_hessian.cc",
    "manifold.cc",
    "minimizer.cc",
    "normal_prior.cc",
    "parallel_for_cxx.cc",
    "parallel_for_nothreads.cc",
    "parallel_for_openmp.cc",
    "parallel_utils.cc",
    "parameter_block_ordering.cc",
    "partitioned_matrix_view.cc",
    "polynomial.cc",
    "preconditioner.cc",
    "preprocessor.cc",
    "problem.cc",
    "problem_impl.cc",
    "program.cc",
    "reorder_program.cc",
    "residual_block.cc",
    "residual_block_utils.cc",
    "schur_complement_solver.cc",
    "schur_eliminator.cc",
    "schur_jacobi_preconditioner.cc",
    "schur_templates.cc",
    "scratch_evaluate_preparer.cc",
    "single_linkage_clustering.cc",
    "solver.cc",
    "solver_utils.cc",
    "sparse_cholesky.cc",
    "sparse_matrix.cc",
    "sparse_normal_cholesky_solver.cc",
    "stringprintf.cc",
    "subset_preconditioner.cc",
    "thread_pool.cc",
    "thread_token_provider.cc",
    "triplet_sparse_matrix.cc",
    "trust_region_minimizer.cc",
    "trust_region_preprocessor.cc",
    "trust_region_step_evaluator.cc",
    "trust_region_strategy.cc",
    "types.cc",
    "visibility_based_preconditioner.cc",
    "visibility.cc",
    "wall_time.cc",
]]

CERES_DEFINES = [
    "CERES_NO_SUITESPARSE",
    "CERES_NO_CXSPARSE",
    "CERES_NO_ACCELERATE_SPARSE",
    "CERES_NO_LAPACK",
    "CERES_USE_EIGEN_SPARSE",
    "CERES_USE_CXX_THREADS",
    "CERES_NO_CUDA",
    "CERES_EXPORT=",
    "CERES_NO_EXPORT=",
]

cc_library(
    name = "ceres",
    srcs = CERES_SRCS + glob([
        "internal/ceres/generated/schur_eliminator_*.cc",
        "internal/ceres/generated/partitioned_matrix_view_*.cc",
        "config/**/*.h",
        "internal/**/*.h",
    ]),
    hdrs = glob([
        "include/ceres/*.h",
        "config/ceres/internal/*.h",
        "include/ceres/internal/*.h",
    ]),
    copts = [
        "-fopenmp",
        "-Wno-deprecated-declarations",
        "-Wno-sign-compare",
    ],
    defines = CERES_DEFINES,
    includes = [
        "config",
        "include",
        "internal",
    ],
    linkopts = [
        "-lgomp",
    ],
    linkstatic = 1,
    visibility = ["//visibility:public"],
    deps = [
        "@eigen",
        "@glog",
    ],
)