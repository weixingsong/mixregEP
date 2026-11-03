# Banknote analysis using the joint MLE method from mixregMLE.R.
# Model: Diagonal ~ 1 + Right, two EP mixture regression components,
# with component-specific mixing proportions, scales and continuous alphas.
# Same data subset as banknote1.R: exclude 138 <= Diagonal <= 139.
# MLE weighted regression, scale regularizer, and Newton updates are preserved.
# The only fitter adaptations are configurable banknote starts, two regression
# coefficients instead of three, centered x for numerical conditioning, and
# banknote-specific labels (higher fitted mean at mean Right = component 1).
# Centering does not rescale residuals or change the scale regularizer.
# Several starts are fitted to the SAME observed data; failed starts are logged.
# No replacement data or simulation is used. Status is not used in fitting.
# Source/run this script beside banknote_data.csv to create banknote_mle_result,
# tables, classifications, a full RDS result, and two PNG figures.

denEp <- function(y, mu, sigma, alpha) {
  cst <- alpha / (2 * sigma * gamma(1 / alpha))
  cst * exp(-(abs(y - mu) / sigma)^alpha)
}

check_loglik <- function(value, stage) {
  if (length(value) != 1L || !is.finite(value)) {
    stop(paste("Non-finite log-likelihood", stage), call. = FALSE)
  }
  invisible(value)
}

update_alpha <- function(tau, residual, sigma, alpha) {
  z <- residual / sigma
  # For alpha > 0, z^alpha*log(z) and z^alpha*log(z)^2 both tend to 0
  # as z tends to 0. Set these terms to their limits at exact zero residuals.
  term1 <- term2 <- numeric(length(z))
  nz <- z > 0
  term1[nz] <- z[nz]^alpha * log(z[nz])
  term2[nz] <- z[nz]^alpha * log(z[nz])^2
  gafun <- 1 / alpha^2 + 2 * digamma(1 / alpha) / alpha^3 +
    trigamma(1 / alpha) / alpha^4
  pd1 <- sum(tau * (1 / alpha + digamma(1 / alpha) / alpha^2)) -
    sum(tau * term1)
  pd2 <- -sum(gafun * tau) - sum(tau * term2)
  # In particular, component 2 uses alpha2 here, correcting the original
  # alp10 in its Newton denominator to alp20.
  0.001 + (alpha - 0.001) * exp(-pd1 / (pd2 * (alpha - 0.001) + pd1))
}

fit_mle_EP <- function(y, X, init_beta, init_prob = c(0.5, 0.5),
                       init_sigma = c(1, 1), init_alpha = c(1, 1),
                       max_iter = 500L, tol = 1e-6) {
  y <- as.numeric(y)
  prob <- init_prob
  beta <- init_beta
  sigma <- init_sigma
  alpha <- init_alpha
  stopifnot(identical(dim(beta), c(2L, ncol(X))),
            length(prob) == 2L, length(sigma) == 2L, length(alpha) == 2L,
            all(is.finite(c(beta, prob, sigma, alpha))),
            all(prob > 0), abs(sum(prob) - 1) < 1e-10,
            all(sigma > 0), all(alpha > 0))
  mu <- X %*% t(beta)
  mix_density <- prob[1] * denEp(y, mu[, 1], sigma[1], alpha[1]) +
    prob[2] * denEp(y, mu[, 2], sigma[2], alpha[2])
  lik0 <- sum(log(mix_density))
  check_loglik(lik0, "at initialization")
  converged <- FALSE
  loglik_trace <- lik0
  
  for (iteration in seq_len(max_iter)) {
    f1 <- denEp(y, mu[, 1], sigma[1], alpha[1])
    f2 <- denEp(y, mu[, 2], sigma[2], alpha[2])
    fsum <- prob[1] * f1 + prob[2] * f2
    if (any(!is.finite(fsum)) || any(fsum <= 0)) {
      stop("Invalid mixture density in E-step", call. = FALSE)
    }
    tau <- cbind(prob[1] * f1 / fsum, prob[2] * f2 / fsum)
    prob_new <- c(mean(tau[, 1]), 1 - mean(tau[, 1]))
    beta_new <- beta
    sigma_new <- alpha_new <- numeric(2)
    
    for (component in 1:2) {
      residual <- abs(y - mu[, component])
      mass <- sum(tau[, component])
      if (!is.finite(mass) || mass <= 0) {
        stop("Empty mixture component", call. = FALSE)
      }
      weights <- sqrt(tau[, component] *
                        (residual + 1e-6)^(alpha[component] - 2))
      if (any(!is.finite(weights))) {
        stop("Non-finite regression weights", call. = FALSE)
      }
      # Row scaling is equivalent to the previous diag(weights) multiplication.
      Xs <- weights * X
      ys <- weights * y
      beta_new[component, ] <- as.numeric(solve(crossprod(Xs), crossprod(Xs, ys)))
      sigma_new[component] <- sqrt(
        (alpha[component] * sum(tau[, component] * residual^alpha[component]) /
           mass)^(2 / alpha[component]) + 0.0001)
      alpha_new[component] <- update_alpha(tau[, component], residual,
                                           sigma_new[component], alpha[component])
    }
    if (any(!is.finite(c(prob_new, beta_new, sigma_new, alpha_new))) ||
        any(prob_new <= 0) || any(sigma_new <= 0) || any(alpha_new <= 0)) {
      stop("Non-finite or invalid updated parameters", call. = FALSE)
    }
    mu_new <- X %*% t(beta_new)
    lik1 <- sum(log(
      prob_new[1] * denEp(y, mu_new[, 1], sigma_new[1], alpha_new[1]) +
        prob_new[2] * denEp(y, mu_new[, 2], sigma_new[2], alpha_new[2])))
    check_loglik(lik1, paste("at iteration", iteration))
    loglik_trace <- c(loglik_trace, lik1)
    converged <- abs(lik1 - lik0) < tol
    
    # Store the latest parameters BEFORE testing the stopping condition.
    # The previous script otherwise returned parameters from the preceding step.
    prob <- prob_new
    beta <- beta_new
    sigma <- sigma_new
    alpha <- alpha_new
    mu <- mu_new
    lik0 <- lik1
    if (converged) break
  }
  # Simulation-truth labels from mixregMLE.R do not apply to banknotes.
  # Label all parameters together after restoring the original x scale.
  fit <- list(prob = prob, beta = beta, sigma = sigma, alpha = alpha)
  c(fit, list(loglik = lik1, iterations = iteration, converged = converged,
              loglik_trace = loglik_trace))
}

# Locate the bundled CSV beside this script (works with source() or Rscript).
banknote_script_dir <- local({
  files <- unlist(lapply(sys.frames(), function(e) e$ofile), use.names = FALSE)
  args <- commandArgs(trailingOnly = FALSE)
  arg_file <- sub("^--file=", "", args[grepl("^--file=", args)])
  file <- if (length(files)) tail(files, 1L) else if (length(arg_file)) arg_file[1L] else NULL
  if (is.null(file)) getwd() else dirname(normalizePath(file, mustWork = TRUE))
})

load_banknote <- function(data_file = file.path(banknote_script_dir, "banknote_data.csv")) {
  if (file.exists(data_file)) return(read.csv(data_file, stringsAsFactors = FALSE))
  if (requireNamespace("mclust", quietly = TRUE)) {
    e <- new.env()
    data("banknote", package = "mclust", envir = e)
    return(e$banknote)
  }
  stop("Place banknote_data.csv beside this script, or install mclust.", call. = FALSE)
}

row_log_sum_exp <- function(A) {
  m <- apply(A, 1L, max)
  if (any(!is.finite(m))) stop("Non-finite posterior log-density", call. = FALSE)
  m + log(rowSums(exp(A - m)))
}

logdenEp <- function(y, mu, sigma, alpha) {
  log(alpha) - log(2) - log(sigma) - lgamma(1 / alpha) -
    (abs(y - mu) / sigma)^alpha
}

posterior_EP <- function(fit, x, y) {
  X <- cbind(1, x)
  mu <- X %*% t(fit$beta)
  log_joint <- cbind(
    log(fit$prob[1]) + logdenEp(y, mu[, 1], fit$sigma[1], fit$alpha[1]),
    log(fit$prob[2]) + logdenEp(y, mu[, 2], fit$sigma[2], fit$alpha[2]))
  log_mix <- row_log_sum_exp(log_joint)
  list(posterior = exp(log_joint - log_mix), loglik = sum(log_mix))
}

analyze_banknote_mle <- function(banknote = load_banknote(), max_iter = 500L,
                                 tol = 1e-6) {
  stopifnot(all(c("Right", "Diagonal", "Status") %in% names(banknote)),
            all(is.finite(banknote$Right)), all(is.finite(banknote$Diagonal)),
            length(max_iter) == 1L, is.finite(max_iter), max_iter >= 1L,
            max_iter == floor(max_iter), length(tol) == 1L, is.finite(tol), tol > 0)
  excluded <- banknote$Diagonal >= 138 & banknote$Diagonal <= 139
  used <- !excluded
  x <- banknote$Right[used]
  y <- banknote$Diagonal[used]
  x_center <- mean(x)
  X <- cbind(1, Right_centered = x - x_center)
  
  original_beta <- rbind(c(143, -0.02), c(75, -2.51))
  empirical_beta <- rbind(c(as.numeric(quantile(y, .75)), 0),
                          c(as.numeric(quantile(y, .25)), 0))
  beta_starts <- list(banknote1 = original_beta, empirical = empirical_beta)
  alpha_starts <- rbind(c(.5, .5), c(1, 1), c(1.5, 1.5), c(2, 2),
                        c(.7, 1.5), c(1.5, .7))
  sigma_starts <- list(separate = c(.34, .16), common = c(.25, .25))
  fits <- list()
  diagnostics <- list()
  start_id <- 0L
  t0 <- proc.time()[["elapsed"]]
  for (b in names(beta_starts)) for (s in names(sigma_starts)) {
    for (a in seq_len(nrow(alpha_starts))) {
      start_id <- start_id + 1L
      beta0 <- beta_starts[[b]]
      beta0[, 1] <- beta0[, 1] + beta0[, 2] * x_center
      warning_messages <- character()
      reason <- ""
      fit <- tryCatch(withCallingHandlers(
        fit_mle_EP(y, X, init_beta = beta0,
                   init_sigma = sigma_starts[[s]], init_alpha = alpha_starts[a, ],
                   max_iter = max_iter, tol = tol),
        warning = function(w) {
          warning_messages <<- c(warning_messages, conditionMessage(w))
          invokeRestart("muffleWarning")
        }), error = function(e) {
          reason <<- conditionMessage(e)
          NULL
        })
      # [id] <- list(fit) keeps the slot even when fit is NULL.
      fits[start_id] <- list(fit)
      diagnostics[[start_id]] <- data.frame(
        start = start_id, beta_start = b, sigma_start = s,
        alpha1_start = alpha_starts[a, 1], alpha2_start = alpha_starts[a, 2],
        success = !is.null(fit), converged = if (is.null(fit)) FALSE else fit$converged,
        iterations = if (is.null(fit)) NA_integer_ else fit$iterations,
        loglik = if (is.null(fit)) NA_real_ else fit$loglik,
        decreases = if (is.null(fit)) NA_integer_ else sum(diff(fit$loglik_trace) < -tol),
        reason = reason, warnings = paste(unique(warning_messages), collapse = "; "),
        stringsAsFactors = FALSE)
    }
  }
  diagnostics <- do.call(rbind, diagnostics)
  ok <- which(diagnostics$success & is.finite(diagnostics$loglik))
  if (!length(ok)) stop("All starting values failed; inspect fitter and data.", call. = FALSE)
  chosen <- ok[which.max(diagnostics$loglik[ok])]
  fit <- fits[[chosen]]
  diagnostics$selected <- diagnostics$start == chosen
  
  # Restore intercepts for raw Right measurements.
  fit$beta[, 1] <- fit$beta[, 1] - fit$beta[, 2] * x_center
  # Identify components by their fitted response at the observed mean predictor,
  # with higher fitted Diagonal first. The known Status is never used here.
  ord <- order(drop(fit$beta %*% c(1, x_center)), decreasing = TRUE)
  fit$beta <- fit$beta[ord, , drop = FALSE]
  fit$prob <- fit$prob[ord]
  fit$sigma <- fit$sigma[ord]
  fit$alpha <- fit$alpha[ord]
  rownames(fit$beta) <- c("Component 1", "Component 2")
  colnames(fit$beta) <- c("Intercept", "Right")
  post <- posterior_EP(fit, x, y)
  stopifnot(abs(post$loglik - fit$loglik) < 1e-7,
            all(abs(rowSums(post$posterior) - 1) < 1e-10))
  fit$loglik <- post$loglik
  fit$posterior <- post$posterior
  fit$selected_start <- chosen
  fit$x_center <- x_center
  fit$runtime_sec <- proc.time()[["elapsed"]] - t0
  
  parameters <- data.frame(component = c("Component 1", "Component 2"),
                           pi = fit$prob, intercept = fit$beta[, 1],
                           slope_Right = fit$beta[, 2], sigma = fit$sigma,
                           alpha = fit$alpha, row.names = NULL)
  # sigma is the EP scale, not the error standard deviation.
  parameters$error_sd <- fit$sigma * exp(
    .5 * (lgamma(3 / fit$alpha) - lgamma(1 / fit$alpha)))
  # 2*(2 beta + 1 sigma + 1 alpha) + 1 free mixing probability = 9.
  npar <- 9L
  summary <- data.frame(
    n_original = nrow(banknote), n_used = sum(used), n_excluded = sum(excluded),
    parameters = npar, loglik = fit$loglik,
    AIC = -2 * fit$loglik + 2 * npar,
    BIC = -2 * fit$loglik + log(sum(used)) * npar,
    iterations = fit$iterations, converged = fit$converged,
    selected_start = chosen, starts_tried = nrow(diagnostics),
    finite_starts = sum(diagnostics$success),
    converged_starts = sum(diagnostics$converged),
    runtime_sec = fit$runtime_sec)
  classified <- data.frame(row_id = which(used), Status = banknote$Status[used],
                           Right = x, Diagonal = y,
                           posterior_component1 = post$posterior[, 1],
                           posterior_component2 = post$posterior[, 2],
                           component = max.col(post$posterior, ties.method = "first"))
  # These names summarize the observed higher/lower Diagonal groups only;
  # they are a post-fit comparison with recorded Status, not supervised fitting.
  classified$predicted_Status <- ifelse(classified$component == 1L, "genuine", "counterfeit")
  confusion <- table(Recorded = classified$Status, Fitted = classified$predicted_Status)
  list(fit = fit, parameters = parameters, summary = summary,
       start_diagnostics = diagnostics, classifications = classified,
       confusion = confusion, banknote = banknote,
       included_rows = which(used), excluded_rows = which(excluded),
       settings = list(max_iter = max_iter, tol = tol,
                       exclusion = "138 <= Diagonal <= 139", x = "Right", y = "Diagonal"))
}

plot_banknote_mle <- function(result) {
  old_par <- par(mar = c(8, 5, 4.2, 1.5))
  on.exit(par(old_par))
  d <- result$banknote
  included <- seq_len(nrow(d)) %in% result$included_rows
  status_col <- ifelse(d$Status == "genuine", "#C6453D", "#2764A5")
  status_col[d$Diagonal < 138] <- "black"
  plot(d$Right[included], d$Diagonal[included],
       col = status_col[included], pch = ifelse(d$Diagonal[included] < 138, 8, 19),
       xlim = range(d$Right), ylim = range(d$Diagonal) + c(-.1, .15),
       xlab = "Right", ylab = "Diagonal",
       main = "Banknote data: two-component EP mixture MLE")
  points(d$Right[!included], d$Diagonal[!included], col = "gray60", pch = 1)
  for (k in 1:2) abline(result$fit$beta[k, ], col = c("#C6453D", "#2764A5")[k], lwd = 2)
  legend("topright", legend = c("Recorded genuine (included)", "Recorded counterfeit (included)",
                                "Included low-diagonal notes (<138)",
                                "Excluded: 138 <= Diagonal <= 139"),
         col = c("#C6453D", "#2764A5", "black", "gray60"), pch = c(19, 19, 8, 1),
         bty = "n", cex = .78)
  p <- result$parameters
  for (k in 1:2) mtext(
    sprintf("Component %d: y = %.3f %+.3f x; alpha = %.3f; pi = %.3f",
            k, p$intercept[k], p$slope_Right[k], p$alpha[k], p$pi[k]),
    side = 1, line = 4.4 + 1.1 * (k - 1L), adj = 0,
    col = c("#C6453D", "#2764A5")[k], cex = .78)
  invisible(result)
}

format_banknote_3dp <- function(tab) {
  ints <- c("component", "start", "iterations", "selected_start", "starts_tried", "finite_starts",
            "converged_starts", "n_original", "n_used", "n_excluded", "parameters", "decreases")
  for (name in names(tab)) if (is.numeric(tab[[name]]) && !name %in% ints) {
    x <- tab[[name]]
    x[is.finite(x) & round(x, 3) == 0] <- 0
    tab[[name]] <- ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = 3))
  }
  tab
}

save_banknote_mle <- function(result, output_dir = banknote_script_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(output_dir, "banknote_MLE_results.rds"))
  write.csv(format_banknote_3dp(result$parameters),
            file.path(output_dir, "banknote_MLE_parameters.csv"), row.names = FALSE)
  write.csv(format_banknote_3dp(result$summary),
            file.path(output_dir, "banknote_MLE_summary.csv"), row.names = FALSE)
  write.csv(format_banknote_3dp(result$start_diagnostics),
            file.path(output_dir, "banknote_MLE_start_diagnostics.csv"), row.names = FALSE)
  # Posterior probabilities retain full precision for further calculations.
  write.csv(result$classifications, file.path(output_dir, "banknote_MLE_classifications.csv"),
            row.names = FALSE)
  png(file.path(output_dir, "banknote_MLE_fit.png"), width = 1800, height = 1400, res = 180)
  tryCatch(plot_banknote_mle(result), finally = dev.off())
  png(file.path(output_dir, "banknote_MLE_convergence.png"), width = 1800, height = 1000, res = 180)
  tryCatch({
    plot(seq_along(result$fit$loglik_trace) - 1L, result$fit$loglik_trace, type = "l",
         xlab = "Iteration", ylab = "Observed log-likelihood",
         main = "Selected MLE starting value: likelihood trace", lwd = 2, col = "#2764A5")
  }, finally = dev.off())
  invisible(result)
}

# Analyze the observed banknote subset once, using deterministic multiple starts.
banknote_mle_result <- analyze_banknote_mle()
print(format_banknote_3dp(banknote_mle_result$parameters), row.names = FALSE)
print(format_banknote_3dp(banknote_mle_result$summary), row.names = FALSE)
print(banknote_mle_result$confusion)
save_banknote_mle(banknote_mle_result)
if (interactive()) plot_banknote_mle(banknote_mle_result)
