"""Exact conditional logistic likelihood for variable-length binary strata."""
import numpy as np
from scipy.special import logsumexp


def _subset_distribution(eta, counts):
    """Return log subset partition and marginal inclusion probabilities per row.

    The forward recurrence takes or skips each cell. Reverse accumulation of
    its take probabilities gives d log(Z) / d eta for every cell. This requires
    O(groups * cells * max_count) work, independent of the number of subsets.
    Padded cells have log weight -inf and can only be skipped.
    """
    n_groups, n_cells = eta.shape
    width = int(counts.max())
    log_partition = np.full((n_groups, width + 1), -np.inf)
    log_partition[:, 0] = 0.
    take = np.zeros((n_cells, n_groups, width))
    for i in range(n_cells):
        # States beyond i+1 successes are still unreachable.
        limit = min(i + 1, width)
        included = log_partition[:, :limit] + eta[:, i, None]
        combined = np.logaddexp(log_partition[:, 1:limit + 1], included)
        transition = take[i, :, :limit]
        possible = np.isfinite(combined)
        np.subtract(included, combined, out=transition, where=possible)
        np.exp(transition, out=transition, where=possible)
        log_partition[:, 1:limit + 1] = combined

    rows = np.arange(n_groups)
    adjoint = np.zeros_like(log_partition)
    adjoint[rows, counts] = 1.
    marginal = np.zeros_like(eta)
    for i in range(n_cells - 1, -1, -1):
        limit = min(i + 1, width)
        included_mass = adjoint[:, 1:limit + 1] * take[i, :, :limit]
        marginal[:, i] = included_mass.sum(axis=1)
        adjoint[:, 1:limit + 1] -= included_mass
        adjoint[:, :limit] += included_mass
    return log_partition[rows, counts], marginal


def conditional_loglike(theta, X, y, lengths, weights=None, return_scores=False):
    """Return negative log likelihood and gradient, optionally weighted scores.

    X[g, :lengths[g]] and y[g, :lengths[g]] are the observed design and
    binary labels of stratum g. Each stratum conditions on its success count.
    Empty, constant-outcome and zero-weight strata contribute exactly zero.
    Scores are derivatives of the *positive* log likelihood by stratum.
    """
    design = np.asarray(X, dtype=float)
    labels = np.asarray(y)
    lengths = np.asarray(lengths)
    theta = np.asarray(theta, dtype=float)
    n_groups, n_cells, p = design.shape
    if labels.shape != (n_groups, n_cells) or lengths.shape != (n_groups,) or theta.shape != (p,):
        raise ValueError("incompatible design, outcome, length, or coefficient dimensions")
    if not np.issubdtype(lengths.dtype, np.integer) or np.any((lengths < 0) | (lengths > n_cells)):
        raise ValueError("stratum lengths must be integers within the allocated cell dimension")
    weight = np.ones(n_groups) if weights is None else np.asarray(weights, dtype=float)
    if weight.shape != (n_groups,) or not np.all(np.isfinite(weight)) or np.any(weight < 0):
        raise ValueError("provide one finite nonnegative weight per stratum")
    valid = np.arange(n_cells) < lengths[:, None]
    labels = np.where(valid, labels, 0)
    if not np.all((labels == 0) | (labels == 1)):
        raise ValueError("outcomes for observed cells must be binary")
    successes = labels.sum(axis=1, dtype=int)
    counts = np.minimum(successes, lengths - successes)
    scores = np.zeros((n_groups, p))
    loss = 0.

    # When there are fewer failures, count excluded sites using reversed
    # logits. A single selected site then reduces to a categorical likelihood.
    for multiple in (False, True):
        rows = np.flatnonzero((counts > 0) & (weight > 0) & ((counts > 1) == multiple))
        if rows.size == 0:
            continue
        end = int(lengths[rows].max())
        active = valid[rows, :end]
        x = np.where(active[:, :, None], design[rows, :end], 0.)
        x -= x[:, :1].copy()
        x *= active[:, :, None]
        reverse = successes[rows] > counts[rows]
        sign = np.where(reverse, -1., 1.)
        target = np.where(reverse[:, None], ~labels[rows, :end].astype(bool),
                          labels[rows, :end]).astype(bool) & active
        eta = np.einsum("gnp,p->gn", x, theta) * sign[:, None]
        if not np.all(np.isfinite(eta)):
            raise ValueError("non-finite predictor in an observed stratum")
        if multiple:
            partition, probability = _subset_distribution(np.where(active, eta, -np.inf), counts[rows])
        else:
            logits = np.where(active, eta, -np.inf)
            partition = logsumexp(logits, axis=1)
            probability = np.exp(logits - partition[:, None])
        loss += np.dot(weight[rows], partition - np.einsum("gn,gn->g", target, eta))
        scores[rows] = np.einsum("gn,gnp->gp", target - probability, x) * (sign * weight[rows])[:, None]
    gradient = -scores.sum(axis=0)
    if return_scores:
        return float(loss), gradient, scores
    return float(loss), gradient
