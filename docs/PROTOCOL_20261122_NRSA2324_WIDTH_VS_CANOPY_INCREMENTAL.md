# After observing additive B+NP+H: distinguish stream width from actual canopy, 2026-11-22

**Explicitly posthoc decomposition**, motivated by the previous positive out-of-basin MMI prediction from adding *width plus canopy together*. Not a confirmatory study.

Exactly same original EPANRSA source-SHA, native UID/visit/date metadata, complete-case method-specific strata: BERW MMI/EPT 1,165 visits/693 HUC8 and BETB MMI/EPT 725 visits/476 HUC8. No change in outcome, flags, region assignments, original date, or sample selection.

Four fixed nested models, all include measured in-stream original TN+TP and same demographic/geographic/season/conductivity controls:
  1. N: base + measured TN+TP
  2. N+W: N+original log(width)
  3. N+C: N+original canopy percent (XCDENMID BERW; XCDENBK BETB)
  4. N+W+C: N+width+canopy (must match original exact baseline MSE).
Four comparisons with full HUC8-disjoint heldout 5fold Ridge10, train-fold-only feature preprocessing, HUC8-equal MSE and 2000 whole-HUC8 fixedOOF bootstrap descriptive95:
  width added after N,
  canopy added after N,
  width added after N+canopy,
  canopy added after N+width.
Keep all four original strata/outcomes and nulls; source STOP if exact reference changed. A positive canopy term only indicates independent prediction, **not** an actual effect of shading or riparian restoration, since light, heat, river flow, urbanization and substrate remain unmeasured.
