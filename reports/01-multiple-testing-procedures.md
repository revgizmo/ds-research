# Multiple-testing correction for fixed-horizon A/B-test reports

## Scope, assumptions, and bottom line

This report assumes that each reported hypothesis has a valid marginal two-sided frequentist p-value at the fixed analysis horizon; the experiment-level family to be protected has not yet been fixed; variants are randomized independently of potential outcomes; and there is no sequential monitoring, alpha spending, online FDR, bandit logic, Bayesian decision rule, or covariate-adjustment issue. I treat a "claim" as whatever declaration the reporting system permits the reader or downstream decision process to act on as statistically established. I also distinguish metrics from metric-by-segment cells because the latter can multiply the effective family dramatically.

A central distinction runs through all four questions. Bonferroni and Holm provide strong family-wise error rate control without assumptions about dependence among valid p-values. The usual generic Sidak correction needs independence, although Sidak's original 1967 Gaussian result is stronger for a particular symmetric multivariate-normal problem. Hochberg and Hommel depend on Simes-type inequalities and therefore require substantially more structure than Holm. BH was proved originally under independence and later under positive regression dependence on a subset, or PRDS. BY obtains FDR control under arbitrary dependence by replacing the BH level \(q\) with \(q/H_m\), where

\[
H_m=\sum_{i=1}^m \frac1i.
\]

Storey's 2002 adaptive/q-value framework was developed under independent p-values; stronger finite-sample statements for adaptive FDR procedures came in later work. For the A/B setting described, "the metrics are all measured on the same units" is not a positive-dependence assumption. Two Bernoulli metrics on the same unit can be positively correlated, independent, or negatively correlated. Mutually exclusive outcomes are negatively correlated. Nested outcomes are positively correlated. Several treatment-control comparisons of the same endpoint share a control and ordinarily have positive correlation, which is the classical Dunnett structure. But combining that with different binary endpoints, two-sided tests, and overlapping segment cells does not in general establish the PRDS or Simes conditions needed by BH, Hochberg, or Hommel. BY explicitly discusses many-to-one comparisons as an important positive-dependence setting, but a very recent counterexample shows that ordinary BH can exceed its nominal FDR for correlated two-sided Gaussian tests, including a finite example with 85 tests. There is also no generally accepted rule that "12 metrics means FWER" and "50 metrics means FDR." The error criterion is primarily a statement about what makes a false claim costly. Regulators frame confirmatory efficacy claims in terms of study-wise Type I error, effectively strong FWER. BH's original motivation was different: if a small fraction of erroneous discoveries is tolerable among a larger set of discoveries, FDR can be a more relevant target. A recent NCI-authored review proposes FWER for a few endpoints, FDX for roughly more than a dozen but fewer than several hundred, and FDR for thousands, while explicitly presenting this as guidance rather than an immutable rule. Experimentation practitioners disagree as well: Microsoft's ExP has moved toward BH-style FDR handling for scorecards with hundreds or thousands of metrics, while Spotify Confidence argues that FWER often better matches product shipping decisions. The evidence therefore supports trade-offs, not a single universally dominant correction. Procedures that need less dependence structure tend to make stronger error guarantees at the price of power; procedures that exploit weaker error criteria or stronger dependence assumptions can gain power, especially when many alternatives are nonnull.

## Error-rate guarantees and dependence conditions

**Plain-language answer.** Bonferroni and Holm are the cleanest with respect to dependence: valid marginal p-values are enough for strong FWER control. Standard Sidak is exact under independence but is not generally safe for arbitrary dependent p-values. Hochberg and Hommel gain power by relying on Simes-type structure, so they do not inherit Holm's arbitrary-dependence robustness. Original BH controls FDR under independence; BY showed that ordinary BH remains valid under PRDS and supplied the harmonic correction for arbitrary dependence. Storey's original adaptive framework assumes independence. The most consequential qualification for this use case is that positive pairwise correlations are not themselves a proof of the required positive-dependence condition, particularly for two-sided p-values.

**Confidence: high** for Bonferroni, Holm, BH, BY, and Storey's original independence assumptions; **medium-high** for the exact historical wording of Hochberg and Hommel because I could access their abstracts but not their full original articles. Confidence there would increase if the 1988 Biometrika full texts were openly accessible. The mathematical modern characterization through the Simes inequality is high confidence.

| Procedure | Error rate controlled | Dependence condition supporting the guarantee | Important limitation or counterexample |
|---|---|---|---|
| **Bonferroni** | Strong FWER: \(P(V\ge1)\le\alpha\) under every configuration of true and false nulls. | No dependence restriction beyond valid, super-uniform null p-values. It follows from the union bound. Holm explicitly starts from this logic and extends it sequentially. | Conservative when rejection events overlap strongly. There is no dependence counterexample to its FWER guarantee if the individual p-values are valid. |
| **Sidak** | Strong FWER at \(\alpha\) for the usual single-step cutoff \(1-(1-\alpha)^{1/m}\) when null p-values are independent. | Independence is the generic p-value condition. Sidak's **original 1967 paper** proves a more specialized result: for multivariate-normal coordinates, a symmetric rectangular simultaneous region constructed as under independence is conservative for arbitrary correlation. | The generic p-value version can fail with negative dependence. For example, with \(p_1=U,p_2=1-U\), \(U\sim U(0,1)\), the two Sidak rejection events are disjoint when the cutoff is small, so their union has probability \(2[1-\sqrt{1-\alpha}]>\alpha\). |
| **Holm** | Strong FWER. | No dependence restriction beyond valid p-values. Holm states protection "for any combination of true hypotheses." | No dependence counterexample under valid p-values. It uniformly dominates single-step Bonferroni in the sense that it can never reject fewer hypotheses for a fixed vector of p-values. |
| **Hochberg** | Strong FWER. | Independence is the classical condition. More generally, validity follows when the relevant intersection-null p-values satisfy the Simes inequality. Later positive-dependence results provide sufficient conditions, including classes based on MTP2. | Simes can fail under negative dependence, so Hochberg is not an arbitrary-dependence procedure. Later work gives explicit negative-dependence failures and modified procedures for that setting. |
| **Hommel** | Strong FWER. | Hommel is a closed-testing construction using Simes local tests. Its guarantee therefore inherits the conditions under which those Simes tests are valid. | It does not have Holm's arbitrary-dependence guarantee. Failure of the underlying Simes inequality can invalidate the claimed level. |
| **Benjamini-Hochberg (BH)** | FDR, \(E[V/\max(R,1)]\). | BH 1995 proves control under independent test statistics. BY 2001 later proves ordinary BH under PRDS on the subset of true-null statistics. | Arbitrary dependence can violate nominal FDR. Most relevant here, Dobriban 2026 constructs correlated **two-sided Gaussian** tests for which BH exceeds nominal FDR, including an \(m=85\) example. |
| **Benjamini-Yekutieli (BY)** | FDR. | Arbitrary dependence, assuming valid null p-values, after harmonic correction. BY Theorem 1.3 uses the BH critical values divided by \(H_m\). | The price is conservatism. There is no dependence counterexample to the theorem under its marginal-p-value assumptions. |
| **Storey 2002 adaptive/q-value procedure** | The paper centers on pFDR and estimation of FDR/pFDR, with q-values and adaptive estimation of the null fraction \(\pi_0\). | Storey explicitly restricts the main 2002 development to independent p-values. | "Storey's procedure" is sometimes used loosely for later adaptive BH variants. The 2002 paper should not be read as a general finite-sample arbitrary-dependence theorem. Storey, Taylor, and Siegmund 2004 supplied later finite-sample/asymptotic FDR-control results for point-estimate-based procedures. |

The exact BY factor is

\[
c(m)=H_m=1+\frac12+\frac13+\cdots+\frac1m.
\]

Thus a BH critical value \(iq/m\) becomes \(iq/[mH_m]\). BY's theorem gives FDR no greater than \((m_0/m)q\), hence no greater than \(q\). For scale, \(H_{12}\approx3.1032\), \(H_{50}\approx4.4992\), and \(H_{100}\approx5.1874\). These are mathematical evaluations of the exact factor, not estimates.

There are two historically important nuances. First, Sidak is often summarized too casually as "requires independence." That is correct for arbitrary p-values, but his 1967 paper proves the multivariate-normal rectangular probability inequality that makes the usual simultaneous Gaussian intervals conservative even under arbitrary Gaussian correlation. Second, statements such as "Hochberg is valid under positive correlation" are too loose. Sarkar's work concerns specific positive-dependence properties such as MTP2 and the validity of the Simes inequality, not every joint distribution having nonnegative pairwise Pearson correlations. A related hierarchy is useful. Independence is sufficient for Simes and original BH. Certain positive-dependence structures are also sufficient. BY's PRDS condition asks, roughly, that conditioning a true-null statistic toward more significant values cannot make increasing events in the remaining statistics less likely. It is a joint distributional condition, not simply a covariance sign. MTP2 is stronger still and supplies Simes-type inequalities in important multivariate-normal families with nonnegative correlations. Sarkar also discusses Dunnett/Sobel-type statistics as examples of positively dependent test statistics covered by this theory. ### Key evidence for the guarantees

| Key claim | Source and location | Passage actually seen |
|---|---|---|
| Holm gives Type I protection regardless of which hypotheses are true. | Holm 1979, abstract and procedure discussion. | "protection against error of the first kind for any combination of true hypotheses" |
| Simes is exact under independence. | Simes 1986, theorem/procedure discussion. | "If all hypotheses are true and the tests are independent, the probability of rejection is exactly alpha." |
| Sidak's original Gaussian result allows arbitrary correlation. | Sidak 1967, theorem on multivariate normal rectangles. | The paper states the result for a multivariate normal distribution with an arbitrary correlation matrix; full passage was seen, but I do not reproduce a longer quotation here. |
| BH's original theorem assumes independence. | Benjamini and Hochberg 1995, Theorem 1. | "For independent test statistics ... the above procedure controls the FDR at \(q^*\)." |
| Ordinary BH is valid under PRDS. | Benjamini and Yekutieli 2001, Theorem 1.2. | "If the joint distribution of the test statistics is PRDS on the subset of test statistics corresponding to true null hypotheses" |
| BY works under arbitrary dependence with the harmonic divisor. | Benjamini and Yekutieli 2001, Theorem 1.3. | The theorem gives \(q/\sum_{i=1}^m1/i\); exact passage seen. |
| Storey 2002 works under independence. | Storey 2002, assumptions preceding Theorem 1. | "This paper will be limited to the situation where the p-values are independent." |
| Hochberg and Hommel rest on Simes validity. | Sarkar 1998, introduction, discussion of Hochberg and Hommel. | "both procedures control the Type I error rate in situations where the Simes method does" |
| Simes need not be conservative under negative dependence. | Block, Savits, and Wang, abstract. | Abstract seen; full text not opened. The reported result is that the Simes inequality is reversed for a broad negatively dependent class. |
| Two-sided correlated Gaussian tests can break ordinary BH. | Dobriban 2026, abstract and main theorem. | "we settle the question negatively, by providing examples of Gaussian factor models for which the BH method fails" |

For Hochberg 1988 and Hommel 1988, the full original Biometrika articles were not accessible during this research. Their abstracts establish the procedures and claimed multiplicity control, but the exact original dependence wording was **not seen**. I therefore rely on Sarkar's later mathematical description for the Simes-based dependence condition rather than inventing wording from the originals. ## Dependence of binary metrics, shared controls, and segment cells

**Plain-language answer.** Binary metrics on the same users or visits do not have a fixed dependence direction. The sign is determined by their joint event structure. Multiple variant-control comparisons of the same metric do have the classical positive shared-control correlation. These are different sources of dependence. Once many different outcomes and overlapping segment cells are combined in a single family, neither the sign nor a strong condition such as PRDS follows merely from randomization or from sharing units.

**Confidence: high** for the covariance statements and shared-control result; **medium-high** for statements about whether an entire real scorecard satisfies PRDS, because that is a property of its full joint distribution. Confidence for a particular implementation would change only after specifying its exact test statistics and joint event/segment structure, or establishing the required inequality analytically.

For two Bernoulli outcomes \(A\) and \(B\) on the same observational unit,

\[
\operatorname{Cov}(1_A,1_B)
=P(A\cap B)-P(A)P(B).
\]

That identity immediately gives three qualitatively different cases:

- If \(A\) and \(B\) are mutually exclusive, \(P(A\cap B)=0\), so the covariance is \(-P(A)P(B)<0\) when both have positive probability.
- If \(B\subseteq A\), then the covariance is \(P(B)[1-P(A)]\ge0\).
- If the events are independent, covariance is zero.

For counts over independent experimental units, the count covariance is \(n\) times the corresponding unit-level covariance; for sample proportions it is the unit covariance divided by \(n\). Thus, common measurement on the same units creates dependence but does not specify its sign.

This matters directly for mutually exclusive funnel or categorical outcomes. "Clicked A" and "clicked B but not A," competing conversion categories, different failure causes, or mutually exclusive message outcomes can be negatively associated even though every outcome comes from exactly the same randomized visit. In contrast, funnel events such as "purchased" nested inside "reached checkout" generally have positive unit-level association by construction. Neither fact alone tells us the dependence among their final two-sided test p-values.

Shared controls behave differently. Suppose \(\hat\theta_i=\hat\mu_i-\hat\mu_0\) is the estimated effect for treatment arm \(i\) versus a common control for one outcome, and treatment arms are otherwise independent. Then for \(i\ne j\),

\[
\operatorname{Cov}(\hat\theta_i,\hat\theta_j)
=\operatorname{Var}(\hat\mu_0)>0.
\]

Under equal allocation and equal outcome variance, the standardized correlations are commonly about \(1/2\). More generally, they depend on allocation and variances. This is the structure Dunnett's 1955 method was designed to exploit when comparing several treatments with a common control. FDA's multiple-endpoint guidance likewise recognizes Dunnett-type procedures as correlation-aware alternatives in suitable many-to-one comparisons. BY specifically notes that PRDS includes important multiple-comparison problems such as comparing several treatments against a single control, as well as positively correlated multivariate-normal and multivariate-\(t\) test statistics under the stated structures. Sarkar's work similarly establishes Simes inequalities for specific positively dependent multivariate distributions, including Dunnett/Sobel-type test statistics. But there are three reasons this does not settle the user's complete reporting system.

First, **different metrics add their own covariance signs**. For treatment-control estimates of metrics \(A\) and \(B\), covariance contains within-arm cross-metric covariance terms. Those can be positive or negative according to the event structures above. A positive common-control contribution for the same endpoint does not force all cross-endpoint statistics to be positively dependent.

Second, **two-sided p-values discard the sign of the statistic**. PRDS arguments that work naturally for one-sided monotone rejection regions do not automatically transfer to \(p_i=2\bar\Phi(|Z_i|)\). Dobriban's October 2026 result is especially relevant: it supplies positively correlated Gaussian factor models with standard two-sided Gaussian p-values in which ordinary BH exceeds nominal FDR. The main constructions are large, but the supplementary finite example uses only 85 tests. The reported excess in that small example is tiny; the paper's larger violations arise for much larger families. That is a counterexample to the idea that "Gaussian correlation plus two-sided tests" by itself guarantees BH.

Third, **overlapping segments produce another dependence layer**. A metric for "US users," "mobile users," and "US mobile users" reuses observations. A parent segment and child segment are nested; partially overlapping combinations share some units; mutually exclusive segment values use disjoint units. Consequently, test statistics can have nested, overlapping, or near-independent relationships. Neither randomization nor the binary nature of the outcomes implies that the entire metric-by-cell p-value vector is PRDS.

The robustness consequences are therefore:

| Procedure | Does same-unit binary dependence threaten its stated guarantee? | Does shared-control dependence threaten it? |
|---|---|---|
| Bonferroni | No, assuming valid marginal p-values. | No. |
| Holm | No. | No. |
| Generic Sidak | Potentially yes, unless independence or a specialized Gaussian rectangle result applies. | Generic p-value independence fails with a shared control. |
| Hochberg | Potentially yes if Simes validity is not established. | Classical positive Dunnett structures can satisfy sufficient positive-dependence conditions, but mixed metrics/two-sided cells need separate justification. |
| Hommel | Same Simes issue. | Same. |
| BH | Original independence is absent. PRDS may hold in special structures but is not implied globally. | Many-to-one structures are among BY's favorable PRDS examples, but two-sided mixed-metric reports are not automatically covered. |
| BY | No dependence restriction on the joint p-values beyond valid marginals. | No. |
| Storey 2002 | Yes, because the original analysis assumes independent p-values. | A shared control violates literal independence. |

### Evidence on dependence

| Key claim | Source and location | Passage actually seen |
|---|---|---|
| Dunnett treats several treatment comparisons with a common control. | Dunnett 1955, title and bibliographic record. | Full article not opened; title/record seen. |
| The common-control correlations are positive in the standard proportions setting. | Modern derivation/tutorial for comparisons of proportions with a shared control. | Formula seen; it gives positive correlations for positive allocation ratios. |
| PRDS includes several-treatment-versus-one-control examples. | Benjamini and Yekutieli 2001, dependence examples following PRDS discussion. | Passage seen; it explicitly includes "many treatments with a single control." |
| Positive dependence sufficient for Simes is stronger than a covariance sign. | Sarkar 1998, sections on positive dependence and MTP2. | The article develops MTP2/Simes conditions and multivariate-normal cases, not merely pairwise correlation. |
| Two-sided Gaussian correlation alone does not guarantee BH. | Dobriban 2026, abstract, Eq. 1, main theorem, discussion. | "The two-sided rejection event ... is not monotone in that common shift" |
| Correlated binary endpoints are a recognized multiple-testing setting in their own right. | Recent methodological literature on multiple correlated binary endpoints. | Full web text/abstract seen; it explicitly treats correlated binary endpoints and correlation-aware multiplicity methods. |

The elementary covariance identities above are algebraic consequences of Bernoulli indicators rather than empirical claims. No assumption about the marginal rates can turn mutually exclusive positive-probability events into positively correlated indicators.

## Published power comparisons

**Plain-language answer.** Power differences are highly scenario-specific. Holm weakly dominates Bonferroni in the set of rejections, but their probability of detecting at least one true effect can be identical because both start at the same \(\alpha/m\) threshold. Hochberg is usually a modest improvement over Holm in small confirmatory families when its dependence assumptions hold. FDR methods can gain much more when many hypotheses are nonnull, and their advantage grows with family size in the original BH simulations. When only one or very few hypotheses are truly nonnull, much of that gain disappears. I found no single peer-reviewed simulation that exactly matches all eight requested procedures, binary A/B metrics, two-sided tests, overlapping cells, and exactly \(m=12\) and \(m=50\); the closest evidence should not be averaged across settings.

**Confidence: high** for the qualitative ordering and the reported results below; **medium** for extrapolation from the studies to the requested experimentation system. Confidence would increase with a simulation using the actual joint binary endpoint and segment-generating process.

### BH's original simulations: family sizes from four to 64

Benjamini and Hochberg simulated independent normal test statistics with \(m=4,8,16,32,64\), comparing Bonferroni, Hochberg, and BH at \(q=.05\), under several fractions of true nulls: \(3m/4\), \(m/2\), \(m/4\), and zero. Thus \(m=8\) and \(16\) bracket the requested "about 12," and \(m=64\) is directly informative for "50 or more," without pretending that these are identical settings. The authors report that BH's power was uniformly higher in their simulations, with the gain increasing both as the number of nonnull hypotheses increased and as the number of tests increased. They explicitly remark that the gain of BH over the FWER methods was much larger than Hochberg's gain over Bonferroni. One example they give is \(m=32\) with no true nulls, where reported power was about 0.42 for Bonferroni versus 0.65 for BH; in a much smaller example with \(m=4\) and half the hypotheses false, the corresponding values were roughly 0.62 and 0.70. **Study setting:** independent normal statistics, fixed family sizes 4 to 64, simulated alternative configurations, FDR target .05. It is not a study of correlated binary A/B metrics.

**Interpretation for sparse versus dense alternatives:** the original simulations support a much larger FDR power advantage in dense-effect configurations than sparse ones. The theoretical intuition is also in the BH paper: when many hypotheses are nonnull, FDR control can permit substantially more rejections than procedures that protect against even one false rejection. ### FWER comparisons in small RCT families

Vickerstaff and colleagues studied multiple primary outcomes in randomized controlled trials, varying the number of outcomes, correlations, effect sizes, and missingness. Their comparison included Bonferroni, Holm, Hochberg, Hommel, and other procedures. An important result is that **Bonferroni and Holm had the same disjunctive power**, meaning the probability of detecting at least one true treatment effect. This is not contradictory to Holm's uniform improvement in rejection sets. The first Holm cutoff is the same \(\alpha/m\) cutoff as Bonferroni, so Holm's extra power appears primarily in making additional rejections after at least one sufficiently small p-value has occurred. For one reported two-outcome scenario with zero correlation, the table gives approximately 72.4% marginal power under Bonferroni, 78.5% under Holm, and 79.2% under Hochberg and Hommel. The precise numbers change with the simulated correlation and design. The paper's overall conclusion is that Hochberg and Hommel offered power gains over Bonferroni, but those gains were small in the studied RCT settings. **Study setting:** only a few multiple primary outcomes, mostly continuous-outcome clinical-trial simulations, not 12 or 50 tests. It is valuable for the requested **Holm vs Bonferroni** and **Hochberg vs Holm** mechanism but should not be numerically extrapolated to 50 correlated binary metrics.

### Eight-metric product-experiment simulation and 1,300 experiments

Spotify Confidence published a 2026 practitioner analysis of about 1,300 recent experiments containing 1,840 treatment-control comparisons, along with an eight-metric simulation. The observational comparison reported ship rates of 23.1% under Bonferroni, 27.7% under Holm, 27.7% under Hommel, and 28.0% under BH. These are **ship rates, not unbiased estimates of statistical power**, because truth is not known in the historical experiments. Its simulation is more directly interpretable. With eight metrics, \(N=1{,}000\), effect size \(\delta=0.10\), and 10,000 replications, it reports:

| Scenario | Hommel gain over Bonferroni | BH gain over Bonferroni |
|---|---:|---:|
| 1 of 8 truly affected, independent | +0.1 percentage point | +0.8 percentage point |
| 1 of 8 truly affected, \(\rho=.95\) | +0.0 | +0.0 |
| all 8 affected, independent | +1.9 | +4.7 |
| all 8 affected, \(\rho=.95\) | +10.2 | +13.0 |

These are the publisher's reported values. **Study setting:** product experimentation, eight metrics, specified common correlation scenarios. It is close to, but smaller than, the requested 12-metric case. The strongest lesson is not a universal percentage but the interaction: gains from stepwise/FDR methods were tiny with one true effect and much larger with many simultaneous effects.

### Positive block correlation at 100 and 2,000 tests

Stevens, Al Masud, and Suyundikov simulated block-correlated positively dependent Gaussian \(Z\) statistics. Their FWER comparisons used \(m=100\), and FDR comparisons used \(m=2{,}000\), with 1,000 simulation repetitions across correlation and alternative configurations. This is larger than the requested 50-test regime but useful as a stress test for high-dimensional positive dependence. The paper also shows that not every nominally FDR-oriented method behaved satisfactorily under its dependence model, reinforcing that "more powerful" and "valid error control" must be evaluated together. **Study setting:** Gaussian block correlation, 100 tests for FWER procedures and 2,000 for FDR procedures. It is not directly comparable numerically to BH 1995 or the RCT simulations, so I do not average their results.

### Correlated binary-endpoint evidence

A published comparison specifically involving correlated binary outcomes found that Hochberg improved over Bonferroni for low to moderate positive correlations in the settings studied, while another correlation-aware procedure could outperform it under high correlation. This is qualitatively relevant to binary A/B metrics but again does not reproduce the full mixed-metric, shared-control, segment-overlap structure.

The general deterministic relationships help explain all of these simulation results:

1. Holm can never reject fewer hypotheses than Bonferroni for the same p-values, but the gain may be zero, particularly under sparse alternatives.
2. When their validity conditions hold, Hochberg's step-up logic can reject hypotheses that Holm's step-down path would not, so its gain is generally small under sparse configurations and potentially larger when several p-values are small.
3. Hommel is generally still less conservative under the Simes-valid conditions. A recent NCI review describes Hochberg as less conservative than Holm and Hommel as generally better than Hochberg. 4. BH's power advantage over FWER control is most pronounced when a material proportion of hypotheses is false, because FDR permits some false rejections as long as their expected share among rejections is controlled. This was an explicit pattern in BH's original simulations. 5. BY can sacrifice substantial power because the nominal FDR level is divided by \(H_m\). That penalty increases slowly with \(m\), approximately logarithmically, and is already a factor of about 4.5 at 50 tests.

### Evidence ledger for power claims

| Key claim | Source and location | Passage actually seen |
|---|---|---|
| BH power advantage increased with number of nonnulls and family size. | BH 1995, simulation section, discussion following tables. | "The advantage increases with the number of false null hypotheses" and the following discussion reports increasing gain with \(m\). |
| BH gain exceeded Hochberg's gain over Bonferroni in those simulations. | BH 1995, simulation discussion. | Passage seen; authors explicitly compare these two gains. |
| Bonferroni and Holm can have equal disjunctive power. | Vickerstaff et al. 2019, Results. | "The Bonferroni and Holm methods provide the same disjunctive power." |
| Hochberg/Hommel gains were small in the RCT simulations. | Vickerstaff et al. 2019, abstract/conclusions. | "The Hochberg and Hommel methods provide power gains ... albeit small" |
| Sparse-effect product simulation gave almost no gain. | Schultzberg/Spotify Confidence 2026, simulation table, lines 42-49. | "When only one metric has a real effect ... the gap collapses to near zero" |
| Dense effects gave much larger gains in that same simulation. | Same source, simulation table, lines 42-49. | Table seen: all 8 affected gives +10.2 pp Hommel and +13.0 pp BH at \(\rho=.95\). |
| High-dimensional positive-dependence study used 100 tests for FWER and 2,000 for FDR. | Stevens et al. 2017, simulation methods. | "control the FDR (when m = 2000) or FWER (when m = 100)" |

The evidence does **not** justify a claim such as "Holm gains X percent at 12 metrics" or "BH gains Y percent at 50 metrics" without specifying the number and magnitude of true effects, correlations, sample size, baseline rates, and chosen power functional. Published studies vary too much on exactly those dimensions.

## Which error rate different communities protect

**Plain-language answer.** The choice between FWER, FDR, and other criteria is driven more by the semantics and cost of a false claim than by whether there are 12 or 50 metrics. Confirmatory regulators protect against any false efficacy claim and therefore use study-wise Type I/FWER logic. Statistical literature permits FDR when the output is a collection of discoveries for which a small false fraction is acceptable. A recent review suggests FDX as an intermediate criterion for dozens to hundreds of hypotheses. Product-experimentation practice is divided because a "false positive" can mean different operational things, such as shipping one ineffective feature, investigating a false lead, or generating a false monitoring alert.

**Confidence: high** for FDA/EMA and the original statistical definitions; **medium-high** for "practitioner practice" because platform policies differ by company and can change. What would change the practical interpretation most is a precise definition of which reported cells can independently trigger a product claim or action.

### Statisticians and methodological reviews

BH's 1995 paper gives perhaps the clearest cost interpretation. It argues that FWER is natural when an overall conclusion is erroneous if even one individual inference is wrong. FDR is attractive when several discoveries are expected and a small proportion of erroneous rejections does not invalidate the scientific conclusion, as in screening or subgroup discovery. This is not a statement that large \(m\) mechanically converts FWER into FDR. It is about loss or consequence.

Formally:

\[
\mathrm{FWER}=P(V\ge1),
\]

so the cost function treats the event "at least one false claim" as the protected bad event, irrespective of whether there is one false claim or twenty. BH defines

\[
\mathrm{FDR}=E\!\left[\frac{V}{\max(R,1)}\right],
\]

so the loss is the expected **share** of claims that are false. When every null is true, FDR reduces to FWER because every rejection is false. Storey's positive FDR is

\[
\mathrm{pFDR}=E[V/R\mid R>0],
\]

which conditions on there being at least one rejection. Storey's 2002 paper argues that pFDR is often the quantity of direct interest in large multiple-testing problems and develops q-values around that perspective. A third criterion is **false discovery exceedance**, or FDX:

\[
P(\mathrm{FDP}>\gamma)\le\alpha.
\]

Instead of controlling only the expectation of the false share, FDX controls the probability that the realized false share exceeds a tolerated fraction \(\gamma\). The 2025 NCI review presents this as a useful intermediate target when the number of hypotheses is too large for stringent FWER but not large enough for the realized FDP to concentrate tightly around its expectation. That review makes an unusually explicit size-based proposal: for only a few primary endpoints, use FWER; for an intermediate number, described as **more than a dozen but less than several hundred**, consider FDX; and for very large families such as thousands of RNA-seq or GWAS hypotheses, FDR may be appropriate. Importantly, the same article says these choices should be scientifically justified and explicitly states that it is not trying to impose an inflexible method in every situation. At about 12 metrics, therefore, classical statistical guidance still often frames the question as FWER unless discoveries are genuinely screening-like. At over 50, there is more methodological support for considering FDX or FDR if false positives are inexpensive and the results are a collection to be triaged rather than individually definitive claims. The threshold itself is not a theorem.

### FDA and EMA

FDA's final 2022 guidance on multiple endpoints treats multiplicity as a threat to erroneous conclusions about treatment effects and requires prospective methods that maintain Type I error control for confirmatory endpoint claims. It discusses Bonferroni, Holm, Hochberg, and other FWER methods, noting that Holm and Hochberg can improve power over Bonferroni where their assumptions are justified. The relevant cost of a false claim is regulatory: an erroneous conclusion supporting efficacy, an endpoint claim, dose, population, or labeling conclusion. The guidance is not built around tolerating, say, an expected 5% share of false efficacy claims among all successful claims. EMA's **current effective** multiplicity document is still *Points to Consider on Multiplicity Issues in Clinical Trials*, CPMP/EWP/908/99, legally effective April 1, 2003. The 2017 document that is sometimes referred to as "the EMA multiplicity guideline" is a draft whose consultation closed and did not replace the current effective document on EMA's present page. This corrects one possible ambiguity in the source list supplied in the request.

EMA likewise frames confirmatory multiplicity around control of study-wise Type I error. Its later 2017 draft gives a particularly explicit definition of strong Type I control as protecting the probability of rejecting at least one true null regardless of which subset of nulls is true. It also makes an instructive exception in emphasis for safety: over-stringent multiplicity correction can conceal safety signals, so false negatives can be more consequential there. The 2017 wording is useful explanatory material, but it should be identified as **draft**, not current effective guidance. For regulators, therefore, the number 12 versus 50 does not itself switch the target from FWER to FDR. If all 50 endpoints are confirmatory routes to a label claim, multiplicity remains a Type I/FWER problem. What changes is that such a design can become prohibitively underpowered, which is often addressed by hierarchy, endpoint prioritization, gatekeeping, or design changes rather than by declaring that false claims are now tolerable.

### Experimentation practitioners

Practice is visibly heterogeneous.

Microsoft's Experimentation Platform reports that modern scorecards can contain hundreds or thousands of metrics. In a July 2026 article it describes moving its large-scale treatment-effect assessment toward an FDR-based approach and reports empirical A/A validation of BH in its own observed scorecard dependence structures. Earlier ExP work on automated alerts likewise discusses hundreds or thousands of metrics and BH-style FDR control as a way to limit false alerts. Here, the practical false-claim cost is often an **alert or candidate metric movement requiring attention**, not necessarily one definitive product efficacy assertion.

Statsig exposes both Bonferroni and BH options. Its documentation describes Bonferroni as controlling the chance of at least one false positive and BH as controlling the expected share of significant findings that are false positives. It allows correction across treatment groups and metric families. This is essentially a menu based on the user's loss function rather than a hard metric-count threshold.

Spotify Confidence takes a different position. In its 2026 analysis, it defines FWER operationally as the probability of at least one false positive **shipping decision** and FDR as the expected fraction of shipped features that do not actually work. It says FWER is typically the more appropriate target for product decisions, while also reporting BH as a comparator. Its argument depends heavily on defining the correction family narrowly around success metrics that can cause a ship decision, rather than including guardrails that only veto shipping. That illustrates why "how many report columns exist?" and "how many hypotheses can cause the false action?" can be very different questions. This disagreement is substantive rather than terminological:

| Decision semantics | Error criterion that directly represents the stated bad outcome | False-claim cost represented |
|---|---|---|
| Any one false efficacy/ship claim is unacceptable at the family level | FWER | Indicator cost for at least one false claim |
| A large list is being screened and some false leads are acceptable | FDR | Expected false fraction among selected claims |
| Some false discoveries are tolerable, but exceeding a fraction \(\gamma\) would be unacceptable | FDX | Tail probability of excessive false-discovery proportion |
| Only analyses with at least one discovery are relevant | pFDR | Expected false fraction conditional on making discoveries |
| Automated scorecard alerts are individually low-cost and triaged | Often FDR in practitioner systems | Expected false-alert share among alerts |
| Each significant cell can independently authorize a consequential action | FWER maps more directly to that event | Probability of at least one erroneous action |

This table states mappings between loss definitions and error criteria, not a recommendation for the requested system.

The presence of a single prespecified primary metric also matters to family definition. FDA/EMA-style confirmatory thinking often distinguishes primary, secondary, and exploratory families. The NCI review similarly argues that different families can warrant different rigor. Consequently, "one primary + 12 secondaries + dozens of segment cells" need not mathematically imply one monolithic family unless the reporting/decision policy treats every one of those cells as interchangeable routes to the same claim. Conversely, splitting a family purely to avoid multiplicity is not justified if all cells can independently trigger the same favorable conclusion.

### Key evidence on error-rate choice

| Key claim | Source and location | Passage actually seen |
|---|---|---|
| FWER is probability of at least one false rejection; FDR is expected false proportion. | BH 1995, introduction/definitions. | FWER is written as \(P\{V\ge1\}\); FDR is defined through \(E(Q)\), with \(Q=V/R\) when \(R>0\). |
| FDR is motivated when some erroneous rejections can be tolerated among many discoveries. | BH 1995, motivation after definition. | The examples include screening situations in which false leads are followed by confirmatory work. |
| A recent review proposes FWER for few, FDX for >12 to <several hundred, FDR for thousands. | Izmirlian et al. 2025, Recommendations, lines 131-141. | "more than a dozen but less than several hundred, one may wish to control FDX" |
| FDA requires suitable multiplicity methods maintaining Type I error. | FDA 2022, appendix/method-selection discussion. | "select the most powerful method that is suitable ... and maintains Type I error rate control" |
| EMA's effective document is the 2003 Points to Consider, not the 2017 draft. | EMA current guideline page, current effective version and history. | "Current effective version" followed by CPMP/EWP/908/99 and legal effective date 01/04/2003. |
| Microsoft uses an FDR-oriented approach for very large scorecards. | Microsoft Research, July 15, 2026. | "A/B tests that can generate hundreds ... or even thousands ... of metrics" and discussion of controlling FDR across a scorecard. |
| Spotify maps FWER to false shipping decisions and FDR to false shipped features. | Spotify Confidence, lines 39-40. | "FDR ... the expected fraction of shipped features that don't actually work." |
| Statsig presents both as different error semantics. | Statsig BH and Bonferroni documentation. | BH controls a share of significant results that are false; Bonferroni controls chance of at least one false positive. |

## Source table and access record

The table distinguishes what I could actually inspect from what was available only as an abstract, search extract, or inaccessible full article. "Full text" means the substantive article or official guidance text was readable through the research interface. For PDFs, I used rendered PDF access where the tool allowed it; an attempted screenshot of the EMA 2003 PDF endpoint was not supported because that endpoint was exposed as HTML rather than `application/pdf`.

| Full reference | Working link | What I actually read |
|---|---|---|
| Holm, Sture. "A Simple Sequentially Rejective Multiple Test Procedure." *Scandinavian Journal of Statistics* 6(2), 65-70, 1979. | | Full PDF text; rendered page inspected. |
| Sidak, Zbynek. "Rectangular Confidence Regions for the Means of Multivariate Normal Distributions." *Journal of the American Statistical Association* 62, 626-633, 1967. | | Full PDF text; rendered page inspected. |
| Simes, R. J. "An Improved Bonferroni Procedure for Multiple Tests of Significance." *Biometrika* 73(3), 751-754, 1986. | | Full PDF text. Screenshot rendering was attempted but timed out. |
| Hochberg, Yosef. "A Sharper Bonferroni Procedure for Multiple Tests of Significance." *Biometrika* 75(4), 800-802, 1988. | | Abstract/bibliographic page only. Full text could not be opened without publisher access. |
| Hommel, Gerhard. "A Stagewise Rejective Multiple Test Procedure Based on a Modified Bonferroni Test." *Biometrika* 75(2), 383-386, 1988. | | Abstract/bibliographic material only. Full article could not be opened. |
| Benjamini, Yoav, and Yosef Hochberg. "Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing." *Journal of the Royal Statistical Society, Series B* 57(1), 289-300, 1995. | | Full PDF text; rendered page inspected. |
| Benjamini, Yoav, and Daniel Yekutieli. "The Control of the False Discovery Rate in Multiple Testing under Dependency." *Annals of Statistics* 29(4), 1165-1188, 2001. | | Full preprint/PDF text; rendered pages inspected. |
| Storey, John D. "A Direct Approach to False Discovery Rates." *Journal of the Royal Statistical Society, Series B* 64(3), 479-498, 2002. | | Full PDF text; rendered page inspected. |
| Storey, John D., Jonathan E. Taylor, and David Siegmund. "Strong Control, Conservative Point Estimation and Simultaneous Conservative Consistency of False Discovery Rates: A Unified Approach." *Journal of the Royal Statistical Society, Series B* 66(1), 187-205, 2004. | | Full PDF summary/open preprint and publisher abstract. |
| Sarkar, Sanat K. "Some Probability Inequalities for Ordered MTP2 Random Variables: A Proof of the Simes Conjecture." *Annals of Statistics* 26(2), 494-504, 1998. | | Full text exposed through web extraction. |
| Sarkar, Sanat K., and Chung-Kuei Chang. "The Simes Method for Multiple Hypothesis Testing with Positively Dependent Test Statistics." *Journal of the American Statistical Association* 92(440), 1601-1608, 1997. | | Abstract only. Full article endpoint returned access failure. |
| Dunnett, Charles W. "A Multiple Comparison Procedure for Comparing Several Treatments with a Control." *Journal of the American Statistical Association* 50(272), 1096-1121, 1955. | | Bibliographic/landing-page material only; original full text could not be opened. |
| Block, Henry W., Thomas H. Savits, and Jian Wang. "Negative Dependence and the Simes Inequality." *Journal of Statistical Planning and Inference*, 2008. | | Abstract/search description only; full article could not be opened. I did not fill in page numbers or DOI because I did not verify them from an opened source. |
| Dobriban, Edgar. "The Benjamini-Hochberg Procedure Can Fail to Control the FDR for Correlated Two-sided Gaussian Tests." Published electronically October 5, 2026; print October 13, 2026. | | Extensive full-text web extraction including abstract, model, theorem, finite \(m=85\) note, and discussion. Direct PMC open subsequently hit a browser challenge. |
| Stevens, John R., Abdullah Al Masud, and Anvar Suyundikov. "A Comparison of Multiple Testing Adjustment Methods with Block-correlation Positively-dependent Tests." *PLoS ONE*, 2017. | | Full/open article and printable-text extracts. |
| Vickerstaff, Victoria, et al. "Methods to Adjust for Multiple Comparisons in the Analysis and Sample Size Calculation of Randomised Controlled Trials with Multiple Primary Outcomes." *BMC Medical Research Methodology* 19:129, 2019. | | Full open web article/search extraction and abstract; numerical table text visible. |
| Izmirlian, Grant, Lev A. Sirota, Vance W. Berger, and Victor Kipnis. "The Fundamentals of Multiplicity Adjustment in Biostatistics." *Journal of the National Cancer Institute Monographs* 2025(68), 10-13, 2025. DOI 10.1093/jncimonographs/lgae050. | | Full PMC text. |
| U.S. Food and Drug Administration. *Multiple Endpoints in Clinical Trials: Guidance for Industry*. Final guidance, 2022. | | Official FDA PDF/text, including multiplicity rationale and appendix methods. |
| European Agency for the Evaluation of Medicinal Products / CPMP. *Points to Consider on Multiplicity Issues in Clinical Trials*, CPMP/EWP/908/99, adopted September 2002, legally effective April 1, 2003. | | Official EMA current-version page and PDF text. The endpoint did not expose itself to the screenshot tool as a PDF. |
| European Medicines Agency. *Draft Guideline on Multiplicity Issues in Clinical Trials*, EMA/CHMP/44762/2017, consultation 2017. | | Full draft PDF text and rendered pages. Not current effective guidance. |
| Qi, Kai, and Momo Jeng. "Treatment Effect Assessment at Scale: Accounting for Correlated Metrics and Metric Relevance in Modern Experimentation." Microsoft Research / Experimentation Platform, July 15, 2026. | | Full web article/search extraction. |
| Agrawal, Ankita, Vaibhav Kumar Bajpai, and Jiannan Lu. "Alerting in Microsoft's Experimentation Platform (ExP)." Microsoft Research. | | Full web article/search extraction. |
| Schultzberg, Marten. "Are Optimal Multiple Testing Corrections Optimal for You?" Confidence by Spotify, April 14, 2026, updated May 12, 2026. | | Full web page, including historical experiment table and simulation table. |
| Statsig. "Bonferroni Correction." Statsig Documentation, accessed November 10, 2026. | | Current documentation text. |
| Statsig. "Benjamini-Hochberg Procedure." Statsig Documentation, accessed November 10, 2026. | | Current documentation text. |
| Comparative literature on multiple correlated binary endpoints. | | Abstract/full web descriptions as returned by indexing services; used only for limited qualitative claims. |

### Sources attempted but not fully open

Every source I attempted and could not inspect in full is listed here. Hochberg 1988 and Hommel 1988 were available only through abstracts/bibliographic pages; I did not bypass the publisher access controls. Sarkar and Chang 1997 returned an abstract but the attempted full-text endpoint returned HTTP 403. Block, Savits, and Wang's negative-dependence article likewise yielded abstract/index text but not an open full article. Dunnett 1955 was available as a bibliographic/publisher landing page but not open full text. The Dobriban 2026 article was initially available through a detailed full-text index extraction, but a subsequent direct PMC open produced a browser/reCAPTCHA challenge; I used only the content that had already been exposed by the research service and did not bypass the challenge. The Simes 1986 PDF itself opened and its text was readable, but an attempted rendered-page screenshot timed out. The EMA 2003 official document text was readable through EMA's endpoint, but the screenshot tool reported that the endpoint was not presented as `application/pdf`; no attempt was made to circumvent that limitation. These are tool-rendering limitations rather than paywall access.
