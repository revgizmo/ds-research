# Confidence intervals alongside multiplicity-adjusted p-values in fixed-horizon A/B testing

## Scope, assumptions, and main findings

This report assumes a single fixed-horizon randomized experiment, with one shared control and one to three treatment variations, two-sided frequentist inference, one primary metric, many secondary metrics, and potentially many overlapping segment cells. It assumes that the inferential target for a binary metric is some treatment-control contrast such as a risk difference, risk ratio, or odds ratio, but it does not assume one particular effect scale. Exact interval details for binary outcomes depend materially on that choice of effect scale and on the underlying test.

I also assume that a "family" of hypotheses has been defined before multiplicity adjustment. The sources reviewed do not settle whether all metrics, variations, segments, and segment combinations in a large product report should constitute one family. That question is consequential, especially when one report contains dozens or hundreds of overlapping cells, but it is separate from how confidence intervals should be generated once the family is defined.

The central findings are:

- Compatible simultaneous intervals do exist for Bonferroni, Holm, Hochberg, Hommel, and substantial classes of closed-testing procedures. Bonferroni is straightforward. Holm and more elaborate procedures require nontrivial inversion or confidence-region constructions, and exact compatibility can come at the price of intervals that are surprisingly uninformative. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/), [2](https://pubmed.ncbi.nlm.nih.gov/18618415/), [3](https://pubmed.ncbi.nlm.nih.gov/22685000/)]
- Benjamini and Yekutieli's 2005 FCR construction is a selective-inference method, not a simultaneous-coverage method. Its general construction uses a selection-dependent quantity \(R_{\min}\), although for several common selectors, including unadjusted testing, Bonferroni selection, and BH selection, it simplifies to using the observed number selected \(R\). Under independent components it controls FCR for any selection rule; under positive dependence additional concordance conditions enter; under arbitrary dependence a harmonic-factor correction is available. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]
- Published methods do satisfy the user's requested decision compatibility: at a fixed alpha, a hypothesis is rejected if and only if its compatible confidence interval excludes the null set. Since an adjusted p-value is the smallest familywise alpha at which that hypothesis is rejected, this gives the desired equivalence with adjusted \(p<\alpha\), subject to equality conventions and complications from discrete tests. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/), [5](https://pmc.ncbi.nlm.nih.gov/articles/PMC4564263)]
- Regulatory and journal guidance does not establish a universal rule that every multiplicity-adjusted p-value must be paired with a multiplicity-adjusted CI. ICH explicitly calls for confidence intervals and consideration of multiplicity in significance and confidence levels. FDA emphasizes prespecified familywise control and discusses resampling but is less explicit about matching CI adjustment. Current JAMA and NEJM policies permit or require unadjusted 95% intervals in some exploratory settings if their multiplicity status is made clear and they are not interpreted as formal hypothesis tests. The strongest explicit statement favoring conservative adjusted CIs when matching intervals are unavailable occurs in the EMA's 2017 draft guideline, which is not the current effective EMA guideline. [[6](https://database.ich.org/sites/default/files/E9_Guideline.pdf), [7](https://jamanetwork.com/journals/jama/pages/instructions-for-authors), [8](https://www.ema.europa.eu/en/scientific-guidelines/multiplicity-issues-clinical-trials), [9](https://www.fda.gov/media/162416/download?action=ReferProject), [10](https://www.nejm.org/author-center/new-manuscripts)]
- Westfall-Young style resampling can exploit dependence among endpoints and shared-control contrasts. Permutation versions can be finite-sample exact when the required randomization or exchangeability invariance is valid and the reference distribution is handled correctly. Bootstrap versions are ordinarily approximate and asymptotic. A bootstrap CI is not obtained merely by converting one adjusted p-value into a confidence level. It is either constructed at a multiplicity-adjusted level or, for exact test-CI compatibility, obtained by inverting the relevant family of multiplicity-adjusted tests. [[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/), [9](https://www.fda.gov/media/162416/download?action=ReferProject), [12](https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf)]

The overlapping-cell structure matters particularly for resampling. If unit-level outcome vectors and treatment assignments are available, one can resample or permute whole units so that cross-metric, shared-control, and overlapping-segment dependence is preserved. If the system retains only separate aggregate counts for each metric-cell pair, those marginal counts generally do not identify the joint dependence required by a joint maxT or minP resampling distribution. That last statement is a statistical identifiability consequence rather than a claim made verbatim in one of the cited papers.

## Simultaneous intervals for Bonferroni, Holm, Hochberg, and closed testing

**Plain-language answer.** Yes. Bonferroni simultaneous two-sided intervals are immediate. Compatible intervals for Holm and broad Bonferroni-based closed tests were published by Strassburger and Bretz and Guilbaud. Guilbaud later extended the construction to closed tests corresponding to Holm, Hochberg, and Hommel procedures. These methods can give exact agreement between rejection decisions and interval exclusion, but there is a genuine information trade-off: for powerful stepwise procedures the compatible intervals can provide almost no quantitative information beyond the reject/do-not-reject decision unless all hypotheses are rejected. [[2](https://pubmed.ncbi.nlm.nih.gov/18618415/), [13](https://pubmed.ncbi.nlm.nih.gov/18932131/), [3](https://pubmed.ncbi.nlm.nih.gov/22685000/)]

**Confidence: high** for existence, construction principles, one-sided Holm/Hochberg/Hommel results, and the informativeness problem. **Confidence: medium** for a blanket claim that every published one-sided construction has an equally simple explicit two-sided rectangular analogue. General two-sided inversion is available in principle, and Hayter and Hsu explicitly studied a two-sided step-down example, but the accessible abstracts of the principal Holm/Hochberg papers emphasize one-sided applications. Full access to those paywalled articles could raise that part to high confidence. [[3](https://pubmed.ncbi.nlm.nih.gov/22685000/), [14](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476453)]

For Bonferroni, suppose \(C_i(\gamma)\) is a valid marginal \(1-\gamma\) two-sided interval for parameter \(\theta_i\), for \(i=1,\ldots,m\). Then

\[
C_i(\alpha/m)
\]

for every \(i\) yields at least \(1-\alpha\) simultaneous coverage by the union bound, regardless of dependence among the tests. Weighted Bonferroni replaces \(\alpha/m\) by prespecified \(\alpha_i\) satisfying \(\sum_i\alpha_i\le\alpha\). If each marginal interval is obtained by inverting the same marginal two-sided test used to compute the p-value, null exclusion and Bonferroni rejection agree at the chosen alpha.

Holm is not obtained safely by simply looking at the observed p-value rank and assigning that hypothesis a post hoc marginal interval at \(1-\alpha/(m-r+1)\). Strassburger and Bretz specifically note that practical ad hoc attempts can fail to control the specified simultaneous confidence level. Their construction works through the closed-testing representation of Holm, using weighted Bonferroni tests for the intersection hypotheses. [[2](https://pubmed.ncbi.nlm.nih.gov/18618415/)]

A general way to understand the construction is inversion of the multiple-testing procedure. For every candidate parameter vector \(\vartheta\), test the corresponding family of shifted null hypotheses at overall level \(\alpha\). The set of parameter vectors not rejected is a \(1-\alpha\) confidence region. Projecting that region on each coordinate produces simultaneous marginal intervals, although projection can enlarge them considerably. A recent paper by Brannath, Kluge, and Scharpenberg states this construction directly and notes that projection of a joint confidence set onto coordinate axes retains at least the joint coverage probability. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/)]

Hochberg is more delicate because it is a step-up procedure rather than Holm's step-down procedure. Guilbaud's 2012 paper explicitly applies a general closed-testing confidence-region construction to procedures related to Holm, Hochberg, and Hommel. The accessible abstract describes the application as one-sided. [[3](https://pubmed.ncbi.nlm.nih.gov/22685000/)]

The informativeness issue is not merely computational. For the Holm, Hochberg, and Hommel constructions studied by Guilbaud, the confidence region can be such that no assertion sharper than the rejection assertion is available unless all null hypotheses are rejected. A rejected one-sided null such as \(\theta_i\le0\) may therefore yield essentially only the statement \(\theta_i>0\), rather than a useful lower bound appreciably above zero. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/), [3](https://pubmed.ncbi.nlm.nih.gov/22685000/)]

More recent work treats this as an explicit trade-off. Brannath and colleagues distinguish "compatible" intervals, which reproduce the original multiple-test decisions exactly, from "informative" intervals, which continue to move as evidence becomes stronger. Their newer graphical procedures can improve informativeness, but in general do not reproduce every rejection of the original, more powerful graphical test. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/)]

Two-sided intervals are entirely straightforward for Bonferroni. Two-sided stepwise confidence constructions also exist in the literature: Hayter and Hsu's 1994 paper includes a two-sided step-down example for two parameters. In a general closed-testing framework, two-sided intervals can be obtained by inverting two-sided shifted tests, but the exact rectangular form and its usefulness depend on the local tests and closure. [[15](https://pmc.ncbi.nlm.nih.gov/articles/PMC4806862/), [14](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476453)]

For binary outcomes there is an additional test-CI matching issue independent of multiplicity. Exact two-sided tests for discrete data do not automatically match every familiar "exact" interval. A confidence interval must be the inversion of the same testing rule if exact null-exclusion equivalence is required; otherwise a discrete p-value can cross alpha while a separately chosen interval still contains the null, or conversely. [[16](https://journal.r-project.org/articles/RJ-2010-008/)]

### Evidence for the key claims

| Key claim | Source and location | Verbatim passage seen |
|---|---|---|
| Compatible Holm intervals exist, and naive alternatives can fail simultaneous coverage. | Strassburger and Bretz 2008, PubMed abstract. [[2](https://pubmed.ncbi.nlm.nih.gov/18618415/)] | "ad hoc solutions used in practice typically do not control the pre-specified simultaneous confidence level" |
| The class includes Holm and other weighted Bonferroni closed tests. | Strassburger and Bretz 2008, abstract. [[2](https://pubmed.ncbi.nlm.nih.gov/18618415/)] | "derive simultaneous confidence intervals that are compatible with a certain class of closed test procedures" |
| Closed-test confidence regions cover Hochberg and Hommel as well as Holm. | Guilbaud 2012, abstract, lines 238-240. [[3](https://pubmed.ncbi.nlm.nih.gov/22685000/)] | "It includes Holm's, Hochberg's, and Hommel's MTPs." |
| Compatibility can produce intervals with little information beyond rejection. | Guilbaud 2012, abstract, lines 238-240. [[3](https://pubmed.ncbi.nlm.nih.gov/22685000/)] | "no confidence assertions sharper than rejection assertions can be made unless all null hypotheses are rejected" |
| Compatibility means rejection iff interval exclusion. | Brannath, Kluge, Scharpenberg 2025, Section 1.1, line 87. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/)] | "a null hypothesis is rejected if and only if it is excluded from the confidence interval" |
| Projection gives simultaneous coverage. | Brannath et al. 2025, Section 2.1, lines 190-196. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/)] | "Since it contains C, its coverage probability is not smaller than 1-alpha." |
| Two-sided step-down confidence sets have been studied. | Hayter and Hsu 1994, abstract. [[14](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476453)] | "for a two-sided step-down decision procedure" |

## False coverage rate intervals after selection

**Plain-language answer.** Benjamini and Yekutieli's 2005 FCR method controls the expected fraction of reported intervals that fail to cover their parameters. For arbitrary selection it is not, in full generality, simply "select \(R\) parameters and report \(1-Rq/m\) intervals." The paper's general construction uses, for each selected parameter, a quantity \(R_{\min}\) describing how few parameters could have been selected while that parameter remained selected as its own statistic varied. For several common selectors \(R_{\min}\) equals the observed \(R\), which yields the familiar \(1-Rq/m\) intervals. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

**Confidence: high.** The full paper was accessible and its definitions and theorems were checked directly. Confidence would change only if a later correction to the 2005 paper materially changed these results; none was identified in the sources reviewed.

Let \(S(T)\) be the data-dependent set of selected parameters and \(R=|S(T)|\). If \(V\) of the \(R\) reported intervals fail to cover, the false coverage proportion is conventionally taken as \(V/R\) when \(R>0\), with zero contribution when none are reported. FCR is the expectation of that false coverage proportion. This is analogous to FDR, but the "false statements" are noncovering confidence intervals rather than rejected true null hypotheses. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

The paper assumes a nested family of marginal intervals \(CI_i(\gamma)\) available at any desired noncoverage level \(\gamma\), with marginal coverage at least \(1-\gamma\). The nesting assumption is that allowing more noncoverage produces a subset of the more conservative interval. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

The exact general construction is:

1. Apply the selection rule \(S(T)\).
2. For each selected \(i\), hold all coordinates other than \(T_i\) fixed and determine \(R_{\min}(T^{(i)})\), the smallest number of selections that can occur over values of \(T_i\) for which \(i\) itself remains selected.
3. Report \(CI_i\!\left(R_{\min}(T^{(i)})q/m\right)\).

Under independent components of \(T\), Theorem 1 states that this controls FCR at \(q\) for **any** selection procedure. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

For common monotone selection mechanisms the construction simplifies. The paper specifically identifies selection by unadjusted tests, Bonferroni tests, and BH tests as cases where \(R_{\min}\) can be replaced by the observed number selected \(R\). Then every selected parameter gets a marginal interval with noncoverage level

\[
\gamma = \frac{Rq}{m},
\]

or equivalently coverage \(1-Rq/m\). If one parameter is selected from \(m\), this is Bonferroni-like. If all \(m\) are selected, the intervals revert to ordinary \(1-q\) marginal intervals. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

For the special BH-linked construction, Definition 1 sorts the p-values, selects \(R\) according to the BH criterion, and then constructs a \(1-Rq/m\) interval for every selected parameter. The paper explicitly treats symmetric confidence intervals selected by two-sided testing and obtains FCR control. It also shows that the BH-linked procedure is less overconservative in some settings than selecting with Bonferroni and then applying FCR adjustment. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

The framework does extend beyond BH selection. Under independence, Theorem 1 is stated for any selection procedure. Under positive regression dependence, the paper gives a further result requiring the interval procedure and selection rule to satisfy a concordance condition. Under completely arbitrary dependence, Theorem 4 bounds the FCR by \(qH_m\), where

\[
H_m=\sum_{j=1}^m \frac1j,
\]

so using \(q/H_m\) in the FCR construction restores a \(q\) bound. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

The guarantee is neither ordinary per-interval conditional coverage after selection nor simultaneous coverage. It controls an expectation over the proportion of failed intervals among those reported. A selected interval can have poor conditional coverage in a particular selection event while the aggregate FCR criterion is still controlled. This distinction is a major reason not to label FCR-adjusted intervals as "95% simultaneous CIs." The BY paper explicitly distinguishes multiplicity adjustment for simultaneous coverage from adjustment for selection. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)]

This matters for the proposed A/B system. If every metric and segment is always displayed, there is no selection event of the kind FCR was designed to address. If the system only surfaces "interesting" cells or metrics based on their observed statistics, FCR becomes directly relevant. That is an application-level inference from the BY framework, not a recommendation that selection be used.

### Evidence for the key claims

| Key claim | Source and location | Verbatim passage seen |
|---|---|---|
| General construction assumes valid nested marginal CIs at arbitrary levels. | Benjamini and Yekutieli 2005, Section 3, Definition 3 setup, p. 75, lines 388-397. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)] | "We assume that we have at our disposal a procedure for constructing marginal CIs at any desired level." |
| Simple \(R\) substitution covers unadjusted, Bonferroni, and BH selection. | BY 2005, Remark 1, p. 75, lines 456-459. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)] | "including selection by unadjusted testing, by Bonferroni testing, and by BH testing" |
| Under independence the procedure works for any selector. | BY 2005, Theorem 1, p. 75, lines 548-550. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)] | "for any selection procedure S(T)" |
| The simplified level is \(Rq/m\). | BY 2005, p. 75, lines 534-543. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)] | "multiply q by the number of parameters selected, divide by the size of the pool" |
| Arbitrary dependence costs the harmonic factor. | BY 2005, Theorem 4, p. 79, lines 1285-1305. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)] | "any dependency structure of the test statistics" |
| Two-sided BH-selected intervals are explicitly treated. | BY 2005, Section 1 and Definition 1, p. 73, lines 259-277. [[4](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf)] | "construction of symmetric CIs for parameters selected by two-sided multiple-hypothesis testing procedures" |

## Exact agreement between adjusted p-values and interval exclusion

**Plain-language answer.** Yes, methods with this property are published. The literature calls them compatible simultaneous confidence intervals or confidence regions. At a specified familywise level \(\alpha\), compatibility means exactly that \(H_i\) is rejected if and only if the interval for parameter \(i\) excludes the entire null set. Because an adjusted p-value is the smallest familywise alpha at which \(H_i\) is rejected, this gives exclusion exactly when the adjusted p-value crosses alpha, apart from the choice of strict versus non-strict inequality and discrete-test boundary effects. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/), [5](https://pmc.ncbi.nlm.nih.gov/articles/PMC4564263)]

**Confidence: high** for Holm and broad closed-testing classes, and high that the concept explicitly exists in published literature. **Confidence: medium-high** for every conceivable implementation of Hochberg with discrete binary tests, because exact equality can depend on how ties, p-value discreteness, and inversion boundaries are defined. Full text of Guilbaud 2012 would improve the implementation-level detail.

The clearest verified examples are Strassburger and Bretz for weighted Bonferroni closed tests including Holm; Guilbaud for Holm and related closed-testing procedures; Guilbaud's later extension covering Holm-, Hochberg-, and Hommel-related procedures; and general closed-test inversion constructions. [[2](https://pubmed.ncbi.nlm.nih.gov/18618415/), [13](https://pubmed.ncbi.nlm.nih.gov/18932131/), [3](https://pubmed.ncbi.nlm.nih.gov/22685000/)]

For a closed testing procedure, an elementary null \(H_i\) is rejected at level \(\alpha\) only if every intersection hypothesis containing \(H_i\) is rejected at that level. Its closure-adjusted p-value is the maximum of the relevant intersection-test p-values. Thus

\[
p_{i,\mathrm{adj}}\le\alpha
\quad\Longleftrightarrow\quad
H_i\text{ is rejected by the closed test at }\alpha.
\]

A compatible confidence set is obtained by inverting those tests over candidate parameter values. Therefore

\[
p_{i,\mathrm{adj}}\le\alpha
\quad\Longleftrightarrow\quad
H_i\cap CI_i(\alpha)=\varnothing
\]

when the interval construction is exactly compatible with that test and the same inequality convention is used. The closed-testing review gives the adjusted-p definition as the maximum local p-value among intersections containing \(i\). [[5](https://pmc.ncbi.nlm.nih.gov/articles/PMC4564263)]

The relevant coverage guarantee in this literature is simultaneous coverage. If the rectangular confidence region has coverage at least \(1-\alpha\), then every individual interval also has unconditional marginal coverage at least \(1-\alpha\), but that marginal statement is weaker than the simultaneous guarantee. If all intervals are constructed before any selection, simultaneous coverage also implies an FCR no greater than alpha for any subset subsequently reported, because the event of a reported noncovering interval is contained in the event that at least one member of the full simultaneous family misses. That FCR statement is a mathematical consequence, not the design target of the compatible-SCI papers.

This does **not** mean that a useful interval can be obtained from the numerical adjusted p-value alone. A single adjusted p-value describes the test at one null value. A confidence interval needs information about the behavior of the test over a continuum of candidate parameter values. Inversion therefore requires either the full family of shifted tests or an equivalent pivotal/parametric representation. This distinction is especially important for stepwise procedures because the effective local level can change as the candidate parameter value changes.

The main negative result is practical rather than existential. Exact compatibility can make intervals effectively collapse to the logical assertion produced by the test. The recent informative-SCI literature emphasizes that, for some stepwise and graphical procedures, demanding both maximal testing power and exact compatible quantitative bounds creates a conflict. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/), [3](https://pubmed.ncbi.nlm.nih.gov/22685000/)]

For discrete binary outcomes, "exactly when" also requires consistency about the boundary. A rule using \(p<\alpha\) need not be identical to one using \(p\le\alpha\). Discrete attainable p-values can land exactly on a threshold with nonzero probability. Moreover, an "exact" CI from one construction need not match an exact test from another construction. The discrete-data literature explicitly documents such test-CI inconsistencies. [[16](https://journal.r-project.org/articles/RJ-2010-008/)]

### Evidence for the key claims

| Key claim | Source and location | Verbatim passage seen |
|---|---|---|
| Compatibility is an iff property, not merely one-way consistency. | Brannath et al. 2025, Section 1.1, line 87. [[1](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/)] | "rejected if and only if it is excluded from the confidence interval" |
| Closure-adjusted p is the maximum intersection p-value. | Westfall et al. review, Section 3, "Overview of Closed Testing." [[5](https://pmc.ncbi.nlm.nih.gov/articles/PMC4564263)] | "The closure adjusted p-value for an individual hypothesis Hi is simply" |
| Holm-compatible simultaneous CIs are published. | Strassburger and Bretz 2008, abstract. [[2](https://pubmed.ncbi.nlm.nih.gov/18618415/)] | "simultaneous confidence intervals that are compatible" |
| Hochberg-compatible closed-test regions are published. | Guilbaud 2012, abstract. [[3](https://pubmed.ncbi.nlm.nih.gov/22685000/)] | "It includes Holm's, Hochberg's, and Hommel's MTPs." |
| Exact two-sided test-CI matching can fail for discrete data if separately constructed. | Fay 2010, R Journal, introduction. [[16](https://journal.r-project.org/articles/RJ-2010-008/)] | "a test-CI inconsistency may occur" |

## Regulators, journals, and reporting guidance

**Plain-language answer.** The guidance is more nuanced than "adjust p-values, therefore adjust every CI." ICH E9 asks for effect estimates with confidence intervals and has a section explicitly titled "Adjustment of Significance and Confidence Levels." FDA expects multiplicity control for confirmatory endpoint families, but its 2022 multiple-endpoints guidance does not provide a general requirement that the displayed CI must be the inversion of the chosen adjusted test. Current journal policies often permit ordinary 95% CIs for exploratory outcomes provided their multiplicity status is stated and they are not treated as formal significance tests. An EMA 2017 draft gave the most explicit advice found here: if compatible intervals are unavailable, simple conservative intervals such as Bonferroni intervals can be used, and decisions should rest on the multiplicity-adjusted test. That document remains a draft rather than the current effective EMA guidance. [[6](https://database.ich.org/sites/default/files/E9_Guideline.pdf), [7](https://jamanetwork.com/journals/jama/pages/instructions-for-authors), [8](https://www.ema.europa.eu/en/scientific-guidelines/multiplicity-issues-clinical-trials), [10](https://www.nejm.org/author-center/new-manuscripts)]

**Confidence: high** on ICH, FDA, JAMA, current EMA document status, and CONSORT. **Confidence: medium** on the exact current NEJM wording because the author-center page was visible in the search extract but a direct fetch returned HTTP 403. **Confidence: medium** on attributing the 2017 EMA interval recommendation to present EMA policy, which I do not do here: it is explicitly labeled draft.

ICH E9 states that effect estimates should be accompanied by confidence intervals whenever possible, that two-sided intervals are generally appropriate for estimating treatment differences, and that when multiplicity remains in confirmatory analyses adjustment should be considered and prespecified. Its Section 5.6 explicitly couples "significance and confidence levels" in the section title, although the body is more explicit about Type I error than about a formula for adjusted CIs. [[6](https://database.ich.org/sites/default/files/E9_Guideline.pdf)]

FDA's October 2022 final guidance, *Multiple Endpoints in Clinical Trials*, emphasizes prospectively specified multiplicity procedures and overall Type I error control. Its resampling section recognizes Westfall-Young style methods as a way to use endpoint correlation while controlling overall Type I error, but cautions that bootstrap and permutation methods rely on important assumptions and may require simulation checks. I did not find in the final document a general instruction saying that every adjusted p-value must be accompanied by a confidence interval adjusted with the identical procedure. [[9](https://www.fda.gov/media/162416/download?action=ReferProject)]

The EMA status requires care. The EMA web page currently identifies CPMP/EWP/908/99, effective April 1, 2003, as the **current effective version**. It separately lists the 2017 replacement guideline as a consultation draft whose consultation is closed. [[8](https://www.ema.europa.eu/en/scientific-guidelines/multiplicity-issues-clinical-trials)] The 2017 draft discusses the problem of nonadjusted intervals having lower-than-nominal simultaneous coverage, acknowledges that confidence regions corresponding to complex multiplicity procedures may be difficult or uninformative, and advises use of simple conservative intervals such as Bonferroni-corrected intervals when needed. Because this text is draft guidance, it should not be represented as binding or current EMA guidance.

NEJM's current author guidance, as returned by its indexed page, says multiplicity-adjusted p-values should be labeled as such for confirmatory analyses. Where no multiplicity method was prespecified, it says secondary and exploratory outcomes should generally be presented as point estimates with 95% CIs, with an explicit statement that the interval widths were not adjusted for multiplicity and should not substitute for hypothesis testing. [[10](https://www.nejm.org/author-center/new-manuscripts)] This is a particularly direct example of a journal explicitly allowing unadjusted intervals next to a multiplicity-aware interpretation, provided their limitations are stated.

JAMA's current author instructions require quantitative estimates and measures of uncertainty rather than p-values alone. They also say that secondary and subgroup analyses should describe how Type I error from multiple comparisons was handled, or else be regarded as exploratory. JAMA does not, in the passage reviewed, state that the confidence levels themselves must always be adjusted. [[7](https://jamanetwork.com/journals/jama/pages/instructions-for-authors)]

CONSORT-Outcomes 2022 similarly emphasizes disclosure: its item 12a.1 asks authors to describe methods used to account for multiplicity in the analysis or interpretation of primary and secondary outcomes. It is a reporting standard, not a prescription for one interval construction. [[17](https://jamanetwork.com/journals/jama/fullarticle/2799401)]

The resulting reporting trade-off is therefore well documented:

- **Multiplicity-adjusted simultaneous intervals** align inferential error control across the whole family, but may be much wider and, under exact closed-test compatibility, can be quantitatively uninformative.
- **Ordinary marginal intervals plus adjusted p-values** retain interpretable effect-size precision for each estimand, but interval exclusion of zero must not be read as the multiplicity-controlled testing decision.
- **Compatible closed-test intervals** make the p-value and CI decision agree exactly, but the literature shows that exact compatibility can sacrifice informativeness.
- **FCR intervals** address data-dependent selection of what gets reported, which is a different target from simultaneous familywise coverage.

Those are trade-offs, not a recommendation among them.

### Evidence for the key claims

| Key claim | Source and location | Verbatim passage seen |
|---|---|---|
| ICH expects estimates plus CIs. | ICH E9, Section 5.5, p. 25, lines 566-578. [[6](https://database.ich.org/sites/default/files/E9_Guideline.pdf)] | "Estimates of treatment effects should be accompanied by confidence intervals, whenever possible" |
| ICH explicitly addresses multiplicity under "Adjustment of Significance and Confidence Levels." | ICH E9, Section 5.6, p. 26, lines 581-584. [[6](https://database.ich.org/sites/default/files/E9_Guideline.pdf)] | "When multiplicity is present... may necessitate an adjustment to the type I error." |
| EMA's current effective document is the 2003 Points to Consider, not the 2017 draft. | EMA multiplicity guideline page, lines 24-44. [[8](https://www.ema.europa.eu/en/scientific-guidelines/multiplicity-issues-clinical-trials)] | "Current effective version" |
| FDA recognizes multiplicity resampling but warns about assumptions. | FDA 2022, Appendix, Section 6, p. 21, lines 252-269. [[9](https://www.fda.gov/media/162416/download?action=ReferProject)] | "require few, albeit important, assumptions about the true distribution of the endpoints" |
| JAMA requires uncertainty measures and disclosure of multiple-comparison handling. | JAMA Instructions for Authors, lines 2961-2965. [[7](https://jamanetwork.com/journals/jama/pages/instructions-for-authors)] | "there should be a description of how the potential for type I error due to multiple comparisons was handled" |
| CONSORT-Outcomes requires reporting how multiplicity was handled. | CONSORT-Outcomes 2022, checklist item 12a.1. [[17](https://jamanetwork.com/journals/jama/fullarticle/2799401)] | "Describe any methods used to account for multiplicity in the analysis or interpretation" |
| NEJM permits unadjusted intervals when multiplicity was not prespecified, with explicit caveats. | NEJM Author Center, indexed current guidance. [[10](https://www.nejm.org/author-center/new-manuscripts)] | "the widths of the intervals have not been adjusted for multiplicity" |

## Resampling, bootstrap intervals, binary counts, Monte Carlo error, and reproducibility

**Plain-language answer.** Westfall-Young style resampling is a familywise multiple-testing framework that estimates a joint null distribution of test statistics and uses maxima, minima, or step-down versions to adjust for multiplicity while retaining dependence. Permutation implementations can be finite-sample exact when the null implies the required randomization or permutation invariance. Bootstrap implementations are generally approximate. Classical computational shortcuts often invoke subset pivotality, although later work shows that full closed testing can obtain strong FWER control under much weaker assumptions if each intersection test is valid. A bootstrap confidence interval is generated from a bootstrap distribution, an adjusted confidence level, or inversion of a bootstrap test. It is not determined by the adjusted p-value at the null alone. [[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/), [9](https://www.fda.gov/media/162416/download?action=ReferProject), [12](https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf)]

**Confidence: high** on the distinction between permutation exactness and bootstrap approximation, the role of subset pivotality, bootstrap interval construction, and Monte Carlo calculations. **Confidence: medium-high** for finite-sample exactness in every A/B implementation because that depends on the exact randomization scheme, null hypothesis, statistic, clustering, and whether permutation invariance actually holds. A concrete data-generating and randomization specification would determine that, but the user asked not to narrow the setting further.

### Westfall-Young and subset pivotality

Westfall and Young's 1993 framework uses resampling to estimate the joint null distribution of a vector of statistics and adjusts for the multiplicity of their extrema. Because the same resampling draw generates all test statistics together, correlation can reduce the penalty relative to Bonferroni. FDA explicitly notes this power advantage for correlated endpoints. [[9](https://www.fda.gov/media/162416/download?action=ReferProject)]

The usual maxT computational shortcut evaluates the distribution of the maximum over a subset under the complete null rather than re-running a separate resampling test for every intersection null. Westfall and Troendle define subset pivotality as the condition that the relevant subset-max statistic has the same null distribution under its own intersection null as under the complete null. [[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/)]

There is an important nuance in the literature. Westfall and Troendle argue that subset pivotality is needed primarily for that computational simplification, not for the closure principle itself. If every intersection hypothesis is supplied with a valid level-alpha resampling test, closed testing gives strong FWER control without requiring the same complete-null distribution for every subset. [[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/)] This is less restrictive in principle but can be exponentially more computationally demanding.

A permutation test can be finite-sample exact when the null makes the data invariant under the permitted transformation group, or when inference is conditioned on the actual randomized assignment mechanism. Randomly sampling only some permutations does not necessarily make the test merely asymptotic: properly constructed Monte Carlo permutation tests can retain exact size. Phipson and Smyth emphasize that the common estimator "number of sampled permutations more extreme divided by \(B\)" is not the right exact Monte Carlo p-value and can even produce impossible zero p-values. [[18](https://pubmed.ncbi.nlm.nih.gov/21044043/)]

Bootstrap tests do not have that same generic finite-randomization exactness. DiCiccio and Efron describe bootstrap CIs as approximate confidence intervals whose attraction is improved asymptotic coverage accuracy. FDA likewise cautions that bootstrap resampling may require large sample sizes and validation of Type I error behavior. [[9](https://www.fda.gov/media/162416/download?action=ReferProject), [12](https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf)]

### How a bootstrap interval relates to an adjusted p-value

There are three conceptually different constructions:

1. **Adjusted confidence level.** For Bonferroni, construct each marginal bootstrap CI at \(1-\alpha/m\), or at a prespecified weighted level. This is simple simultaneous control, assuming each bootstrap CI has the claimed marginal coverage approximation.

2. **Simultaneous bootstrap quantile.** Resample all outcomes jointly, calculate a studentized vector or max statistic on each resample, estimate a common critical value from the max distribution, and use that critical value in every interval. This is the interval analogue of a single-step maxT test.

3. **Test inversion.** For every candidate value \(\theta_i=\theta_0\), construct the multiplicity-adjusted bootstrap or permutation test and include \(\theta_0\) in the interval exactly when that null is not rejected. This gives the strongest test-CI compatibility, but is computationally more expensive.

A numerical adjusted p-value by itself is insufficient to recover the interval. The value \(p_{\mathrm{adj}}(0)\) tells you the adjusted test evidence at the null value zero. It does not tell you \(p_{\mathrm{adj}}(\theta_0)\) for every other candidate \(\theta_0\), which is what test inversion requires.

Percentile bootstrap intervals are different again. They use empirical quantiles of the bootstrap estimator distribution. DiCiccio and Efron survey percentile-related, BCa, bootstrap-t, ABC, and calibrated procedures and emphasize that bootstrap CIs are approximate inferential procedures in their own right. [[12](https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf)] A percentile CI therefore should not be made "multiplicity compatible" by setting its nominal confidence to \(1-p_{\mathrm{adj}}\). That operation has no general inferential justification.

### Extreme percentile levels and Monte Carlo error

Suppose a bound uses the 99.5th percentile of a bootstrap distribution. Its upper-tail probability is \(p=0.005\). With \(B\) independent bootstrap replicates, the number falling beyond the target population quantile has approximately a binomial distribution with mean \(Bp\). Jonsson and Nyberg use the corresponding binomial standard error

\[
SE(\hat p)=\sqrt{\frac{p(1-p)}{B}}
\]

to quantify the numerical precision of estimated percentile limits. [[19](https://pmc.ncbi.nlm.nih.gov/articles/PMC9197539/), [20](https://pubmed.ncbi.nlm.nih.gov/35353958/)]

For \(p=0.005\), this gives:

| Bootstrap replicates \(B\) | Expected replicates in 0.5% tail | SE of estimated tail probability | Relative SE \(SE/p\) | Approximate 95% half-width on tail probability |
|---:|---:|---:|---:|---:|
| 1,000 | 5 | 0.00223 | 44.6% | 0.00437 |
| 5,000 | 25 | 0.000998 | 20.0% | 0.00196 |
| 10,000 | 50 | 0.000705 | 14.1% | 0.00138 |
| 20,000 | 100 | 0.000499 | 10.0% | 0.000978 |
| 50,000 | 250 | 0.000315 | 6.31% | 0.000618 |
| 100,000 | 500 | 0.000223 | 4.46% | 0.000437 |
| 500,000 | 2,500 | 0.0000997 | 1.99% | 0.000195 |

These values are direct calculations from the binomial formula above. They describe Monte Carlo uncertainty in the **tail probability or rank**, not directly the uncertainty of the CI endpoint in percentage points, odds-ratio units, or another effect scale. Endpoint-value error also depends on how densely the bootstrap distribution is concentrated near the 99.5th percentile. A flat or sparse tail can produce substantial endpoint movement even when the tail-probability error looks small.

A convenient inversion of the same formula is

\[
B \approx \frac{1-p}{p r^2}
\]

for target relative standard error \(r\). At \(p=.005\), roughly 19,900 replicates give a 10% relative SE in the tail probability, about 79,600 give 5%, and about 497,500 give 2%. These are precision calculations, not universal minimum sample-size prescriptions.

Multiplicity can require tails much more extreme than 0.5%. For example, a two-sided Bonferroni interval with \(m=50\) hypotheses and familywise \(\alpha=.05\) allocates only \(.0005\) to each tail if equal allocation is used. Then even 100,000 replicates produce only about 50 observations in the relevant tail. The computational burden therefore grows quickly with large multiplicity families.

### Binary counts and overlapping metrics

For a single binary endpoint in a two-arm test, aggregate successes and denominators can be enough for many standard tests and CIs. For a **joint multiplicity resampling procedure**, however, the dependence structure across metrics matters.

If the same visit or send contributes to many metrics, a valid joint resampling implementation would ordinarily resample the unit as a vector:

\[
(Y_{i1},Y_{i2},\ldots,Y_{im},S_i,Z_i),
\]

where \(Z_i\) is treatment assignment and \(S_i\) contains segment membership. All metric statistics, all variation-control contrasts, and all applicable segment-cell statistics are then recomputed from the same resampling draw. This preserves the empirically observed dependence among metrics and the dependence caused by a shared control.

If instead the input consists solely of separate records such as "control: 1,200 successes of 10,000; treatment: 1,260 of 10,100" for every metric and segment, the cross-metric joint outcomes of the same units are not recoverable from those margins. Many different joint data sets have identical marginal counts but different correlations and therefore different maxT reference distributions. Thus a dependence-exploiting Westfall-Young procedure cannot generally be reconstructed exactly from marginal count tables alone.

For permutation inference, the randomization unit is critical. If users are randomized but multiple visits or sends from one user appear in the counts, permuting visits independently would not reproduce the original assignment mechanism. The randomization or cluster unit must be preserved. This follows from randomization-test logic rather than from the binary nature of the outcome itself.

### Reproducibility and auditability

Exact enumeration of the relevant assignment or permutation space is deterministic and avoids simulation error when feasible. Random bootstrap or Monte Carlo permutation implementations need more records for reproducibility.

For a reproducible statistical audit, the computational record should identify at least the input data version or hash, randomization unit, resampling unit, any strata or blocking, test statistic, effect scale, null-enforcing transformation, multiplicity family, maxT/minP or step-down algorithm, \(B\), treatment of ties, Monte Carlo p-value formula, RNG algorithm, seed or stream state, software and package versions, source-code version, and treatment of failed or degenerate resamples. FDA's computerized-systems guidance emphasizes retaining enough source data, audit trail, and computational logic to reconstruct reported results, while ICH E9 states that statistical software and its validation should be documented. [[6](https://database.ich.org/sites/default/files/E9_Guideline.pdf), [21](https://www.fda.gov/inspections-compliance-enforcement-and-criminal-investigations/fda-bioresearch-monitoring-information/guidance-industry-computerized-systems-used-clinical-trials)]

A fixed seed makes a pseudorandom analysis repeatable only together with the RNG algorithm, software version, and parallel-stream convention. Parallel execution can otherwise change the mapping from random numbers to resamples even under the same nominal seed. Exact enumeration, when practical, is the most straightforward deterministic alternative.

For random permutation p-values, the implementation should not report zero merely because no sampled permutation exceeded the observed statistic. Phipson and Smyth show that correct Monte Carlo treatment produces a valid discrete p-value for any number of random permutations. [[18](https://pubmed.ncbi.nlm.nih.gov/21044043/)] This becomes important in multiplicity settings because tiny p-values are precisely where a one-count Monte Carlo discrepancy can affect adjusted decisions.

### Evidence for the key claims

| Key claim | Source and location | Verbatim passage seen |
|---|---|---|
| Resampling can exploit endpoint dependence while controlling overall Type I error. | FDA 2022, Appendix Section 6, p. 21, lines 252-259. [[9](https://www.fda.gov/media/162416/download?action=ReferProject)] | "one general statistical approach that can provide more power" |
| Subset pivotality is tied to the maxT computational shortcut. | Westfall and Troendle 2008, Section 2.3, lines 94-99. [[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/)] | "assume a model that implies 'subset pivotality'" |
| Full resampling closed testing can have strong FWER control with much weaker assumptions. | Westfall and Troendle 2008, Summary and Introduction, lines 71-83. [[11](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/)] | "no assumptions whatsoever on the data-generating process are required" |
| Bootstrap CIs are approximate rather than generically exact. | DiCiccio and Efron 1996, Abstract and Section 1, p. 189. [[12](https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf)] | "bootstrap methods for producing good approximate confidence intervals" |
| Random-permutation p-values need exact Monte Carlo treatment. | Phipson and Smyth 2010, abstract, lines 239-242. [[18](https://pubmed.ncbi.nlm.nih.gov/21044043/)] | "permutation should be viewed as generating an exact discrete null distribution" |
| Percentile precision depends on the number of bootstrap samples. | Jonsson and Nyberg 2022, abstract and indexed methods extract. [[19](https://pmc.ncbi.nlm.nih.gov/articles/PMC9197539/), [20](https://pubmed.ncbi.nlm.nih.gov/35353958/)] | "trade-off between the number of bootstrap samples... and the increase in precision" |
| FDA requires reconstruction and retention of computational methods in regulated systems. | FDA computerized-systems guidance, audit/data-retention sections. [[21](https://www.fda.gov/inspections-compliance-enforcement-and-criminal-investigations/fda-bioresearch-monitoring-information/guidance-industry-computerized-systems-used-clinical-trials)] | "source documents should be retained to enable a reconstruction and evaluation of the trial" |

## Source inventory, access status, and limitations

The table distinguishes what I actually inspected from material that was available only through abstracts, indexed extracts, or metadata. "Full text" means I could inspect the body of the article or guidance. For PDFs used analytically, I also rendered pages to verify the visual text where relevant.

| Full reference | Working link | What I actually read |
|---|---|---|
| Benjamini Y, Yekutieli D. *False Discovery Rate-Adjusted Multiple Confidence Intervals for Selected Parameters.* Journal of the American Statistical Association. 2005;100(469):71-81. DOI 10.1198/016214504000001907. | [Author-hosted PDF](https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf) | **Full text PDF**, including Definitions 1 and 3, Remark 1, Theorems 1, 3 and 4, and the two-sided BH discussion. |
| Strassburger K, Bretz F. *Compatible simultaneous lower confidence bounds for the Holm procedure and other Bonferroni-based closed tests.* Statistics in Medicine. 2008;27(24):4914-4927. DOI 10.1002/sim.3338. | [PubMed](https://pubmed.ncbi.nlm.nih.gov/18618415/) | **Abstract only.** Publisher full article not opened. |
| Guilbaud O. *Simultaneous Confidence Regions Corresponding to Holm's Step-Down Procedure and Other Closed-Testing Procedures.* Biometrical Journal. 2008;50(5):678-692. DOI 10.1002/bimj.200710449. | [PubMed](https://pubmed.ncbi.nlm.nih.gov/18932131/) | **Abstract/search extract only.** Publisher full article not opened. |
| Guilbaud O. *Simultaneous confidence regions for closed tests, including Holm-, Hochberg-, and Hommel-related procedures.* Biometrical Journal. 2012;54(3):317-342. DOI 10.1002/bimj.201100123. | [PubMed](https://pubmed.ncbi.nlm.nih.gov/22685000/) | **Full PubMed abstract only.** Wiley full text not opened. |
| Hayter AJ, Hsu JC. *On the Relationship between Stepwise Decision Procedures and Confidence Sets.* Journal of the American Statistical Association. 1994;89(425):128-136. DOI 10.1080/01621459.1994.10476453. | [Taylor & Francis abstract](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476453) | **Abstract/search extract only.** Direct publisher open returned 403. |
| Brannath W, Kluge L, Scharpenberg M. *Informative simultaneous confidence intervals for graphical test procedures.* Statistical Methods in Medical Research. Published 2025;35(1):101-117. DOI 10.1177/09622802251393666. | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/) | **Full text.** |
| Westfall PH, Young SS. *Resampling-Based Multiple Testing: Examples and Methods for p-Value Adjustment.* Wiley; 1993. | [Wiley search/catalog record](https://www.wiley.com/) | **Publisher metadata/description and indexed references only. Full book not opened.** |
| Westfall PH, Troendle JF. *Multiple Testing with Minimal Assumptions.* Biometrical Journal. 2008;50(5):745-755. DOI 10.1002/bimj.200710456. | [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/) | **Full text.** |
| DiCiccio TJ, Efron B. *Bootstrap Confidence Intervals.* Statistical Science. 1996;11(3):189-228. DOI 10.1214/ss/1032280214. | [Accessible PDF](https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf) | **Full text PDF**, with page rendering checked. |
| Phipson B, Smyth GK. *Permutation P-values Should Never Be Zero: Calculating Exact P-values When Permutations Are Randomly Drawn.* Statistical Applications in Genetics and Molecular Biology. 2010;9:Article 39. DOI 10.2202/1544-6115.1585. | [PubMed](https://pubmed.ncbi.nlm.nih.gov/21044043/) | **Full abstract only.** |
| Jonsson EN, Nyberg J. *A quantitative approach to the choice of number of samples for percentile estimation in bootstrap and visual predictive check analyses.* CPT: Pharmacometrics & Systems Pharmacology. 2022;11(6):673-686. DOI 10.1002/psp4.12790. | [PubMed](https://pubmed.ncbi.nlm.nih.gov/35353958/) | **Full abstract plus indexed full-text extracts.** PMC direct open hit a browser challenge. |
| U.S. Food and Drug Administration. *Multiple Endpoints in Clinical Trials: Guidance for Industry.* October 2022. | [FDA PDF](https://www.fda.gov/media/162416/download) | **Full guidance PDF**, including resampling section, with page rendering checked. |
| International Council for Harmonisation. *E9 Statistical Principles for Clinical Trials.* Step 4, 5 February 1998. | [ICH PDF](https://database.ich.org/sites/default/files/E9_Guideline.pdf) | **Full guidance PDF**, including Sections 5.5 and 5.6, with pages rendered. |
| European Medicines Agency. *Multiplicity issues in clinical trials - Scientific guideline.* Current page for CPMP/EWP/908/99. | [EMA guideline page](https://www.ema.europa.eu/en/scientific-guidelines/multiplicity-issues-clinical-trials) | **Full web page.** The linked current PDF returned HTTP 429 in this session. |
| European Medicines Agency. *Draft guideline on multiplicity issues in clinical trials.* EMA/CHMP/44762/2017. | [EMA draft PDF](https://www.ema.europa.eu/en/documents/scientific-guideline/draft-guideline-multiplicity-issues-clinical-trials_en.pdf) | **Full PDF text inspected in the research session.** It is a draft, not the current effective version. |
| New England Journal of Medicine. *New Manuscripts: Statistical Guidelines.* Current Author Center guidance. | [NEJM Author Center](https://www.nejm.org/author-center/new-manuscripts) | **Indexed current-page extract read; direct open returned HTTP 403.** |
| JAMA. *Instructions for Authors.* Current statistical reporting instructions. | [JAMA Instructions for Authors](https://jamanetwork.com/journals/jama/pages/instructions-for-authors) | **Full web page**, including multiple comparisons and data-presentation sections. |
| Butcher NJ, Monsour A, Mew EJ, et al. *Guidelines for Reporting Outcomes in Trial Reports: The CONSORT-Outcomes 2022 Extension.* JAMA. 2022;328(22):2252-2264. DOI 10.1001/jama.2022.21022. | [JAMA article](https://jamanetwork.com/journals/jama/fullarticle/2799401) | **Web/indexed full-text extracts**, including checklist item 12a.1. |
| FDA. *Guidance for Industry: Computerized Systems Used in Clinical Trials.* | [FDA](https://www.fda.gov/inspections-compliance-enforcement-and-criminal-investigations/fda-bioresearch-monitoring-information/guidance-industry-computerized-systems-used-clinical-trials) | **Full web guidance sections on audit trails, reconstruction, and system documentation.** |

The sources I could **not open in full** were Strassburger and Bretz 2008, Guilbaud 2008, Guilbaud 2012, Hayter and Hsu 1994, Westfall and Young 1993, Phipson and Smyth 2010, and the current EMA CPMP/EWP/908/99 PDF. The current NEJM author-center page also blocked direct opening, although its indexed current-page text was available. Jonsson and Nyberg's PMC full page presented a browser challenge, so I relied on PubMed plus indexed article extracts. I did not bypass any paywall, login, 403 response, 429 response, or browser challenge.

A further limitation is that some of the strongest published compatible-interval results are presented in the accessible abstracts as one-sided. The report therefore distinguishes verified one-sided constructions from the more general statement that two-sided confidence sets can be obtained by inversion of two-sided tests. The strongest direct two-sided stepwise source I could verify from its abstract was Hayter and Hsu's two-parameter example. [[14](https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476453)]

Finally, regulatory and journal sources are used here as reporting precedents, not as claims that a commercial A/B-testing system is legally governed by pharmaceutical trial regulation. The relevant common issue is statistical presentation of multiple frequentist inferences, not regulatory jurisdiction.

## Linked citation sources

1. https://pmc.ncbi.nlm.nih.gov/articles/PMC12783375/
2. https://pubmed.ncbi.nlm.nih.gov/18618415/
3. https://pubmed.ncbi.nlm.nih.gov/22685000/
4. https://www.math.tau.ac.il/~yekutiel/papers/JASA%20FCR%20prints.pdf
5. https://pmc.ncbi.nlm.nih.gov/articles/PMC4564263
6. https://database.ich.org/sites/default/files/E9_Guideline.pdf
7. https://jamanetwork.com/journals/jama/pages/instructions-for-authors
8. https://www.ema.europa.eu/en/scientific-guidelines/multiplicity-issues-clinical-trials
9. https://www.fda.gov/media/162416/download?action=ReferProject
10. https://www.nejm.org/author-center/new-manuscripts
11. https://pmc.ncbi.nlm.nih.gov/articles/PMC3117234/
12. https://www.math.wustl.edu/~kuffner/AlastairYoung/LeeYoung1996b.pdf
13. https://pubmed.ncbi.nlm.nih.gov/18932131/
14. https://www.tandfonline.com/doi/abs/10.1080/01621459.1994.10476453
15. https://pmc.ncbi.nlm.nih.gov/articles/PMC4806862/
16. https://journal.r-project.org/articles/RJ-2010-008/
17. https://jamanetwork.com/journals/jama/fullarticle/2799401
18. https://pubmed.ncbi.nlm.nih.gov/21044043/
19. https://pmc.ncbi.nlm.nih.gov/articles/PMC9197539/
20. https://pubmed.ncbi.nlm.nih.gov/35353958/
21. https://www.fda.gov/inspections-compliance-enforcement-and-criminal-investigations/fda-bioresearch-monitoring-information/guidance-industry-computerized-systems-used-clinical-trials