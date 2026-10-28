# The Rubin Causal Model: potential outcomes, randomization, Wald and LATE

This follows the lecture "Randomized Experiment, Wald Estimate, and Local Average Treatment Effect" and uses its notation. The running example is Miguel and Kremer (2004): deworming pills in Kenyan schools and their effect on school attendance.

**Notation (from the slides):**

| Symbol | Meaning |
|---|---|
| $D_i \in \{0,1\}$ | treatment actually received (took the deworming pill) |
| $Y_i(1),\ Y_i(0)$ | potential outcomes: $i$'s attendance *if* treated / *if* untreated |
| $Y_i$ | the observed outcome |
| $Z_i \in \{0,1\}$ | random assignment (school assigned to the treatment group) |
| $D_i(Z_i{=}1),\ D_i(Z_i{=}0)$ | potential treatment status under each assignment (the slides also write $D_{1i}, D_{0i}$) |
| $\pi_c, \pi_n, \pi_a$ | population shares of compliers, never-takers, always-takers |

Where a regression appears I use your notation $y = \mathbf{X}\beta + e$ with hats for estimates.

---

## 1. The core idea: two potential outcomes per person

The Rubin Causal Model (RCM) defines a causal effect as a comparison of two states of the world for the **same** unit. Every child $i$ has two potential outcomes:

$$
Y_i = \begin{cases} Y_i(1) & \text{if } D_i = 1 \\ Y_i(0) & \text{if } D_i = 0 \end{cases}
\qquad\Longleftrightarrow\qquad
Y_i = Y_i(0) + D_i\,[\,Y_i(1) - Y_i(0)\,].
$$

The causal effect of deworming for child $i$ is $Y_i(1) - Y_i(0)$.

**The fundamental problem of causal inference:** we only ever see one of the two. A child who took the pill shows us $Y_i(1)$, and $Y_i(0)$ is a counterfactual that never happened. So individual effects are never observed. That is why the slides move straight to *averages*, which we can sometimes learn.

> The model also quietly assumes that $i$'s outcome depends only on $i$'s own treatment, not on whether classmates were treated (this is called SUTVA). The slides don't name it, but keep it in mind for section 8, because deworming is a case where it fails.

## 2. What we want to estimate

**Average treatment effect (ATE):** the average effect over everyone.
$$
\text{ATE} = E[\,Y_i(1) - Y_i(0)\,]
$$

**Average treatment effect on the treated (ATET):** the average effect among those who actually got treated.
$$
\text{ATET} = E[\,Y_i(1) - Y_i(0) \mid D_i = 1\,] = E[Y_i(1)\mid D_i=1] - E[Y_i(0)\mid D_i=1]
$$

Neither is directly observable. In the ATET, $E[Y_i(0)\mid D_i=1]$ is "how often would the treated kids have attended if they had *not* been treated", which is a counterfactual.

ATE and ATET differ when effects are heterogeneous and the people who get treated are not a random slice of the population.

## 3. The naive comparison and selection bias

The obvious thing to compute is the difference in mean outcomes between treated and untreated. Both pieces are observed:
$$
E[Y_i \mid D_i = 1] - E[Y_i \mid D_i = 0] = E[Y_i(1)\mid D_i=1] - E[Y_i(0)\mid D_i=0].
$$

Add and subtract the counterfactual $E[Y_i(0)\mid D_i=1]$:

$$
\underbrace{E[Y_i(1)\mid D_i=1] - E[Y_i(0)\mid D_i=0]}_{\text{naive difference}}
= \underbrace{E[Y_i(1) - Y_i(0)\mid D_i=1]}_{\text{ATET}}
+ \underbrace{E[Y_i(0)\mid D_i=1] - E[Y_i(0)\mid D_i=0]}_{\text{selection bias}}
$$

Selection bias is the difference in *untreated* outcomes between the two groups. It is about who they are, not about the treatment. In the deworming setting, pills could only be given to children who were at school that day. Children who attend more are healthier or richer to begin with, so $E[Y_i(0)\mid D_i=1] > E[Y_i(0)\mid D_i=0]$, and the naive comparison **overstates** the effect.

**Regression view.** The naive difference is exactly the OLS slope in
$$
y_i = \beta_1 + \beta_2 D_i + e_i,\qquad \hat\beta_2 = \bar y_{D=1} - \bar y_{D=0}.
$$
Selection bias is the regression problem you already know as $E[e_i \mid D_i] \ne 0$: $D_i$ is correlated with the error, which contains $Y_i(0)$.

## 4. Random assignment removes selection bias

If $D_i$ is randomly assigned, it is **independent of the potential outcomes**:
$$
\{Y_i(1), Y_i(0)\} \perp D_i .
$$
Then $E[Y_i(0)\mid D_i=1] = E[Y_i(0)\mid D_i=0] = E[Y_i(0)]$, so selection bias is zero. Also $E[Y_i(1)-Y_i(0)\mid D_i=1] = E[Y_i(1)-Y_i(0)]$, so ATET = ATE, and

$$
E[Y_i\mid D_i=1] - E[Y_i\mid D_i=0] = E[Y_i(1) - Y_i(0)] = \text{ATE}.
$$

In regression terms, randomization makes $D_i$ uncorrelated with $e_i$, so $\hat\beta_2$ is unbiased for the ATE.

**Balance checks (Table I).** Randomization makes the groups equal *in expectation*. Miguel and Kremer check this by comparing pre-treatment (1998 baseline) characteristics across treatment and comparison schools. If the groups look similar on things measured before the program, that supports the claim that they are also similar on $Y_i(0)$, which we cannot see. In R this is just a series of `lm(x ~ z)` regressions (or t-tests) on baseline covariates.

## 5. The catch: assignment is not treatment


In the deworming study, the thing randomized was the **school** ($Z_i$), not who took the pill ($D_i$). Table III shows the gap. In 1998, among Group 1 (treatment) pupils who were girls under 13 or boys, only **0.78** received any medical treatment, and for girls 13 and older only **0.19** did (older girls were not given one of the drugs for safety reasons). Comparison schools had 0. In 1999, Group 3 (still comparison) had about **0.01** treated.

If everyone followed their assignment ($D_i = Z_i$), the slide shows we would be done:
$$
E[Y_i(1)\mid D_i=1] - E[Y_i(0)\mid D_i=0] = E[Y_i(1)\mid Z_i=1] - E[Y_i(0)\mid Z_i=0] = E[Y_i(1)-Y_i(0)].
$$
With imperfect take-up, the slides ask which comparison is right. There are three candidates:

| Comparison | What it compares | Problem |
|---|---|---|
| **As-treated** | took the pill vs did not | This is the naive comparison again, with selection bias. Who takes the pill is not random. |
| **Per-protocol** | took the pill (in treatment schools) vs everyone in comparison schools | Still biased: the pill-takers inside treatment schools are a selected group (they showed up). |
| **Intention to treat (ITT)** | everyone *assigned* to treatment vs everyone assigned to comparison | Unbiased, because $Z_i$ is random. But it measures the effect of *being offered* the pill, diluted by non-takers. |

The ITT is $E[Y_i\mid Z_i=1] - E[Y_i\mid Z_i=0]$. The rest of the lecture shows how to rescale it into a treatment effect.

## 6. Compliance types

Each child has two **potential treatment statuses**, $D_i(Z_i{=}1)$ and $D_i(Z_i{=}0)$: whether they would take the pill if their school were assigned to treatment, and if it were not. That splits the population into four types:

| Type | $D_i(Z_i{=}1)$ | $D_i(Z_i{=}0)$ | Deworming story |
|---|---|---|---|
| **Compliers** | 1 | 0 | take the pill only when offered |
| **Never-takers** | 0 | 0 | never treated (e.g. absent on pill day) |
| **Always-takers** | 1 | 1 | treated regardless (e.g. got deworming elsewhere) |
| **Defiers** | 0 | 1 | do the opposite of assignment |

Like potential outcomes, the type is never directly observed for an individual. A treated child in a treatment school could be a complier or an always-taker.

## 7. The Wald estimate and the LATE theorem

### Assumptions (A1 to A4 in the slides)

1. **Independence.** $Z_i$ is independent of the potential outcomes and potential treatments: $\{Y_i(\cdot), D_i(1), D_i(0)\} \perp Z_i$. Randomization delivers this.
2. **Exclusion.** $Z_i$ affects $Y$ only through $D$. Being in a treatment school does not change attendance except via the pill.
3. **First stage.** $E[D_i\mid Z_i=1] - E[D_i\mid Z_i=0] > 0$. Assignment actually moves take-up. The share of people who take the treatment ($D_i = 1$) must be higher among those assigned to it ($Z_i = 1$) than among those not assigned ($Z_i = 0$). 
4. **Monotonicity (no defiers).** Nobody is pushed *out* of treatment by being assigned to it.

### Step 1: recover the type shares

Under independence and no defiers:

- In comparison schools, only always-takers are treated: $E[D_i\mid Z_i=0] = \pi_a$.
- In treatment schools, always-takers and compliers are treated: $E[D_i\mid Z_i=1] = \pi_a + \pi_c$.

So
$$
\pi_a = E[D_i\mid Z_i=0],\qquad
\pi_c = E[D_i\mid Z_i=1] - E[D_i\mid Z_i=0],\qquad
\pi_n = 1 - E[D_i\mid Z_i=1].
$$

### Step 2: decompose the ITT by type

Split $E[Y_i\mid Z_i=1]$ by type. Independence means the type shares are the same in both arms, and exclusion means each type's outcome depends only on its $D$:
$$
E[Y_i\mid Z_i=1] = E[Y_i(1)\mid \text{complier}]\,\pi_c + E[Y_i(0)\mid \text{never-taker}]\,\pi_n + E[Y_i(1)\mid \text{always-taker}]\,\pi_a
$$
$$
E[Y_i\mid Z_i=0] = E[Y_i(0)\mid \text{complier}]\,\pi_c + E[Y_i(0)\mid \text{never-taker}]\,\pi_n + E[Y_i(1)\mid \text{always-taker}]\,\pi_a
$$

Subtract. Never-takers and always-takers have the same $D$ in both arms, so their terms cancel. Only compliers change treatment status, so only they contribute:
$$
E[Y_i\mid Z_i=1] - E[Y_i\mid Z_i=0] = E[\,Y_i(1) - Y_i(0)\mid \text{complier}\,]\cdot \pi_c .
$$

### Step 3: divide by the first stage

$$
\boxed{\;\text{Wald estimate} = \beta_{IV} = \frac{E[Y_i\mid Z_i=1] - E[Y_i\mid Z_i=0]}{E[D_i\mid Z_i=1] - E[D_i\mid Z_i=0]} = E[\,Y_i(1) - Y_i(0)\mid D_{1i} > D_{0i}\,] = \text{LATE}\;}
$$

In words: **ITT divided by the share of compliers equals the average effect for compliers.** This is the **local average treatment effect** (Imbens and Angrist, 1994).

### Is it the ATE? No.

The data say nothing about the effect for never-takers (we never see them treated) or always-takers (we never see them untreated). The LATE is "local" to the compliers. If deworming helps chronically absent children the most, and those children are mostly never-takers, the LATE can be well below the ATE.

Two notes on the slides' wording:

- Slide 18 calls the result "ATET". It is the effect for **compliers**, which equals the ATET only when there are no always-takers. That happens to be the 1998 deworming case: comparison schools had take-up of 0 (Table III), so $\pi_a = 0$, every treated child is a complier, and LATE = ATET there. In general they differ.
- The sample Wald estimator is a ratio of sample means, so strictly it is *consistent* for the LATE rather than exactly unbiased in finite samples.

### Regression form

The slides note that the ITT is the OLS slope from regressing $Y$ on a constant and $Z$. In full:

- **Reduced form:** $y_i = \gamma_1 + \gamma_2 Z_i + u_i$, with $\hat\gamma_2 = \bar y_{Z=1} - \bar y_{Z=0}$ (the ITT).
- **First stage:** $D_i = \delta_1 + \delta_2 Z_i + v_i$, with $\hat\delta_2 = \bar D_{Z=1} - \bar D_{Z=0}$ (estimates $\pi_c$).
- **Wald / IV:** $\hat\beta_{IV} = \hat\gamma_2 / \hat\delta_2$.

With one binary instrument and no covariates, this equals the 2SLS estimate of $\beta_2$ in $y_i = \beta_1 + \beta_2 D_i + e_i$ using $Z_i$ as the instrument for $D_i$. In matrix form with $\mathbf{X} = [\mathbf{1}, D]$ and $\mathbf{Z} = [\mathbf{1}, Z]$: $\hat\beta_{IV} = (\mathbf{Z}'\mathbf{X})^{-1}\mathbf{Z}'y$.

## 8. A caveat specific to deworming: spillovers

Finding 2 in the slides is that deworming improved health and attendance among *untreated* children in treatment schools and in neighboring schools. Worms spread, so treating your classmates helps you. This breaks the "my outcome depends only on my treatment" assumption from section 1, and it also bends exclusion: being in a treatment school helps a never-taker even without the pill. It is the reason Miguel and Kremer randomized at the **school** level rather than the pupil level. Within-school comparisons of takers vs non-takers would understate the effect, because the non-takers were partly protected too.

---

## 9. R simulation: see every estimand side by side

The simulation builds a population where we know everyone's type and both potential outcomes, so we can compute the true ATE, ATET and LATE, then check what each estimator recovers. The outcome is an attendance rate.

Design (population values):

| Type | Share | $E[Y(0)]$ | Effect $Y(1)-Y(0)$ |
|---|---|---|---|
| Complier | 0.6 | 0.75 | 0.10 |
| Never-taker | 0.3 | 0.55 | 0.20 |
| Always-taker | 0.1 | 0.85 | 0.02 |

Never-takers have low baseline attendance (they miss pill day) and would benefit most. Always-takers are healthy and gain little.

```r
set.seed(2004)
n <- 200000

type <- sample(c("complier", "never", "always"), n, replace = TRUE,
               prob = c(0.6, 0.3, 0.1))

# potential outcomes
y0 <- c(complier = 0.75, never = 0.55, always = 0.85)[type] + rnorm(n, sd = 0.05)
te <- c(complier = 0.10, never = 0.20, always = 0.02)[type]
y1 <- y0 + te

# random assignment and potential treatment status
z  <- rbinom(n, 1, 0.5)
d1 <- as.integer(type %in% c("complier", "always"))  # D_i(Z = 1)
d0 <- as.integer(type == "always")                    # D_i(Z = 0)
d  <- ifelse(z == 1, d1, d0)

# observed outcome: we only ever see one potential outcome
y  <- ifelse(d == 1, y1, y0)
dat <- data.frame(y, d, z, type)

# ---- true parameters (only knowable in a simulation) ----
c(ATE  = mean(y1 - y0),
  ATET = mean((y1 - y0)[d == 1]),
  LATE = mean((y1 - y0)[type == "complier"]))
#  ATE ~ 0.122, ATET ~ 0.080, LATE ~ 0.100

# ---- estimators ----
as_treated   <- coef(lm(y ~ d, data = dat))["d"]            # naive: ~0.205
per_protocol <- mean(y[z == 1 & d == 1]) - mean(y[z == 0])   # ~0.151
itt          <- coef(lm(y ~ z, data = dat))["z"]             # reduced form: ~0.060
first_stage  <- coef(lm(d ~ z, data = dat))["z"]             # pi_c: ~0.60
wald         <- itt / first_stage                            # ~0.100 = LATE

round(c(as_treated = unname(as_treated), per_protocol = per_protocol,
        itt = unname(itt), first_stage = unname(first_stage),
        wald = unname(wald)), 3)

# type shares recovered from the data alone
c(pi_a = mean(d[z == 0]),
  pi_c = mean(d[z == 1]) - mean(d[z == 0]),
  pi_n = 1 - mean(d[z == 1]))
```

What to notice:

- **As-treated** (~0.205) is far too big, because of selection bias: takers have high $Y(0)$.
- **Per-protocol** (~0.151) is still biased for the same reason.
- **ITT** (~0.060) is unbiased for the effect of *assignment*, but diluted: $0.10 \times 0.6$.
- **Wald** (~0.100) recovers the **LATE** exactly, not the ATE (0.122) or the ATET (0.080).

### The same thing as 2SLS

```r
# install.packages("AER")
library(AER)
iv <- ivreg(y ~ d | z, data = dat)   # instrument d with z
coef(iv)["d"]                        # identical to `wald` above
summary(iv)                          # gives a proper IV standard error

# or with fixest:
# library(fixest); feols(y ~ 1 | d ~ z, data = dat)
```

Use the IV regression rather than dividing two `lm` coefficients when you need standard errors. In a real cluster-randomized study like deworming you would also cluster by school, e.g. `feols(y ~ 1 | d ~ z, data = dat, cluster = ~school)`.

### The Wald estimate in DuckDB

Since it is just four conditional means, it is a one-query job:

```r
library(duckdb)
con <- dbConnect(duckdb())
duckdb_register(con, "dat", dat)

dbGetQuery(con, "
  SELECT
    AVG(y) FILTER (WHERE z = 1) - AVG(y) FILTER (WHERE z = 0) AS itt,
    AVG(d) FILTER (WHERE z = 1) - AVG(d) FILTER (WHERE z = 0) AS first_stage,
    itt / first_stage                                         AS wald
  FROM dat
")
dbDisconnect(con, shutdown = TRUE)
```

### Try this

1. Set the always-taker share to 0 (one-sided non-compliance, like deworming in 1998). Check that the true ATET and the LATE become equal, and that Wald recovers both.
2. Make all three effects equal (say 0.1). Now ATE = ATET = LATE, and only the as-treated and per-protocol comparisons are wrong. Heterogeneity is what makes "local" matter.
3. Add some defiers (`D(1) = 0, D(0) = 1`) with a different effect. Watch the Wald estimate drift away from the complier effect. That is why monotonicity (A4) is needed.

---

## One-paragraph summary

The RCM defines a causal effect as $Y_i(1) - Y_i(0)$, which is never observed for any one person, so we target averages such as the ATE and ATET. Comparing treated with untreated gives ATET plus selection bias. Random assignment makes treatment independent of potential outcomes, which kills selection bias, so a simple difference in means (the OLS slope on $D$) estimates the ATE. When people do not follow their assignment, compare by *assignment* instead (ITT) and divide by the first-stage difference in take-up. Under independence, exclusion, a non-zero first stage, and no defiers, that Wald ratio equals the average effect for **compliers**, the LATE, and not the ATE.
