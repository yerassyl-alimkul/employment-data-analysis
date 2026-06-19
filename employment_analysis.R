# 0. SETUP: IMPORT AND PREPARE THE DATA

emp <- read.table("employment.txt",
                  header     = TRUE,
                  sep        = "\t",
                  dec        = ",",
                  na.strings = "NA")

emp$empl_f    <- factor(emp$empl,        levels = c(0, 1), labels = c("Unemployed", "Employed"))
emp$sex_f     <- factor(emp$sex,         levels = c(0, 1), labels = c("Male", "Female"))
emp$cit_f     <- factor(emp$citizenship, levels = c(0, 1), labels = c("Non-Italian", "Italian"))
emp$pension_f <- factor(emp$pension,     levels = c(0, 1), labels = c("No pension", "Pension"))

cat(" DATA OVERVIEW\n")
cat("Number of observations          :", nrow(emp), "\n")
cat("Employed subjects (with salary) :", sum(emp$empl == 1), "\n")
cat("Unemployed subjects (salary NA) :", sum(emp$empl == 0), "\n\n")
str(emp[, 1:6])

# 0b. HELPER FUNCTIONS BUILT FROM THE LECTURE FORMULAS

association_indices <- function(tab) {
  N    <- sum(tab)
  nhat <- outer(rowSums(tab), colSums(tab)) / N   # independence table
  chi2 <- sum((tab - nhat)^2 / nhat)              # chi-squared index
  Phi  <- sqrt(chi2 / N)                          # Phi coefficient
  V    <- Phi / sqrt(min(dim(tab)) - 1)           # Cramer's V
  cat("Independence table (expected frequencies under independence):\n")
  print(round(nhat, 1))
  cat("Chi-squared index =", round(chi2, 3),
      "| Phi =", round(Phi, 4),
      "| Cramer's V =", round(V, 4), "\n")
  invisible(list(chi2 = chi2, Phi = Phi, V = V))
}

describe <- function(x, label) {
  x  <- x[!is.na(x)]
  N  <- length(x)
  m  <- mean(x); s <- sd(x)
  skew <- mean(((x - m) / s)^3)
  cat(sprintf("%-28s n = %5d | mean = %8.2f | median = %8.2f | sd = %7.2f\n",
              label, N, m, median(x), s))
  cat(sprintf("%-28s range = [%g; %g] | IQR = %.2f | CV = %.3f | skewness = %.3f\n",
              "", min(x), max(x), IQR(x), s / m, skew))
}

least_squares <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  b1 <- cov(x, y) / var(x)
  b0 <- mean(y) - b1 * mean(x)
  r  <- cor(x, y)
  mpe <- sd(y) * sqrt(1 - r^2)
  cat("Regression line:  y =", round(b0, 2), "+", round(b1, 3), "* x\n")
  cat("Coefficient of determination r^2 =", signif(r^2, 3), "\n")
  cat("Mean prediction error =", round(mpe, 2), "\n")
  invisible(list(b0 = b0, b1 = b1, r2 = r^2, mpe = mpe))
}

normal_interval95 <- function(x, label) {
  x <- x[!is.na(x)]
  m <- mean(x); s <- sd(x)
  z <- qnorm(0.975)
  lo <- m - z * s; hi <- m + z * s
  cat(sprintf("%s: n = %d, mean = %.2f, sd = %.2f\n", label, length(x), m, s))
  cat(sprintf("Central 95%% interval under the normal model: [%.2f; %.2f]\n", lo, hi))
  cat(sprintf("Empirical check: %.1f%% of the observed values fall inside it\n",
              100 * mean(x >= lo & x <= hi)))
  invisible(c(lo, hi))
}

# 0c. DESCRIPTIVE OVERVIEW OF THE DATASET (Lectures 2-6)
cat(" DESCRIPTIVE STATISTICS\n")
describe(emp$age,    "Age (years)")
describe(emp$salary, "Salary (EUR, employed only)")

cat("\nFrequency distributions of the categorical variables:\n")
for (v in c("empl_f", "sex_f", "cit_f", "pension_f")) {
  tb <- table(emp[[v]])
  cat("\n", v, ":\n"); print(tb)
  cat("Relative frequencies (%):\n"); print(round(prop.table(tb) * 100, 1))
}

dir.create("figures", showWarnings = FALSE)
png("figures/00_descriptives.png", width = 1200, height = 900, res = 120)
par(mfrow = c(2, 2))
hist(emp$age, breaks = 30, col = "steelblue", border = "white",
     main = "Distribution of age", xlab = "Age (years)")
hist(emp$salary, breaks = 40, col = "darkseagreen", border = "white",
     main = "Distribution of salary (employed)", xlab = "Salary (EUR)")
barplot(table(emp$empl_f), col = c("indianred", "seagreen"),
        main = "Employment status", ylab = "Frequency")
barplot(table(emp$sex_f), col = c("lightskyblue", "lightpink"),
        main = "Sex", ylab = "Frequency")
par(mfrow = c(1, 1))
dev.off()

#  Q1. IS THERE A LINK BETWEEN GENDER AND OCCUPATION (EMPLOYMENT STATUS)?

cat(" Q1. Gender vs employment status (contingency table analysis)\n")

tab1 <- table(Gender = emp$sex_f, Employment = emp$empl_f)
cat("Observed contingency table:\n"); print(tab1)
cat("\nConditional relative frequency distributions of employment given gender (%):\n")
print(round(prop.table(tab1, margin = 1) * 100, 2))
cat("\n")
association_indices(tab1)

png("figures/q1_gender_employment.png", width = 900, height = 600, res = 120)
barplot(t(prop.table(tab1, 1)) * 100, beside = TRUE,
        col = c("indianred", "seagreen"), legend.text = colnames(tab1),
        ylab = "%", main = "Conditional distribution of employment by gender",
        ylim = c(0, 100))
dev.off()

#  Q2. IS THERE A CORRELATION BETWEEN AGE AND SALARY?

cat(" Q2. Correlation between age and salary (employed subjects)\n")
ok <- !is.na(emp$salary)
r_pearson  <- cor(emp$age[ok], emp$salary[ok])
r_spearman <- cor(emp$age[ok], emp$salary[ok], method = "spearman")
cat("Pearson correlation coefficient r :", round(r_pearson, 4), "\n")
cat("Spearman rank correlation         :", round(r_spearman, 4), "\n")
cat("Reminder of the scale: |r| close to 0 = no linear correlation,\n")
cat("|r| close to 1 = perfect linear correlation.\n")

png("figures/q2_age_salary_scatter.png", width = 900, height = 700, res = 120)
plot(emp$age[ok], emp$salary[ok], pch = 16, cex = 0.3, col = rgb(0, 0, 1, 0.10),
     xlab = "Age (years)", ylab = "Monthly salary (EUR)",
     main = "Scatterplot of age vs salary with the least squares line")
fit_tmp <- least_squares(emp$age, emp$salary)
abline(fit_tmp$b0, fit_tmp$b1, col = "red", lwd = 2)
dev.off()

#  Q3. CAN AGE BE USED TO PREDICT SALARY?

cat(" Q3. Least squares regression: salary on age\n")
fit3 <- least_squares(emp$age, emp$salary)
cat("\nInterpretation:\n")
cat(" - Slope b1 =", round(fit3$b1, 3),
    "EUR: average change in salary per additional year of age\n")
cat(" - Example: predicted salary at age 60 =",
    round(fit3$b0 + fit3$b1 * 60, 2), "EUR\n")
cat(" - r^2 =", signif(fit3$r2, 3),
    ": share of salary variability explained by the regression line\n")
cat(" - Mean prediction error =", round(fit3$mpe, 2),
    "EUR, barely smaller than sd(salary) =",
    round(sd(emp$salary, na.rm = TRUE), 2), "EUR\n")

#  Q4. DOES GENDER AFFECT SALARY?

cat(" Q4. Salary by gender (mean dependence)\n")

describe(emp$salary[emp$sex == 0], "Salary, men")
describe(emp$salary[emp$sex == 1], "Salary, women")
dm <- mean(emp$salary[emp$sex == 0], na.rm = TRUE) -
      mean(emp$salary[emp$sex == 1], na.rm = TRUE)
cat("Difference of conditional means (men - women):", round(dm, 2), "EUR\n")

png("figures/q4_salary_by_gender.png", width = 800, height = 600, res = 120)
boxplot(salary ~ sex_f, data = emp, col = c("lightskyblue", "lightpink"),
        xlab = "", ylab = "Monthly salary (EUR)", main = "Salary by gender")
dev.off()

#  Q5. IS CITIZENSHIP ASSOCIATED WITH EMPLOYMENT STATUS AMONG PEOPLE 65+?

cat(" Q5. Citizenship vs employment among people aged 65+\n")
emp65 <- subset(emp, age >= 65)
cat("Subjects aged 65+ :", nrow(emp65), "\n")

tab5 <- table(Citizenship = emp65$cit_f, Employment = emp65$empl_f)
cat("Observed contingency table:\n"); print(tab5)
cat("\nConditional distributions of employment given citizenship (%):\n")
print(round(prop.table(tab5, 1) * 100, 2))
cat("\n")
association_indices(tab5)


#  Q6. DO PENSION RECIPIENTS EARN A DIFFERENT AVERAGE SALARY THAN NON-RECIPIENTS?

cat(" Q6. Salary by pension receipt (mean dependence)\n")
describe(emp$salary[emp$pension == 0], "Salary, no pension")
describe(emp$salary[emp$pension == 1], "Salary, pension")
dp <- mean(emp$salary[emp$pension == 1], na.rm = TRUE) -
      mean(emp$salary[emp$pension == 0], na.rm = TRUE)
cat("Difference of conditional means (pension - no pension):", round(dp, 2), "EUR\n")

png("figures/q6_salary_by_pension.png", width = 800, height = 600, res = 120)
boxplot(salary ~ pension_f, data = emp, col = c("wheat", "darkseagreen"),
        ylab = "Monthly salary (EUR)", main = "Salary by pension receipt")
dev.off()

#  Q7. IS AGE ASSOCIATED WITH RECEIVING A PENSION?

cat(" Q7. Age and pension receipt\n")
cat("(a) Conditional distributions of age:\n")
describe(emp$age[emp$pension == 0], "Age, no pension")
describe(emp$age[emp$pension == 1], "Age, pension")

cat("\n(b) Pension receipt vs the age-65 threshold:\n")
tab7 <- table(Pension = emp$pension_f, "Age 65+" = emp$age >= 65)
print(tab7)
association_indices(tab7)
cat("Age range of pension recipients:",
    paste(range(emp$age[emp$pension == 1]), collapse = " - "), "\n")
cat("Age range of non-recipients    :",
    paste(range(emp$age[emp$pension == 0]), collapse = " - "), "\n")
cat("-> Cramer's V = 1: pension receipt is COMPLETELY determined by age\n")
cat("   (everyone aged 65+ receives a pension, nobody younger does).\n")

png("figures/q7_age_by_pension.png", width = 800, height = 600, res = 120)
boxplot(age ~ pension_f, data = emp, col = c("wheat", "darkseagreen"),
        ylab = "Age (years)", main = "Age by pension receipt")
dev.off()

#  Q8. A 95% INTERVAL FOR THE SALARY OF EMPLOYED PEOPLE AGED 65 AND OLDER

cat(" Q8. 95% normal interval: salary of employed people aged 65+\n")
sal65 <- emp$salary[emp$empl == 1 & emp$age >= 65]
normal_interval95(sal65, "Salary of employed subjects aged 65+")

png("figures/q8_salary65_normal.png", width = 900, height = 600, res = 120)
hist(sal65, breaks = 30, freq = FALSE, col = "darkseagreen", border = "white",
     main = "Salary of employed people 65+ with fitted normal curve",
     xlab = "Monthly salary (EUR)")
curve(dnorm(x, mean(sal65), sd(sal65)), add = TRUE, col = "red", lwd = 2)
abline(v = mean(sal65) + c(-1, 1) * qnorm(0.975) * sd(sal65),
       col = "blue", lwd = 2, lty = 2)
dev.off()

#  Q9. DO EMPLOYED AND UNEMPLOYED INDIVIDUALS DIFFER IN THEIR AVERAGE AGE?

cat(" Q9. Age by employment status (mean dependence)\n")
describe(emp$age[emp$empl == 1], "Age, employed")
describe(emp$age[emp$empl == 0], "Age, unemployed")
da <- mean(emp$age[emp$empl == 0]) - mean(emp$age[emp$empl == 1])
cat("Difference of conditional means (unemployed - employed):",
    round(da, 2), "years\n")

png("figures/q9_age_by_employment.png", width = 800, height = 600, res = 120)
boxplot(age ~ empl_f, data = emp, col = c("indianred", "seagreen"),
        ylab = "Age (years)", main = "Age by employment status")
dev.off()


#  Q10. DO MEN EARN MORE THAN WOMEN? AND ARE THERE DIFFERENCES WHEN IT COMES TO PENSION?

cat(" Q10. Salary by gender, overall and within pension groups\n")
cat("Overall conditional means of salary by gender:\n")
print(round(tapply(emp$salary, emp$sex_f, mean, na.rm = TRUE), 2))

cat("\nConditional means of salary by gender WITHIN pension groups:\n")
m10 <- tapply(emp$salary, list(emp$sex_f, emp$pension_f), mean, na.rm = TRUE)
print(round(m10, 2))
cat("\nGender gap (men - women) inside each pension group:\n")
print(round(m10["Male", ] - m10["Female", ], 2))

png("figures/q10_salary_sex_pension.png", width = 950, height = 600, res = 120)
boxplot(salary ~ sex_f + pension_f, data = emp,
        col = c("lightskyblue", "lightpink"),
        names = c("M / no pens.", "F / no pens.", "M / pension", "F / pension"),
        ylab = "Monthly salary (EUR)", main = "Salary by sex and pension")
dev.off()

#  Q11. IS CITIZENSHIP ASSOCIATED WITH PENSION RECEIPT?

cat(" Q11. Citizenship vs pension receipt (contingency table analysis)\n")
tab11 <- table(Citizenship = emp$cit_f, Pension = emp$pension_f)
cat("Observed contingency table:\n"); print(tab11)
cat("\nConditional distributions of pension given citizenship (%):\n")
print(round(prop.table(tab11, 1) * 100, 2))
cat("\n")
association_indices(tab11)

#  Q12. DO ITALIAN CITIZENS HAVE A DIFFERENT AVERAGE SALARY THAN NON-ITALIANS?

cat(" Q12. Salary by citizenship (mean dependence)\n")
describe(emp$salary[emp$citizenship == 0], "Salary, non-Italians")
describe(emp$salary[emp$citizenship == 1], "Salary, Italians")
dc <- mean(emp$salary[emp$citizenship == 1], na.rm = TRUE) -
      mean(emp$salary[emp$citizenship == 0], na.rm = TRUE)
cat("Difference of conditional means (Italians - non-Italians):",
    round(dc, 2), "EUR\n")

png("figures/q12_salary_by_citizenship.png", width = 800, height = 600, res = 120)
boxplot(salary ~ cit_f, data = emp, col = c("khaki", "palegreen3"),
        ylab = "Monthly salary (EUR)", main = "Salary by citizenship")
dev.off()

#  Q13. DOES AGE AFFECT SALARY?

cat(" Q13. Mean dependence of salary on age classes\n")
emp$age_class <- cut(emp$age, breaks = c(30, 40, 50, 60, 70, 90),
                     labels = c("31-40", "41-50", "51-60", "61-70", "71+"),
                     include.lowest = TRUE)
cat("Conditional means of salary by age class:\n")
print(round(tapply(emp$salary, emp$age_class, mean, na.rm = TRUE), 2))
cat("\nConditional standard deviations:\n")
print(round(tapply(emp$salary, emp$age_class, sd, na.rm = TRUE), 2))
cat("\nNumber of employed subjects per class:\n")
print(tapply(!is.na(emp$salary), emp$age_class, sum))
cat("\nFor comparison, the overall mean salary is",
    round(mean(emp$salary, na.rm = TRUE), 2), "EUR and the regression of Q3\n")
cat("gives slope b1 =", round(fit3$b1, 3), "EUR/year with r^2 =",
    signif(fit3$r2, 3), "\n")

png("figures/q13_salary_age_class.png", width = 900, height = 600, res = 120)
boxplot(salary ~ age_class, data = emp, col = "lightsteelblue",
        xlab = "Age class", ylab = "Monthly salary (EUR)",
        main = "Salary by age class")
dev.off()

#  Q14. A 95% INTERVAL FOR THE AGE OF YOUNGER EMPLOYED PEOPLE

cat(" Q14. 95% normal interval: age of employed people younger than 65\n")
age_young <- emp$age[emp$empl == 1 & emp$age < 65]
normal_interval95(age_young, "Age of employed subjects younger than 65")

png("figures/q14_age_young_normal.png", width = 900, height = 600, res = 120)
hist(age_young, breaks = 30, freq = FALSE, col = "steelblue", border = "white",
     main = "Age of employed people < 65 with fitted normal curve",
     xlab = "Age (years)")
curve(dnorm(x, mean(age_young), sd(age_young)), add = TRUE, col = "red", lwd = 2)
abline(v = mean(age_young) + c(-1, 1) * qnorm(0.975) * sd(age_young),
       col = "blue", lwd = 2, lty = 2)
dev.off()

cat(" ANALYSIS COMPLETE - all figures saved in ./figures/\n")