# matstats
Assessing material test data for A or B basis.
## Disclaimer
**This software is provided "as is", without warranty of any kind**, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. 

1. **No Liability**: In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in an action of contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.
2. **Engineering Judgment**: This tool is intended for use by qualified engineers as a supplement to, not a replacement for, professional engineering judgment. The user is solely responsible for verifying the accuracy of the output and its compliance with relevant certification requirements (e.g., FAA, EASA, CMH-17).
3. **Statistical Validity**: The results provided by this program are dependent on the quality and quantity of the input data. Use of this program does not guarantee regulatory approval.

# Composite Material Basis Value Analysis Tool
### According to CMH-17-1G Statistics

This tool automates the calculation of **A-Basis** (99% reliability) and **B-Basis** (90% reliability) values for composite material data. It follows the official hierarchy for distribution selection.

---

## 1. Installation & Dependencies
Ensure you have Python 3 installed along with the following libraries:
`pip install pandas numpy scipy`

---

## 2. File Architecture
The program requires three specific input files to be present in the same folder:

| File | Description |
| :--- | :--- |
| **Input.csv** | Contains `DATA VALUES` and `DATA SET NO.` (for ANOVA). |
| **InputHK.csv** | Factors ($r$, $k_B$, $k_A$) for Hanson-Koopmans Non-Parametric logic. |
| **InputWeibull.csv** | Factors ($V_B$, $V_A$) for Weibull MLE Basis calculation. |

---

## 3. Statistical Hierarchy
The program does not simply pick the lowest number. It follows the **CMH-17 Selection Logic**:

1. **ANOVA (Pooled)**: Checked first. If batches are significantly different ($p \le 0.05$), the Pooled ANOVA method is prioritized.
2. **Normal**: If data is poolable, it is tested for Normality via Anderson-Darling. If it passes, this is the selected method.
3. **Lognormal**: Tested if Normal fails.
4. **Weibull**: Tested if Lognormal fails. Valid only if Shape ($\beta$) > 2.0.
5. **Non-Parametric (HK)**: The default method if no parametric distribution fits.

---

## 4. Key Calculation Methods

### ANOVA Pooled Method
Uses the Mean Square Between (MSB) and Mean Square Error (MSE) to calculate a Pooled Standard Deviation ($S_{pool}$):

$$S_{pool} = \sqrt{\max\left(0, \frac{MSB - MSE}{\bar{n}} + MSE\right)}$$

### Weibull (MLE) Method
Calculates basis values using the Scale ($\alpha$) and Shape ($\beta$) parameters combined with CMH-17 tolerance factors ($V$):
$$Result = \alpha \cdot [-\ln(R)]^{1/\beta} \cdot \exp\left(\frac{-V}{\beta \cdot \sqrt{n}}\right)$$

---

## 5. Interpreting the Output
The results are saved to **`final_basis_results.csv`**. 

* **Validity Column**: States "Pass" or "Fail" for each distribution test.
* **Selected Method**: Found at the top of the file; this is the statistically correct method to use for design according to the handbook.



---

## 6. Author Note
This program was developed to match the specific quantile and correction term requirements of CMH-17, Section 8.3.4.2.2.
