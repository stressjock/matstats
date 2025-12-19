import pandas as pd
import numpy as np
from scipy import stats
import csv

"""
DISCLAIMER:
This program is provided for educational and engineering analysis purposes only. 
The author assumes NO RESPONSIBILITY for the accuracy, completeness, or 
reliability of the calculations provided. Users must verify all results against 
official CMH-17 tables and standards. Use of this software constitutes 
acceptance of all risks associated with its output.
"""

def get_normal_k(n, basis='B'):
    """Calculates normal tolerance factors based on MIL-HDBK-17 approximations."""
    if n > 15:
        if basis == 'B':
            return 1.282 + np.exp(0.958 - 0.52 * np.log(n) + 3.19 / n)
        else: # Basis A
            return 2.326 + np.exp(1.34 - 0.522 * np.log(n) + 3.87 / n)
    else:
        p, gamma = (0.90, 0.95) if basis == 'B' else (0.99, 0.95)
        delta = stats.norm.ppf(p) * np.sqrt(n)
        return stats.nct.ppf(gamma, df=n-1, nc=delta) / np.sqrt(n)

def get_hk_factors(n):
    """Retrieves H-K factors from InputHK.csv."""
    try:
        hk_df = pd.read_csv('InputHK.csv')
        a_row = hk_df[hk_df.iloc[:, 4] == n]
        ka = a_row.iloc[0, 5] if not a_row.empty else None
        b_row = hk_df[hk_df.iloc[:, 0] == n]
        rb = int(b_row.iloc[0, 1]) if not b_row.empty and not np.isnan(b_row.iloc[0, 1]) else None
        kb = b_row.iloc[0, 2] if not b_row.empty and not np.isnan(b_row.iloc[0, 2]) else None
        return rb, kb, ka
    except Exception:
        return None, None, None

def get_weibull_factors(n):
    """Retrieves Weibull V factors from InputWeibull.csv (n, V_b, V_a)."""
    try:
        w_df = pd.read_csv('InputWeibull.csv')
        row = w_df[w_df.iloc[:, 0] == n]
        if not row.empty:
            return row.iloc[0, 1], row.iloc[0, 2]
        return None, None
    except Exception:
        return None, None

def run_analysis(csv_file='Input.csv', value_col='DATA VALUES', set_col='DATA SET NO.'):
    try:
        df = pd.read_csv(csv_file)
        clean_df = df.dropna(subset=[value_col])
        data = np.sort(clean_df[value_col].values)
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    n = len(data)
    mean_val, sd_val = np.mean(data), np.std(data, ddof=1)
    cov = (sd_val / mean_val) * 100 if mean_val != 0 else 0
    
    # 1. ANOVA CMH-17 Pooled Method
    msb, mse, s_pool, anova_p = 0, 0, 0, None
    tb, ta = 1.93, 3.30 
    b_anova, a_anova = None, None
    if set_col in clean_df.columns and clean_df[set_col].nunique() > 1:
        groups = [group[value_col].values for name, group in clean_df.groupby(set_col)]
        k = len(groups)
        n_i = np.array([len(g) for g in groups])
        n_bar = (n - np.sum(n_i**2) / n) / (k - 1)
        msb = sum(len(g) * (np.mean(g) - mean_val)**2 for g in groups) / (k - 1)
        mse = sum(sum((x - np.mean(g))**2 for x in g) for g in groups) / (n - k)
        f_stat, anova_p = stats.f_oneway(*groups)
        s_pool = np.sqrt(max(0, (msb - mse) / n_bar + mse))
        b_anova = mean_val - (tb * s_pool)
        a_anova = mean_val - (ta * s_pool)

    # 2. Normal Analysis
    ad_norm = stats.anderson(data, dist='norm')
    pass_norm = ad_norm.statistic < ad_norm.critical_values[2]
    knb, kna = get_normal_k(n, 'B'), get_normal_k(n, 'A')
    b_norm, a_norm = mean_val - knb * sd_val, mean_val - kna * sd_val

    # 3. Lognormal Analysis
    log_data = np.log(data)
    ad_log = stats.anderson(log_data, dist='norm')
    pass_log = ad_log.statistic < ad_log.critical_values[2]
    b_log = np.exp(np.mean(log_data) - knb * np.std(log_data, ddof=1))
    a_log = np.exp(np.mean(log_data) - kna * np.std(log_data, ddof=1))

    # 4. Weibull Analysis (MLE and Basis)
    shape_hat, loc_hat, scale_hat = stats.weibull_min.fit(data, floc=0)
    vb, va = get_weibull_factors(n)
    if vb is not None and va is not None:
        q_b, q_a = (-np.log(0.90))**(1/shape_hat), (-np.log(0.99))**(1/shape_hat)
        corr_b = np.exp(-vb / (shape_hat * np.sqrt(n)))
        corr_a = np.exp(-va / (shape_hat * np.sqrt(n)))
        b_weibull = scale_hat * q_b * corr_b
        a_weibull = scale_hat * q_a * corr_a
    else:
        b_weibull, a_weibull = None, None
    pass_weibull = True if shape_hat > 2.0 else False

    # 5. Nonparametric (HK)
    rb, kb, ka = get_hk_factors(n)
    b_hk = data[rb-1] * (data[0] / data[rb-1])**kb if rb else None
    a_hk = data[n-1] * (data[0] / data[n-1])**ka if ka else None

    # SELECTION LOGIC
    if anova_p is not None and anova_p <= 0.05:
        selected_method = "ANOVA (Pooled)"
    elif pass_norm:
        selected_method = "Normal"
    elif pass_log:
        selected_method = "Lognormal"
    elif pass_weibull:
        selected_method = "Weibull"
    else:
        selected_method = "Non-Parametric (HK)"

    # --- CSV FORMATION ---
    row1 = ['N', 'Mean', 'SD', 'COV', 'SELECTED METHOD:', selected_method, '', '']
    row2 = [n, round(mean_val, 4), round(sd_val, 4), round(cov, 4), '', '', '', '']
    row3 = ['', '', '', '', '', 'B-Basis', '', 'A-Basis', 'Validity', '']
    
    # ANOVA Row
    row5 = ['ANOVA CMH-17 Analysis', '', '', '', '', '', '', '','']
    row6 = ['', 'MSB', 'MSE', 'S_pool', 'Tb', 'B_Val (ANOVA)', 'Ta', 'A_Val (ANOVA)', 'Pass (Poolable)' if anova_p and anova_p > 0.05 else 'Fail']
    row7 = ['', round(msb, 2), round(mse, 2), round(s_pool, 2), tb, round(b_anova, 2) if b_anova else '', ta, round(a_anova, 2) if a_anova else '', '']

    # Normal Row
    row9 = ['Normal Analysis', '', '', '', '', '', '', '','']
    row10 = ['', '','','', 'Normal_Kb', 'Normal_B_Val', 'Normal_Ka', 'Normal_A_Value', 'AD_Stat', 'Pass' if pass_norm else 'Fail']
    row11 = ['', '','','', round(knb, 4), round(b_norm, 4), round(kna, 4), round(a_norm, 4), round(ad_norm.statistic, 4), '']
    
    # Lognormal Row
    row13 = ['Lognormal Analysis', '', '', '', '', '', '', '','']
    row14 = ['', '','','', 'Log_Kb', 'Log_B_Val', 'Log_Ka', 'Log_A_Value', 'AD_Stat', 'Pass' if pass_log else 'Fail']
    row15 = ['', '','','', round(knb, 4), round(b_log, 4), round(kna, 4), round(a_log, 4), round(ad_log.statistic, 4), '']

    # Non-Parametric HK Row
    row17 = ['Non Parametric (HK) Analysis', '', '', '', '', '', '', '','']
    row18 = ['', '','HK_rb', 'HK_r', 'HK_kb', 'HK_B_Val', 'HK_ka', 'HK_A_Value', '', '']
    row19 = ['', '', rb, data[rb-1] if rb else '', round(kb, 4) if kb else '', round(b_hk, 4) if b_hk else '', ka, round(a_hk, 4) if a_hk else '', '', '']
    
    # Weibull Row
    row21 = ['Weibull Analysis', '', '', '', '', '', '', '']
    row22 = ['', '', 'Shape', 'Scale', 'W_VB', 'B_Val (Weibull)', 'W_VA', 'A_Val (Weibull)', 'Validity']
    row23 = ['', '', round(shape_hat, 4), round(scale_hat, 4), vb, round(b_weibull, 4) if b_weibull else '', va, round(a_weibull, 4) if a_weibull else '', 'Pass' if pass_weibull else 'Fail (Shape<2)']

    final_rows = [row1, row2, row3, [], row5, row6, row7, [], row9, row10, row11, [], row13, row14, row15, [], row17, row18, row19, [], row21, row22, row23]
    
    with open('final_basis_results.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(final_rows)
    
    print(f"Full analysis complete. Selected: {selected_method}")

if __name__ == "__main__":
    run_analysis()