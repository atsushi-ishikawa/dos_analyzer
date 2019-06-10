import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

df = pd.read_csv("tmpout.csv")
# df = df[ (df["system"].str.contains("0.50")) | (df["system"].str.contains("1.00")) ]

df  = df[ (df["system"].str.contains("0.50")) ]
df  = pd.concat([df, df["system"].str.split("0.50", expand=True)], axis=1)

df  = df.rename(columns={0:"tmp1", 1:"tmp2"})
df1 = df.loc[:,["tmp1","tmp2","CO-atop"]]
df2 = df.loc[:,["tmp2","tmp1","CO-atop"]]
df1 = df1.rename(columns={"tmp1":"elem1", "tmp2":"elem2"})
df2 = df2.rename(columns={"tmp2":"elem1", "tmp1":"elem2"})

df  = pd.concat([df1,df2])
df2 = pd.pivot_table(df, values="CO-atop", index="elem1", columns="elem2")
df2.fillna(0.0)
sb.set(font="Helvetica")
sb.set_context("notebook", font_scale=1.8)
heatmap = sb.heatmap(df2, linewidth=0.5, cmap="Blues_r", vmin=-2.1, vmax=0.0, \
                        square=True, cbar_kws={"label": "CO adsorption energy (eV)", "format": "%.1f"})
plt.rcParams['figure.figsize'] = 9,9
plt.rcParams['font.family'] = 'Helvetica'
plt.xlabel("")
plt.ylabel("")
plt.show()
