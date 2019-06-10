import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np

#from pylab import rcParams

plt.rcParams['figure.figsize'] = 19,4
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.size'] = 10
plt.rcParams['xtick.labelsize'] = 10-2

df  = pd.read_csv("tmpout.csv")
df2 = df.sort_values("CO-atop")
df2.index = df2["system"].str.split("_",expand=True)[0]
df2.plot.bar(y="CO-atop", legend=None, color="cadetblue")
plt.xlabel("")
plt.ylabel("CO adsorption energy (eV)")
plt.xticks(rotation=60,ha="right")
plt.yticks(np.linspace(-2.0, 0, 5))
plt.ylim([-2.2, 0.0])
plt.savefig("all_adsorption.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
plt.show()
df_all  = df2.loc[:,["system","CO-atop"]]
df_all.to_csv("all_adsorption.csv")

# making pure metal csv file
df_pure = df2[ df2["system"].str.contains("1.00") ]
df_pure = df_pure.loc[:,["system","CO-atop"]]
df_pure["system"] = df_pure["system"].str.replace("_111","")
df_pure.set_index("system")
del df_pure["system"]
df_pure.to_csv("puremetals.csv")
