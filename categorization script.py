import pandas as pd



# for i, row in data.iterrows():
#      if row["Allele_Frequency"] <= 0.001:
#          row["Allele_Frequency"]=1
#      elif row['Allele_Frequency'] > 0.001 and row['Allele_Frequency'] <= 0.005:
#          row["Allele_Frequency"] = 2
#      elif row['Allele_Frequency'] > 0.005 and row['Allele_Frequency'] <= 0.05:
#          row["Allele_Frequency"] = 3
#      else :
#          row["Allele_Frequency"] = 4
#
# data['Allele_Frequency'] = data['Allele_Frequency'].astype('category')
# #data.apply(categorize, axis=1).head()
#
# category = { "Allele_frequency" : {1.0:"VERY RARE", 2.0:"RARE",3.0:"LOW",4.0:"COMMON"}}
# data.replace(category,inplace=True)


data = pd.read_csv("Allele Frequency.csv")
maxval=data['Allele_Frequency'].max()

labels=["RARE","LOW","COMMON"]
cut_points=[0,0.005,0.05] +[maxval]

data = pd.cut(data['Allele_Frequency'],bins=cut_points,labels=labels,include_lowest=True)
data.to_csv("Allele frequency_categorized", sep='\t')