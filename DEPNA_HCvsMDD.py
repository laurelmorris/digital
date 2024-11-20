#!/usr/bin/env python

# Dependency Network Analysis (DEPNA) digital data
# 
# Jacob, Y., Y. Winetraub, G. Raz, E. Ben-Simon, H. Okon-Singer, K. Rosenberg-Katz, T. 
# Hendler and E. Ben-Jacob. Dependency Network Analysis (DEPNA) Reveals Context 
# Related Influence of Brain Network Nodes. Scientific Reports 6: 27444 (2016).
# 
#
# (c) 2024 Yael Jacob


import os
import math
import numpy as np
import pandas as pd
import seaborn as sns  
import matplotlib.pyplot as plt
import scipy.stats as stats
import networkx as nx
import pingouin as pg
from scipy.stats import ttest_ind
from statsmodels.stats import multitest
from matplotlib.lines import Line2D  # Import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable

np.random.seed(42)  # You can choose any seed value 

# choose the columns that you want to include from the spreadsheet
# inputsTypes = ['active','steps_active','passiveTotals_active','screenTotals_active','passive_active','screen_active']
# inputsTypes = ['active_new','steps_active_new','passiveTotals_active_new','screenTotals_active_new','passive_active_new','screen_active_new']
inputsTypes = ['active_sparse','steps_active_sparse','passiveTotals_active_sparse','screenTotals_active_sparse','passive_active_sparse','screen_active_sparse']

for inputsType in inputsTypes:


	if inputsType == 'active':
		indexes = list(range(13, 21)) #8 items
	if inputsType == 'steps_active':
		indexes = list(range(12, 21)) #9 items
	if inputsType == 'passiveTotals_active':
		indexes = list(range(10, 21)) #11 items
	if inputsType == 'screenTotals_active':
		indexes = list(range(10, 12)) + list(range(13, 21))#11 items
	if inputsType == 'passive_active':
		indexes = list(range(2, 10)) + list(range(12, 21)) #17 items 'passive_active'
	if inputsType == 'screen_active':
		indexes = list(range(2, 10)) + list(range(13, 21)) #16 items 'screen_active'

	if inputsType == 'active_new':
		indexes = list(range(13, 18)) + list(range(19, 20)) #8 items
	if inputsType == 'steps_active_new':
		indexes = list(range(12, 18))+ list(range(19, 20))#9 items
	if inputsType == 'passiveTotals_active_new':
		indexes = list(range(10, 18))+ list(range(19, 20)) #11 items
	if inputsType == 'screenTotals_active_new':
		indexes = list(range(10, 12)) + list(range(13,18))+ list(range(19, 20))#11 items
	if inputsType == 'passive_active_new':
		indexes = list(range(2, 10)) + list(range(12, 18))+ list(range(19, 20)) #17 items 'passive_active'
	if inputsType == 'screen_active_new':
		indexes = list(range(2, 10)) + list(range(13, 18))+ list(range(19, 20)) #16 items 'screen_active'

	if inputsType == 'active_sparse':
		indexes = list(range(13, 18))  #5 items
	if inputsType == 'steps_active_sparse':
		indexes = list(range(12, 18))#6 items
	if inputsType == 'passiveTotals_active_sparse':
		indexes = list(range(10, 18)) #7 items
	if inputsType == 'screenTotals_active_sparse':
		indexes = list(range(10, 12)) + list(range(13,18))
	if inputsType == 'passive_active_sparse':
		indexes = list(range(2, 10)) + list(range(12, 18)) 
	if inputsType == 'screen_active_sparse':
		indexes = list(range(2, 10)) + list(range(13, 18)) 


	# FolderPath='DEPNAdigital/DEPNA'
	FolderPath = 'NewProcessed'
	group1 = 'HC' # 'Healthy','PTSD','GAD','SAD',AnxNOS', 'MDD','Anx'
	group2 = 'MA' # 'Healthy','PTSD','MDD','Anx','Patients'

	pdf_file_name = 'DEPNA_LAMP_' + group1 + 'vs' + group2 + '_' + inputsType + '_new.pdf'
	# print(pdf_file_name)


	def calculateSubjAdj(path,indexes):
		# SubjFiles = [i for i in os.listdir(path) if '10' in i ]
		SubjFiles = [i for i in os.listdir(path) if 'U' in i and not i.startswith('.')]

		# SubjFiles = [i for i in os.listdir(path) if os.path.isdir(os.path.join(path, i)) in i]

		print(SubjFiles)
		si=0
		adj_list = []
		adj_batch =[]
		x_list = []
		x_batch =[]
		n=0
		
		
		# nettsdata = [j for j in os.listdir(path+'/'+SubjFiles[0]) if '_000.netts' in j]
		# print(nettsdata)
		for s in SubjFiles:
			print(s)
			# n=n+1
		# for s in range(0,5):
			nettsdata_name = []
			# nettsdata_name = [j for j in os.listdir(path) if  '_screentime_steps_ema.csv' in j and not j.startswith('.')]

			nettsdata_name = s
			# nettsdata_name = [j for j in os.path.isfile(path) if  '_screentime_steps_ema.csv' in j]

			if not nettsdata_name:
				continue

			# n=n+1
			# print(nettsdata_name)
			# read text file into pandas DataFrame
			# nettsdata = pd.read_csv(path+'/'+SubjFiles[s]+ '/' + nettsdata_name[0], sep=" ", header=None)
			nettsdata = []
			nettsdata = pd.read_csv(path+'/'+ nettsdata_name)


			
			
			ROIsLables2use2 = nettsdata.columns[indexes].tolist()

			# Dictionary mapping old strings to new strings
			replacements = {'screen time total': 'Screentime', 'checks total': 'Checks Total', 'steps': 'Steps','anxious':'Anxiety','depressed':'Depression','distressed': 'Distressed','ext_mot': 'Extrinsic Motivation','int_mot': 'Intrinsic Motivation'}

			# Use list comprehension to replace old strings with new strings
			ROIsLables2use = [replacements.get(x, x) for x in ROIsLables2use2]



			n=n+1
			# print(n)
				# if ROIsTC2use.shape[0]==84:
			cols=len(indexes)
					# cols=len(ROIsTC2use)
			# print(cols)

			# ROIsTC2use_transpose=ROIsTC2use.transpose()

			if n==1 :
				SLI_all= np.zeros((len(SubjFiles),cols))
				SLD_all= np.zeros((len(SubjFiles),cols))
				# SubjList= np.zeros((len(SubjFiles)))
				SubjList = [''] * len(SubjFiles)
				allSubjectsDEPNA = np.zeros((len(SubjFiles), cols, cols))
	 
					
			# print(ROIsTC2use_transpose)
			# df =[]
			df = nettsdata.iloc[:,indexes].copy() 
			df.dropna(inplace=True)
			# df = pd.DataFrame(ROIsTC2use_transpose)
					# print(df)

			if len(df) > 4:

				# Z-score normalization
				df_normalized = (df - df.mean()) / df.std()

				corrMatrix = np.zeros((cols,cols))
				corrMatrix = df_normalized.corr()
				corrMatrix[corrMatrix>=1]=0.99
				# corrMatrix[corrMatrix==1]=0.99
				corrMatrix[corrMatrix<=(-1)]=-0.99
				np.fill_diagonal(corrMatrix.values, 1)
				# print('corrMatrix:',corrMatrix)
				corrMatrix2 = corrMatrix.fillna(0)
				# print('corrMatrix2:',corrMatrix2)

				# corrMatrix = np.nan_to_num(corrMatrix)
				# sns.heatmap(corrMatrix, annot=True)
				# plt.show()
				# print(corrMatrix)
				# R = np.zeros((cols,cols))
				R = corrMatrix2
				# R = corrMatrix2.values
				# print(R)
				# print('max:',R.max().max())
				# print('min:',R.min().min())
				# print('R.shape:')
				# print(R.shape)
				# print('cols:')
				# print(cols)
				
				# sns.heatmap(R, annot=False)
				# plt.show()

				important_variable = np.zeros((cols,cols))
				important_variable_std = np.zeros((cols,cols))
				important_variable_mean = np.zeros((cols,cols))
				# important_variable = [[0 for k in xrange(cols)] for j in xrange(cols)]
				Do_z = 1
				diff = 0
				#i = 0 # the exact number of the influencing ROI 
				#j = 1 # the exact number of the influenced ROI
				for i in range(0, cols):
					for j in range(0, cols):
						for ind in range(0, cols):
									
							if (ind!=i) & (ind!=j) & (i!=j):
								corr=0

								# corr=(R(i,j)-R(i,ind)*R(j,ind))/sqrt((1-R(i,ind)^2)*(1-R(j,ind)^2));
								# print('R[i][ind]:',R.iloc[i,ind])
								# print('R[i][j]:',R.iloc[i,j])
								# print('R[j][ind]:',R.iloc[j,ind])
								# first = (R.iloc[i,j]-R.iloc[i,ind]*R.iloc[j,ind])
								# # print('first:',first)
								# second = (1-R.iloc[i,ind]**2)*(1-R.iloc[j,ind]**2)
								# print('second:',second)
								corr=(R.iloc[i,j]-R.iloc[i,ind]*R.iloc[j,ind])/(math.sqrt((1-R.iloc[i,ind]**2)*(1-R.iloc[j,ind]**2)))
								# print('corr:',corr)
								if np.isinf(corr):
									corr = 0.99 if corr > 0 else -0.99

								if corr>=1:
									corr = 0.99
								if corr <=(-1):
									corr = -0.99

								data_z=0
								corr_z=0
								diff=0
								# print(corr)
								if Do_z==1:
									# data_z=.5*math.log((1+R[i][j])/(1-R[i][j]))
									# print('R.iloc[i,j]:',R.iloc[i,j])
									data_z=math.atanh(R.iloc[i,j]) #0.5*(log(1+d_R) - log(1-d_R)) #atanh (R)
									# print('data_z:',data_z)
									# print('corr2:',corr)
									corr_z=math.atanh(corr) #0.5*(log(1+corr) - log(1-corr)) #atanh (R)
									# corr_z=.5*math.log((1+corr)/(1-corr))
									diff=data_z-corr_z
									#print(diff)
								else:
									diff=R[i,j]-corr
									# diff=corr;
								if diff < 0:
									# important_variable[i][j][ind] = 0
									important_variable[j,ind] = 0
								else:
									important_variable[j,ind] = diff
									# important_variable[i][j][ind] = diff
							#print(ind)
					for t in range(0, cols):
						if (i!=t):
							temp = important_variable[:,t]
							# print(temp)
							# temp = np.delete(temp, [i,t])
							# print(temp)
							important_variable_std[i,t]=np.nanstd(temp)
							important_variable_mean[i,t]=np.nanmean(temp)
						else:
							important_variable_std[i,t]=0
							important_variable_mean[i,t]=0	 

						# print(important_variable)
						# print('important_variable[1][2][4]')
						# print(important_variable[1][2][4])
						# print('important_variable[2][1][4]')
						# print(important_variable[2][1][4])
						# print('important_variable[4][1][2]')
						# print(important_variable[4][1][2])
						# print('important_variable[1][4][2]')
						# print(important_variable[1][4][2])

				np.fill_diagonal(important_variable_mean, 0)
				print('important_variable.shape:')
				print(important_variable.shape)



						# important_variable_std = np.zeros((cols,cols))
						# important_variable_mean = np.zeros((cols,cols))
						# for i in range(0, cols):
						# 	for t in range(0, cols):
						# 		if (i!=t):
						# 		#print(important_variable[i][:][t])
						# 			temp = important_variable[i][:][t]
						# 			# print(temp)
						# 			temp = np.delete(temp, [i,t])
						# 			# print(temp)
						# 			important_variable_std[i][t]=np.nanstd(temp)
						# 			important_variable_mean[i][t]=np.nanmean(temp)
						# 		else:
						# 			important_variable_std[i][t]=0
						# 			important_variable_mean[i][t]=0
						# sns.heatmap(important_variable_mean, annot=False)
						# plt.show()
				print(important_variable_mean)
				print('important_variable_mean.shape:')
				print(important_variable_mean.shape)
				SLI = np.nansum(important_variable_mean, axis=0)
				# SLI = SLI.transpose
				# print('SLI:',SLI)

				print(SLI.shape)
				# np.append(SLI_all,SLI)
				SLI_all[n-1,:]=np.copy(SLI)
				# print(SLI_all)

				SLD = np.nansum(important_variable_mean.transpose(),axis=0)
				# SLD = SLD.transpose
				# print('SLD:',SLD)
				# np.append(SLD_all,SLD)
				SLD_all[n-1,:]=np.copy(SLD)
				# print('n:', n)
				SubjList[n-1]=s

				allSubjectsDEPNA[n-1, :, :] = important_variable_mean


				# df_depna = pd.DataFrame(important_variable_mean)
				# df_depna.columns=ROIsLables2use

				# # Create a directed graph
				# G = nx.DiGraph()

				# # Add nodes to the graph
				# G.add_nodes_from(df_depna.columns)

				# # Add directed edges based on the top 10% of highest connections
				# edges = []
				# for i in range(len(df_depna.columns)):
				# 	for j in range(len(df_depna.columns)):
				# 		if i != j:  # Exclude self-loops
				# 			correlation = df_depna.iloc[i, j]
				# 			edges.append((df_depna.columns[i], df_depna.columns[j], correlation))

				# # Sort edges based on the absolute value of weights in descending order
				# edges.sort(key=lambda x: abs(x[2]), reverse=True)

				# # Select the top 10% of edges
				# top_10_percent = int(0.1 * len(edges))
				# selected_edges = edges[:top_10_percent]

				# # Add selected edges to the graph
				# G.add_weighted_edges_from(selected_edges)

				# # Draw the directed graph with edge thickness based on weight
				# pos = nx.circular_layout(G)  # You can choose different layout algorithms
				# labels = {node: node for node in G.nodes()}
				# edge_weights = [abs(G[edge[0]][edge[1]]['weight']) * 3 for edge in G.edges()]  # Adjust the multiplier for edge thickness

				# nx.draw(G, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color='red', connectionstyle='arc3,rad=0.1')

				# plt.title('DEPNA Directed Graph with Top 10% of Highest Connections')
				# plt.show()
			else:
				# Move to the next subject
				print("Number of rows not larger than 4, moving to the next subject...")


		return(SLI_all, SLD_all,ROIsLables2use, SubjList,allSubjectsDEPNA)






	[SLI_all_1, SLD_all_1,ROIsLables2use, SubjList_1,allSubjectsDEPNA_1] = calculateSubjAdj(FolderPath +'/'+ group1,indexes)

	SLI_all=pd.DataFrame(SLI_all_1)
	SLI_all.columns = ROIsLables2use

	SLI_all.insert(0, 'SubjList', SubjList_1)
	SLI = SLI_all
	SLI_all_1 = SLI_all_1[(SLI.iloc[:, 1:] != 0).any(axis=1)]
	SLI_all = SLI_all[(SLI_all.iloc[:, 1:] != 0).any(axis=1)]
	# SLI_all = SLI_all[SLI_all.iloc[:, 0].notnull()]
	SLI_all.to_csv('SLI_totals_'+ group1 + '_' + inputsType + '_new.csv', index=False)

	SLD_all=pd.DataFrame(SLD_all_1)
	SLD_all.columns = ROIsLables2use
	SLD_all.insert(0, 'SubjList', SubjList_1)
	SLD = SLD_all
	SLD_all_1 = SLD_all_1[(SLD.iloc[:, 1:] != 0).any(axis=1)]
	SLD_all = SLD_all[(SLD_all.iloc[:, 1:] != 0).any(axis=1)]
	# SLD_all = SLD_all[SLD_all.iloc[:, 0].notnull()]
	SLD_all.to_csv('SLD_totals_'+ group1 + '_' + inputsType + '_new.csv', index=False)

	allSubjectsDEPNA_1= allSubjectsDEPNA_1[(SLI.iloc[:, 1:] != 0).any(axis=1)]
	# allSubjectsDEPNA_1 = allSubjectsDEPNA_1[SLI_all.iloc[:, 0].notnull()]


	[SLI_all_2, SLD_all_2,ROIsLables2use, SubjList_2,allSubjectsDEPNA_2] = calculateSubjAdj(FolderPath +'/'+ group2,indexes)


	SLI_all=pd.DataFrame(SLI_all_2)
	SLI_all.columns = ROIsLables2use

	SLI_all.insert(0, 'SubjList', SubjList_2)
	SLI = SLI_all
	SLI_all_2 = SLI_all_2[(SLI.iloc[:, 1:] != 0).any(axis=1)]
	# SLI_all = SLI_all[(SLI_all[:, 1:] != 0).any(axis=1)]
	SLI_all = SLI_all[(SLI_all.iloc[:, 1:] != 0).any(axis=1)]
	# print(SLI_all.iloc[1, 0].notnull())
	# SLI_all = SLI_all[SLI_all.iloc[:, 0].notnull()]

	SLI_all.to_csv('SLI_totals_'+ group2 + '_' + inputsType + '_new.csv', index=False)

	# SLI_all_MDD=pd.DataFrame(SLI_all_2)
	# SLI_all_MDD.columns = ROIsLables2use
	# SLI_all_MDD.to_csv('SLI_all_MDD.csv', index=False)

	SLD_all=pd.DataFrame(SLD_all_2)
	SLD_all.columns = ROIsLables2use
	SLD_all.insert(0, 'SubjList', SubjList_2)
	SLD = SLD_all
	SLD_all_2 = SLD_all_2[(SLD.iloc[:, 1:] != 0).any(axis=1)]
	SLD_all = SLD_all[(SLD_all.iloc[:, 1:] != 0).any(axis=1)]
	# SLD_all = SLD_all[SLD_all.iloc[:, 0].notnull()]
	SLD_all.to_csv('SLD_totals_'+ group2 + '_' + inputsType + '_new.csv', index=False)

	# SLD_all_MDD=pd.DataFrame(SLD_all_2)
	# SLD_all_MDD.columns = ROIsLables2use
	# SLD_all_MDD.to_csv('SLD_all_MDD.csv', index=False) 

	allSubjectsDEPNA_2= allSubjectsDEPNA_2[(SLI.iloc[:, 1:] != 0).any(axis=1)]
	# allSubjectsDEPNA_2 = allSubjectsDEPNA_2[SLI_all.iloc[:, 0].notnull()]

	# Assuming you have 3D matrices for healthy and patient groups
	# Replace these with your actual data
	healthy_data = allSubjectsDEPNA_1  # Replace with your actual dimensions
	patient_data = allSubjectsDEPNA_2  # Replace with your actual dimensions

	# Get the dimensions of the matrices
	num_subjects, i_size, j_size = healthy_data.shape

	# Initialize an array to store t_values and p-values for each element [i, j]
	t_values = np.zeros((i_size, j_size))
	p_values = np.zeros((i_size, j_size))
	cohen_d = np.zeros((i_size, j_size))

	# Perform t-test for each element [i, j]
	for i in range(i_size):
	    for j in range(j_size):
	        healthy_values = healthy_data[:, i, j]
	        patient_values = patient_data[:, i, j]

			# Perform t-test
	        t_statistic, p_value = ttest_ind(healthy_values, patient_values)

	        # Store the t-value and p-value
	        t_values[i, j] = t_statistic
	        p_values[i, j] = p_value
	        # Calculate Cohen's d value
	        cohen_d[i, j] = pg.compute_effsize(healthy_values, patient_values, eftype='cohen')



	# average_DEPNA = np.mean(allSubjectsDEPNA_1, axis=0)


	
	labelsNames = np.array([ROIsLables2use])
	df_depna = pd.DataFrame(t_values)
	df_depna.columns=ROIsLables2use

	df_depna_t = pd.DataFrame(t_values)
	df_depna_t.columns=ROIsLables2use
	df_depna_t.insert(0, 'Labels', labelsNames.T.flatten())
	df_depna_t.to_csv('df_depna_t_'+ group1 + 'vs' + group2 + '_' + inputsType + '_new.csv', index=False)


	df_depna_p = pd.DataFrame(p_values)
	df_depna_p.columns=ROIsLables2use
	df_depna_p.insert(0, 'Labels', labelsNames.T.flatten())
	df_depna_p.to_csv('df_depna_p_'+ group1 + 'vs' + group2 + '_' + inputsType + '_new.csv', index=False)

	df_depna_d = pd.DataFrame(cohen_d)
	df_depna_d.columns=ROIsLables2use
	df_depna_d.insert(0, 'Labels', labelsNames.T.flatten())
	df_depna_d.to_csv('df_depna_d_'+ group1 + 'vs' + group2 + '_' + inputsType + '_new.csv', index=False)


	# Add directed edges based on the p-values smaller than 0.05
	# for i in range(i_size):
	# 	print(f'p_values  {i} :', p_values[i, :])
	# 	# Perform FDR correction (Benjamini-Hochberg procedure)
	# 	rejected, corrected_p_values, _, _ = multitest.multipletests(p_values[i, :], method='fdr_bh')
	# 	print(f'corrected_p_values  {i} :' ,rejected)
	# 	for j in range(j_size):
	# 		# if i != j and p_values[i, j] < (0.05/j_size): #Bonferroni
	# 		if i != j and p_values[i, j] < 0.05: #uncorrected
	# 		# if i != j and rejected[j] == 'True':
	# 			t_value = t_values[i, j]
	# 			# Determine the direction of the arrow based on the sign of the t-value
	# 			if t_value > 0:
	# 				arrow_color = 'blue'  # Healthy > Patients
	# 			else:
	# 				arrow_color = 'red'   # Patients > Healthy

	# 			G.add_edge(df_depna.columns[i], df_depna.columns[j], weight=t_value, arrow_color=arrow_color)


	# for i in range(i_size):
	#     # print(f'p_values {i}:', p_values[i, :])

	#     # Identify NaN indices
	#     not_nan_indices = ~np.isnan(p_values[i, :])

	#     # Perform FDR correction (Benjamini-Hochberg procedure)
	#     rejected, _, _, _ = multitest.multipletests(p_values[i, not_nan_indices], method='fdr_bh')

	#     # Update rejected array in place for non-NaN positions
	#     rejected_with_nan = np.full_like(p_values[i, :], False)
	#     rejected_with_nan[not_nan_indices] = rejected
	#     print(f'rejected {i} :', rejected)
	#     print(f'rejected_with_nan {i} :', rejected_with_nan)

	#     for j in range(j_size):
	#         # Skip NaN values
	#         if i != j and not np.isnan(p_values[i, j]) and rejected_with_nan[j]:
	#             t_value = t_values[i, j]

	#             # Determine the direction of the arrow based on the sign of the t-value
	#             arrow_color = 'blue' if t_value > 0 else 'red'

	#             G.add_edge(df_depna.columns[i], df_depna.columns[j], weight=t_value, arrow_color=arrow_color)


	with PdfPages(pdf_file_name) as pdf:

		# Create a directed graph
		seed = 13648  # Seed random number generators for reproducibility

		G = nx.DiGraph()

		# Add nodes to the graph
		G.add_nodes_from(df_depna.columns)
		plt.figure(figsize=(10, 6))  # Adjust the figure size as needed

		for i in range(i_size):
			not_nan_indices = ~np.isnan(p_values[i, :])
			# Perform FDR correction (Benjamini-Hochberg procedure)
			rejected, _, _, _ = multitest.multipletests(p_values[i, not_nan_indices], method='fdr_bh')
			# Update rejected array in place for non-NaN positions
			rejected_with_nan = np.full_like(p_values[i, :], False)
			rejected_with_nan[not_nan_indices] = rejected

			for j in range(j_size):
				if i != j and not np.isnan(p_values[i, j]) and rejected_with_nan[j] and t_values[i, j] <0 :
					t_value = abs(t_values[i, j])
					# Determine the direction of the arrow based on the sign of the t-value
					# arrow_color = 'blue' if t_value > 0 else 'red'
					arrow_color = 'dimgray'
					# G.add_edge(df_depna.columns[i], df_depna.columns[j], weight=t_value, arrow_color=arrow_color) #Old
					G.add_edge(df_depna.columns[j], df_depna.columns[i], weight=t_value, arrow_color=arrow_color) #New



		# pos = nx.kamada_kawai_layout(G)
		# pos = nx.spring_layout(G, seed=seed)
		pos = nx.circular_layout(G)  # You can choose different layout algorithms
		labels = {node: node for node in G.nodes()}
		edge_weights = [abs(G[edge[0]][edge[1]]['weight']) * 1 for edge in G.edges()]  # Adjust the multiplier for edge thickness
		arrow_colors = [G[edge[0]][edge[1]]['arrow_color'] for edge in G.edges()]


		nx.draw_networkx_edges(G, pos, alpha=1, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')
		# # nx.draw(G, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')
		
		for ii in range(i_size):
			for jj in range(j_size):
				if ii != jj and p_values[ii, jj] <= 0.05 and t_values[ii, jj] < 0 :
					t_value = abs(t_values[ii, jj])
					# Determine the direction of the arrow based on the sign of the t-value
					# arrow_color = 'lightblue' if t_value > 0 else 'lightpink'
					arrow_color = 'lightgray'
					# G.add_edge(df_depna.columns[ii], df_depna.columns[jj], weight=t_value, arrow_color=arrow_color) #Old
					G.add_edge(df_depna.columns[jj], df_depna.columns[ii], weight=t_value, arrow_color=arrow_color) #New
		

		# Draw the directed graph with edge thickness based on weight and arrow color
		# pos = nx.circular_layout(G)  # You can choose different layout algorithms
		labels = {node: node for node in G.nodes()}
		edge_weights = [abs(G[edge[0]][edge[1]]['weight']) * 1 for edge in G.edges()]  # Adjust the multiplier for edge thickness
		arrow_colors = [G[edge[0]][edge[1]]['arrow_color'] for edge in G.edges()]
		
		nx.draw_networkx_edges(G, pos, alpha=0.5, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')
		
		# nx.draw(G, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10)




		SLI_t_values = np.zeros((i_size))
		SLI_p_values = np.zeros((i_size))
		SLI_cohen_d = np.zeros((i_size))
		SLI_healthy_Av_values = np.zeros((i_size))
		SLI_patient_Av_values = np.zeros((i_size))
		SLI_healthy_STD_values = np.zeros((i_size))
		SLI_patient_STD_values = np.zeros((i_size))

		# Perform t-test for each element [i]
		for i in range(i_size):
			SLI_healthy_values = SLI_all_1[:, i]
			SLI_patient_values = SLI_all_2[:, i]
			SLI_healthy_Av_values[i] = np.nanmean(SLI_all_1[:, i])
			SLI_patient_Av_values[i] = np.nanmean(SLI_all_2[:, i])
			SLI_healthy_STD_values[i] = np.nanstd(SLI_all_1[:, i])
			SLI_patient_STD_values[i] = np.nanstd(SLI_all_2[:, i])
			# Perform t-test
			SLI_t_statistic, SLI_p_value = ttest_ind(SLI_healthy_values, SLI_patient_values)
			# Store the t-value and p-value
			SLI_t_values[i] = SLI_t_statistic
			SLI_p_values[i] = SLI_p_value
			# Calculate Cohen's d value
			SLI_cohen_d[i] = pg.compute_effsize(SLI_healthy_values, SLI_patient_values, eftype='cohen')



		# Perform FDR correction (Benjamini-Hochberg procedure)
		SLI_rejected, SLI_corrected_p_values, _, _ = multitest.multipletests(SLI_p_values, method='fdr_bh')
		
		labelsNames = np.array([ROIsLables2use])
		print(labelsNames.T)

		
		# print('labelsNames:',labelsNames.T.flatten())
		# print('SLI_healthy_Av_values:', SLI_healthy_Av_values.T.flatten())
		# print('SLI_healthy_STD_values:', SLI_healthy_STD_values.T.flatten())
		# print('SLI_patient_Av_values:', SLI_healthy_STD_values.T.flatten())
		# print('SLI_patient_STD_values:', SLI_healthy_STD_values.T.flatten())
		# print('SLI_t_values:', SLI_t_values.T.flatten())
		# print('SLI_p_values:', SLI_p_values.T.flatten())
		# print('SLI_corrected_p_values:', SLI_corrected_p_values.T.flatten())



		df_SLI_t = pd.DataFrame({
    		'Label': labelsNames.T.flatten(),  # Assuming labelsNames is a 2D array
    		'mean healthy': SLI_healthy_Av_values.T.flatten(),
    		'std healthy': SLI_healthy_STD_values.T.flatten(),
    		'mean patients': SLI_patient_Av_values.T.flatten(),
    		'std patients': SLI_patient_STD_values.T.flatten(),
    		't-value': SLI_t_values.T.flatten(),
    		'p-value': SLI_p_values.T.flatten(),
    		'FDR adjusted p-value': SLI_corrected_p_values.T.flatten(),
    		'Cohen d': SLI_cohen_d.T.flatten()
		})


		# df_SLI_t.columns=['Label','mean HC','std HC','mean patients','std patients','t-value','p-value', 'FDR adjusted p-value']
		df_SLI_t.to_csv('df_SLI_t_'+ group1 + 'vs' + group2 + '_' + inputsType + '_new.csv', index=False)

		# Calculate significant t-tests

		# SLI_significant_ttests = SLI_p_values < (0.05/i_size) #Bonferroni
		SLI_FDRsignificant_ttests = SLI_corrected_p_values < 0.05 #FDR
		SLI_significant_ttests = SLI_p_values < 0.05 #Uncorrected

		# Draw the nodes

		# node_colors = SLI_t_values * 1000
		colorsT=(-1)*SLI_t_values
		# Apply colormap with transparency
		cmap = plt.cm.jet  # Choose your desired colormap
		norm = plt.Normalize(colorsT.min(), colorsT.max())  # Normalize the values
		rgba_colors = cmap(norm(colorsT))  # Get RGBA colors from the colormap
		# rgba_colors = cmap(abs(norm))  # Get RGBA colors from the colormap

		# Set transparency level (adjust as needed)
		alpha = 1 #0.7  # Change this value to set the desired transparency level

		# Modify the RGBA colors to include transparency
		rgba_colors[:, 3] = alpha

		# Draw the graph with nodes colored using the modified RGBA colors
		nx.draw_networkx_nodes(G, pos, node_size=400, node_color=rgba_colors)



		# nx.draw_networkx_nodes(G, pos, node_size=500, node_color=abs(SLI_t_values), cmap=plt.cm.jet)
		
		# Specify offset for label positions (adjust this value as needed)
		# label_pos = {node: (x, y + 0.05) for node, (x, y) in pos.items()}

		# # Draw the labels outside of the nodes
		# label_pos = {k: (v[0], v[1] + 0.1) for k, v in pos.items()}  # Adjust the label positions
		# nx.draw_networkx_labels(G, label_pos, labels=labels, font_color='black', font_size=10)

		# Divide nodes into upper and lower halves
		upper_nodes = set(n for n, p in nx.circular_layout(G).items() if p[1] > 0)
		lower_nodes = set(G.nodes()) - upper_nodes

		# # Draw the nodes with a circular layout
		# pos = nx.circular_layout(G)

		# # Draw the nodes
		# nx.draw_networkx_nodes(G, pos)

		# # Draw the edges
		# nx.draw_networkx_edges(G, pos)

		# Draw the labels above for upper half and below for lower half
		label_pos = {}
		for node, (x, y) in pos.items():
		    if node in upper_nodes:
		        label_pos[node] = (x, y + 0.15)  # Adjust label position above
		    else:
		        label_pos[node] = (x, y - 0.15)  # Adjust label position below

		nx.draw_networkx_labels(G, label_pos, labels=labels, font_color='black', font_size=16)




		# # Draw the labels
		# nx.draw_networkx_labels(G, pos, labels=labels, font_color='black', font_size=10)

		# nx.draw(G, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')


		# # Create custom legend entries
		# legend_entries = [Line2D([0], [0], color='blue', lw=2, label='Healthy > Patients'),
		#                   Line2D([0], [0], color='red', lw=2, label='Patients > Healthy')]



		# # nx.draw(G, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')

		# plt.legend(handles=legend_entries, loc='upper left', bbox_to_anchor=(1, 1))  # Add legend
		plt.title(f'{group2}  >  {group1}')
		plt.axis('off')

		# Create a fake ScalarMappable object for the colorbar
		sm = ScalarMappable(cmap=cmap, norm=norm)
		sm.set_array([])  # Set dummy array

		# Show the colorbar
		cbar = plt.colorbar(sm, fraction=0.05, pad=0.05)  # Adjust fraction and pad as needed
		cbar.ax.yaxis.set_tick_params(labelsize=12)  # Change font size
		# Set custom tick labels to represent transparency levels
		# cbar.ax.set_yticklabels(['Low', 'High'])

		# Add a label to the colorbar
		cbar.set_label('Influencing Degree t-value')

		# plt.colorbar(sm)
		# plt.colorbar()
		# plt.tight_layout()
		pdf.savefig()
		# plt.show()
		plt.close()

		# Create a directed graph Healthy > Patients
		seed = 13648  # Seed random number generators for reproducibility

		G2 = nx.DiGraph()

		# Add nodes to the graph
		G2.add_nodes_from(df_depna.columns)
		
		plt.figure(figsize=(10, 6))  # Adjust the figure size as needed

		for i in range(i_size):
			not_nan_indices = ~np.isnan(p_values[i, :])
			# Perform FDR correction (Benjamini-Hochberg procedure)
			rejected, _, _, _ = multitest.multipletests(p_values[i, not_nan_indices], method='fdr_bh')
			# Update rejected array in place for non-NaN positions
			rejected_with_nan = np.full_like(p_values[i, :], False)
			rejected_with_nan[not_nan_indices] = rejected

			for j in range(j_size):
				if i != j and not np.isnan(p_values[i, j]) and rejected_with_nan[j] and t_values[i, j] > 0:
					t_value = t_values[i, j]
					# Determine the direction of the arrow based on the sign of the t-value
					# arrow_color = 'blue' if t_value > 0 else 'red'
					arrow_color = 'dimgray'
					G.add_edge(df_depna.columns[i], df_depna.columns[j], weight=t_value, arrow_color=arrow_color)


		# pos = nx.kamada_kawai_layout(G2)
		# pos = nx.spring_layout(G2, seed=seed)
		pos = nx.circular_layout(G2)  # You can choose different layout algorithms
		labels = {node: node for node in G2.nodes()}
		edge_weights = [abs(G2[edge[0]][edge[1]]['weight']) * 1 for edge in G2.edges()]  # Adjust the multiplier for edge thickness
		arrow_colors = [G2[edge[0]][edge[1]]['arrow_color'] for edge in G2.edges()]


		nx.draw_networkx_edges(G2, pos, alpha=1, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')
		# # nx.draw(G2, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')
		
		
		for ii in range(i_size):
			for jj in range(j_size):
				if ii != jj and p_values[ii, jj] <= 0.05 and t_values[ii, jj] > 0 :
					t_value = t_values[ii, jj]
					# Determine the direction of the arrow based on the sign of the t-value
					arrow_color = 'lightgray'
					G2.add_edge(df_depna.columns[ii], df_depna.columns[jj], weight=t_value, arrow_color=arrow_color)



		# Draw the directed graph with edge thickness based on weight and arrow color
		# pos = nx.circular_layout(G2)  # You can choose different layout algorithms
		labels = {node: node for node in G2.nodes()}
		edge_weights = [abs(G2[edge[0]][edge[1]]['weight']) * 1 for edge in G2.edges()]  # Adjust the multiplier for edge thickness
		arrow_colors = [G2[edge[0]][edge[1]]['arrow_color'] for edge in G2.edges()]
		
		nx.draw_networkx_edges(G2, pos, alpha=0.5, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')
		
		# nx.draw(G2, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10)




		SLI_t_values = np.zeros((i_size))
		SLI_p_values = np.zeros((i_size))

		# Perform t-test for each element [i]
		for i in range(i_size):
			SLI_healthy_values = SLI_all_1[:, i]
			SLI_patient_values = SLI_all_2[:, i]
			# Perform t-test
			SLI_t_statistic, SLI_p_value = ttest_ind(SLI_healthy_values, SLI_patient_values)
			# Store the t-value and p-value
			SLI_t_values[i] = SLI_t_statistic
			SLI_p_values[i] = SLI_p_value

		# Perform FDR correction (Benjamini-Hochberg procedure)
		SLI_rejected, SLI_corrected_p_values, _, _ = multitest.multipletests(SLI_p_values, method='fdr_bh')

		# Calculate significant t-tests

		# SLI_significant_ttests = SLI_p_values < (0.05/i_size) #Bonferroni
		SLI_FDRsignificant_ttests = SLI_corrected_p_values < 0.05 #FDR
		SLI_significant_ttests = SLI_p_values < 0.05 #Uncorrected

		# Draw the nodes

		# node_colors = SLI_t_values * 1000
		colorsT=SLI_t_values
		# Apply colormap with transparency
		cmap = plt.cm.jet  # Choose your desired colormap
		norm = plt.Normalize(colorsT.min(), colorsT.max())  # Normalize the values
		rgba_colors = cmap(norm(colorsT))  # Get RGBA colors from the colormap
		# rgba_colors = cmap(abs(norm))  # Get RGBA colors from the colormap

		# Set transparency level (adjust as needed)
		alpha = 1 #0.7  # Change this value to set the desired transparency level

		# Modify the RGBA colors to include transparency
		rgba_colors[:, 3] = alpha

		# Draw the graph with nodes colored using the modified RGBA colors
		nx.draw_networkx_nodes(G2, pos, node_size=400, node_color=rgba_colors)



		# nx.draw_networkx_nodes(G2, pos, node_size=500, node_color=abs(SLI_t_values), cmap=plt.cm.jet)
		
		# Specify offset for label positions (adjust this value as needed)
		# label_pos = {node: (x, y + 0.05) for node, (x, y) in pos.items()}
		
		# # Draw the labels outside of the nodes
		# label_pos = {k: (v[0], v[1] + 0.1) for k, v in pos.items()}  # Adjust the label positions
		# nx.draw_networkx_labels(G2, label_pos, labels=labels, font_color='black', font_size=10)

		# Divide nodes into upper and lower halves
		upper_nodes = set(n for n, p in nx.circular_layout(G2).items() if p[1] > 0)
		lower_nodes = set(G.nodes()) - upper_nodes

		# # Draw the nodes with a circular layout
		# pos = nx.circular_layout(G)

		# # Draw the nodes
		# nx.draw_networkx_nodes(G, pos)

		# # Draw the edges
		# nx.draw_networkx_edges(G, pos)

		# Draw the labels above for upper half and below for lower half
		label_pos = {}
		for node, (x, y) in pos.items():
		    if node in upper_nodes:
		        label_pos[node] = (x, y + 0.15)  # Adjust label position above
		    else:
		        label_pos[node] = (x, y - 0.15)  # Adjust label position below

		nx.draw_networkx_labels(G2, label_pos, labels=labels, font_color='black', font_size=16)


		# Draw the labels
		# nx.draw_networkx_labels(G2, pos, labels=labels, font_color='black', font_size=10)

		# nx.draw(G2, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')


		# # Create custom legend entries
		# legend_entries = [Line2D([0], [0], color='blue', lw=2, label='Healthy > Patients'),
		#                   Line2D([0], [0], color='red', lw=2, label='Patients > Healthy')]



		# # nx.draw(G2, pos, with_labels=True, labels=labels, node_size=700, node_color='skyblue', font_color='black', font_size=10, arrowsize=15, width=edge_weights, edge_color=arrow_colors, connectionstyle='arc3,rad=0.1')

		# plt.legend(handles=legend_entries, loc='upper left', bbox_to_anchor=(1, 1))  # Add legend

		plt.title(f'{group1}  >  {group2}')
		plt.axis('off')

		# Create a fake ScalarMappable object for the colorbar
		sm = ScalarMappable(cmap=cmap, norm=norm)
		sm.set_array([])  # Set dummy array

		# Show the colorbar
		cbar = plt.colorbar(sm, fraction=0.05, pad=0.05)  # Adjust fraction and pad as needed
		cbar.ax.yaxis.set_tick_params(labelsize=12)  # Change font size
		# Set custom tick labels to represent transparency levels
		# cbar.ax.set_yticklabels(['Low', 'High'])

		# Add a label to the colorbar
		cbar.set_label('Influencing Degree t-value')

		# plt.colorbar(sm)
		# plt.colorbar()
		# plt.tight_layout()
		pdf.savefig()
		# plt.show()
		plt.close()



		print(SLI_significant_ttests)
		# Plot the bar plot with averages and mark significant t-tests
		bar_width = 0.35
		index = np.arange(i_size)

		# Calculate mean and SEM for each group
		SLI_mean_healthy = np.mean(SLI_all_1, axis=0)
		SLI_sem_healthy = np.std(SLI_all_1, axis=0) #/ np.sqrt(SLI_all_1.shape[0])

		SLI_mean_patient = np.mean(SLI_all_2, axis=0)
		SLI_sem_patient = np.std(SLI_all_2, axis=0) #/ np.sqrt(SLI_all_2.shape[0])

		# plt.bar(index, SLI_mean_healthy, yerr=SLI_sem_healthy, width=bar_width,capsize=5, label=f'{group1} (N = {len(SLI_all_1)})', color='lightskyblue')
		# plt.bar(index + bar_width, SLI_mean_patient, yerr=SLI_sem_patient, width=bar_width, capsize=5,label=f'{group2} (N = {len(SLI_all_2)})', color='lightcoral')
		plt.bar(index, SLI_mean_healthy, width=bar_width, label=f'{group1} (N = {len(SLI_all_1)})', color='#999999')
		plt.bar(index + bar_width, SLI_mean_patient, width=bar_width, label=f'{group2} (N = {len(SLI_all_2)})', color='#99CCFF')


		# Customize error bar line width
		plt.errorbar(index, SLI_mean_healthy, yerr=SLI_sem_healthy, fmt='none', linewidth=0.5, capsize=3, capthick=0.5, color='black')
		plt.errorbar(index + bar_width, SLI_mean_patient, yerr=SLI_sem_patient, fmt='none', linewidth=0.5, capsize=3, capthick=0.5, color='black')

		# Mark significant t-tests with '*'
		# for i, significant in enumerate(SLI_significant_ttests):
		#     if significant:
		for i, (significant, FDRsignificant) in enumerate(zip(SLI_significant_ttests, SLI_FDRsignificant_ttests)):
			if significant and not FDRsignificant:
				plt.text(index[i] + bar_width / 2, max(SLI_mean_healthy[i] + SLI_sem_healthy[i], SLI_mean_patient[i] + SLI_sem_patient[i]) + 0.05, '*', ha='center', va='bottom', color='black')

		# Mark FDR significant t-tests with '**'
		for i, significant in enumerate(SLI_FDRsignificant_ttests):
		    if significant:
		        plt.text(index[i] + bar_width / 2, max(SLI_mean_healthy[i] + SLI_sem_healthy[i], SLI_mean_patient[i] + SLI_sem_patient[i]) + 0.05, '**', ha='center', va='bottom', color='black')
		
		plt.ylim(0, 1)
		# plt.xlabel('Element [i]')
		plt.ylabel('Influencing Degree')
		#plt.title('Influencing Degree')
		plt.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

		# Set x-axis labels
		plt.xticks(index + bar_width / 2, [f'{i}' for i in ROIsLables2use], rotation=45, ha='right')  # Replace with your actual labels
		# Remove the box around the bars and leave only x and y axis
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.rcParams.update({'font.size': 14})  # Change 12 to your desired font size
		plt.tight_layout()
		pdf.savefig()
		# plt.show()
		plt.close()



		SLD_t_values = np.zeros((i_size))
		SLD_p_values = np.zeros((i_size))
		SLD_cohen_d = np.zeros((i_size))
		SLD_healthy_Av_values = np.zeros((i_size))
		SLD_patient_Av_values = np.zeros((i_size))
		SLD_healthy_STD_values = np.zeros((i_size))
		SLD_patient_STD_values = np.zeros((i_size))

		# Perform t-test for each element [i]
		for i in range(i_size):
			SLD_healthy_values = SLD_all_1[:, i]
			SLD_patient_values = SLD_all_2[:, i]
			SLD_healthy_Av_values[i] = np.nanmean(SLD_all_1[:, i])
			SLD_patient_Av_values[i]= np.nanmean(SLD_all_2[:, i])
			SLD_healthy_STD_values[i] = np.nanstd(SLD_all_1[:, i])
			SLD_patient_STD_values[i] = np.nanstd(SLD_all_2[:, i])
		
			# Perform t-test
			SLD_t_statistic, SLD_p_value = ttest_ind(SLD_healthy_values, SLD_patient_values)
			# Store the t-value and p-value
			SLD_t_values[i] = SLD_t_statistic
			SLD_p_values[i] = SLD_p_value
			# Calculate Cohen's d value
			SLD_cohen_d[i] = pg.compute_effsize(SLD_healthy_values, SLD_patient_values, eftype='cohen')


		# Perform FDR correction (Benjamini-Hochberg procedure)
		SLD_rejected, SLD_corrected_p_values, _, _ = multitest.multipletests(SLD_p_values, method='fdr_bh')
		print('SLD corrected_p_values :' ,SLD_corrected_p_values)

		# df_SLD_t = pd.DataFrame([labelsNames.T, 
		# 	SLD_healthy_Av_values.T,
		# 	SLD_healthy_STD_values.T,
		# 	SLD_patient_Av_values.T,
		# 	SLD_patient_STD_values.T,
		# 	SLD_t_values.T, 
		# 	SLD_p_values.T,
		# 	SLD_corrected_p_values.T])
		
		
		# df_SLD_t.columns=['Label','mean HC','std HC','mean patients','std patients','t-value','p-value', 'FDR adjusted p-value']
		df_SLD_t = pd.DataFrame({
    		'Label': labelsNames.T.flatten(),  # Assuming labelsNames is a 2D array
    		'mean healthy': SLD_healthy_Av_values.T.flatten(),
    		'std healthy': SLD_healthy_STD_values.T.flatten(),
    		'mean patients': SLD_patient_Av_values.T.flatten(),
    		'std patients': SLD_patient_STD_values.T.flatten(),
    		't-value': SLD_t_values.T.flatten(),
    		'p-value': SLD_p_values.T.flatten(),
    		'FDR adjusted p-value': SLD_corrected_p_values.T.flatten(),
    		'Cohen d': SLD_cohen_d.T.flatten()
    	})

		df_SLD_t.to_csv('df_SLD_t_'+ group1 + 'vs' + group2 + '_' + inputsType + '_new.csv', index=False)





		# Calculate significant t-tests
		# SLD_significant_ttests = SLD_p_values < (0.05/i_size) #Bonferroni
		SLD_FDRsignificant_ttests = SLD_corrected_p_values < 0.05 #FDR
		SLD_significant_ttests = SLD_p_values < 0.05 #Uncorrected

		print(SLD_significant_ttests)
		# Plot the bar plot with averages and mark significant t-tests
		bar_width = 0.35
		index = np.arange(i_size)

		# Calculate mean and SEM for each group
		SLD_mean_healthy = np.mean(SLD_all_1, axis=0)
		SLD_sem_healthy = np.std(SLD_all_1, axis=0) #/ np.sqrt(SLD_all_1.shape[0])

		SLD_mean_patient = np.mean(SLD_all_2, axis=0)
		SLD_sem_patient = np.std(SLD_all_2, axis=0) #/ np.sqrt(SLD_all_2.shape[0])

		# plt.bar(index, SLD_mean_healthy, yerr=SLD_sem_healthy, width=bar_width,capsize=5,label=f'{group1} (N = {len(SLD_all_1)})', color='lightskyblue')
		# plt.bar(index + bar_width, SLD_mean_patient, yerr=SLD_sem_patient, width=bar_width, capsize=5, label=f'{group2} (N = {len(SLD_all_2)})', color='lightcoral')

		plt.bar(index, SLD_mean_healthy, width=bar_width,label=f'{group1} (N = {len(SLD_all_1)})', color='#999999')
		plt.bar(index + bar_width, SLD_mean_patient, width=bar_width, label=f'{group2} (N = {len(SLD_all_2)})', color='#99CCFF')


		# Customize error bar line width
		plt.errorbar(index, SLD_mean_healthy, yerr=SLD_sem_healthy, fmt='none', linewidth=0.5, capsize=3, capthick=0.5, color='black')
		plt.errorbar(index + bar_width, SLD_mean_patient, yerr=SLD_sem_patient, fmt='none', linewidth=0.5, capsize=3, capthick=0.5, color='black')

		# Mark significant t-tests with '*'
		#for i, significant in enumerate(SLD_significant_ttests):
		    #if significant:
		for i, (significant, FDRsignificant) in enumerate(zip(SLD_significant_ttests, SLD_FDRsignificant_ttests)):
			if significant and not FDRsignificant:
				plt.text(index[i] + bar_width / 2, max(SLD_mean_healthy[i] + SLD_sem_healthy[i], SLD_mean_patient[i] + SLD_sem_patient[i]) + 0.05, '*', ha='center', va='bottom', color='black')

		# Mark significant t-tests with '**'
		for i, significant in enumerate(SLD_FDRsignificant_ttests):
		    if significant:
		        plt.text(index[i] + bar_width / 2, max(SLD_mean_healthy[i] + SLD_sem_healthy[i], SLD_mean_patient[i] + SLD_sem_patient[i]) + 0.05, '**', ha='center', va='bottom', color='black')

		plt.ylim(0, 1)

		# plt.xlabel('Element [i]')
		plt.ylabel('Influenced Degree')
		#plt.title('Influenced Degree')
		plt.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

		# Set x-axis labels
		plt.xticks(index + bar_width / 2, [f'{i}' for i in ROIsLables2use], rotation=45, ha='right')  # Replace with your actual labels
		# Remove the box around the bars and leave only x and y axis
		plt.gca().spines['top'].set_visible(False)
		plt.gca().spines['right'].set_visible(False)
		plt.rcParams.update({'font.size': 14})  # Change 12 to your desired font size
		plt.tight_layout()
		pdf.savefig()
		# plt.show()
		plt.close()

