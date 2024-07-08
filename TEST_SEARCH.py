#!/usr/bin/env python
# coding: utf-8

# In[34]:


import csv

#Set the working directory to our sample dataset

file_path = '/home/eoin/Working_Directory/sorted_data.tsv'

#Produce a function that searches through our test dataset, obtains the matches for the EGGNOG annotations and appends them to an array.
def search_EGGNOG_matches(file_path, column_index, search_value, context_column_index = 15, window_size = 5):
    results = []
    all_rows = []
    
    with open(file_path, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        header = next(reader)
        
        all_rows = [row for row in reader]
#This for loop specifically designates the search window surrounding our anchor gene.
    for i, row in enumerate(all_rows):
        if search_value in row[column_index]:
            start = max(0, i - window_size)
            end = min(len(all_rows), i + window_size + 1)
            search_entry_context_value = row[context_column_index]
#This for loop limits the results being appended to only include those on the same contig.
            for j in range(start, end):
                context_row = all_rows[j]
                if context_row[context_column_index] == search_entry_context_value:
                    results.append(context_row)
  
    return results

# Query user for input
column_index = 12  # 13th column is where the EGGNOG annotations are kept for each gene. 
search_value = input("Which EGGNOG annotation are you querying?")  # This queries the script user for which specific EGGNOG ID they would like to search for in the dataset.
context_column_index = 15 

#Printing the initial search results
initial_search = search_EGGNOG_matches(file_path, column_index, search_value)
for row in initial_search:
    print(row)



# In[ ]:




