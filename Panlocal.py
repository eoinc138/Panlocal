#!/usr/bin/env python
# coding: utf-8

# In[4]:


#Importing packages for later use.
import csv
import re
import os
import pandas as pd
import matplotlib
from collections import defaultdict

# Querying user for input
column_index = 12  # 13th column is where the EGGNOG annotations are kept for each gene.
search_value = input("Which EGGNOG annotation are you querying?") + "@1"  # This queries the script user for which specific EGGNOG ID they would like to search for in the dataset. The @1 ensures the LUCA EGGNOG annotation is the only one considered.
context_column_index = 15 # This is the column incorporating the contig ID. This ensures the contig ID is considered in the search window.
window_size = int(input("How large of a search window around the anchor gene do you intend to define?")) # This input queries the user on how large of a search window they would like to designate around their anchor gene.

#This function is designed to take the first column of the dataset(biosample/genome ID) and count for each entry the number of instances it occurs.
#This will later be used to dictate the chunk size for reading in the incredibly large dataset. As each chunk will be specific to the genome it contains,
# the search window will not be impeeded from detecting genes outside the chunk.
def count_first_column_entries(input_file_path, output_file_path):
    counts = defaultdict(int) #The dictionary holds the counts for each entry found.
    
    with open(input_file_path, 'r', newline='') as tsvfile: #Tsv file is opened and a reader created to go through it.
        reader = csv.reader(tsvfile, delimiter='\t')

        next(reader, None) #The header is skipped to ensure it is not counted as a chunk size.
        
        for row in reader: #In each row the occurences of a unique value is added up.
            if row: 
                counts[row[0]] += 1
    
    with open(output_file_path, 'w') as outfile: #For later use the counts in order of occurance are written to a text file.
        for entry, count in counts.items():
            outfile.write(f"{count}\n")

#This function is used to assign the counts of each unique Genome ID to a chunk size to be called later. As such each chunk size will correspond
# exactly to the genome it contains.
def get_chunk_sizes(chunk_sizes_file_path):
    with open(chunk_sizes_file_path, 'r') as file: #Opening the text file containing the counts.
        chunk_sizes = [int(line.strip()) for line in file] #Counting the value in each line.
    return chunk_sizes

#This function is specifically designed to read in the dataset in the chunk sizes as outlined by the "counts" file. 
def read_file_in_chunks(file_path, chunk_sizes):
    with open(file_path, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t') #Reading the tsv file in an delimiting it appropriately
        header = next(reader)  #Skipping the initial header of the dataset
        start_index = 0
        for chunk_size in chunk_sizes: 
            chunk = []
            for i, row in enumerate(reader):
                chunk.append(row)
                if len(chunk) == chunk_size:
                    yield chunk
                    chunk = []
                    start_index += chunk_size
                    break
            if chunk:
                yield chunk

def search_EGGNOG_matches(file_path, search_value, window_size, column_index, context_column_index, chunk_sizes):
    results = []
    added_indices = set()

    def find_window(chunk, search_value, window_size, column_index, context_column_index, start_index):
        local_results = []
        local_added_indices = set()
        for i, row in enumerate(chunk):
            global_index = start_index + i
            if search_value in row[column_index]:
                start = max(0, i - window_size)
                end = min(len(chunk), i + window_size + 1)
                search_entry_context_value = row[context_column_index]
                for j in range(start, end):
                    context_row = chunk[j]
                    global_context_index = start_index + j
                    if context_row[context_column_index] == search_entry_context_value and global_context_index not in added_indices and global_context_index not in local_added_indices:
                        local_results.append((global_context_index, context_row))
                        local_added_indices.add(global_context_index)
        return local_results

    start_index = 0
    for chunk in read_file_in_chunks(file_path, chunk_sizes):
        local_results = find_window(chunk, search_value, window_size, column_index, context_column_index, start_index)
        results.extend([context_row for _, context_row in local_results])
        added_indices.update([index for index, _ in local_results])
        start_index += len(chunk)
  
    return results


# Set the working directory to our sample dataset
current_directory = os.getcwd()
file_path = current_directory + '/Finished_test_dataset.tsv'
count_first_column_entries(current_directory + '/Finished_test_dataset.tsv', 'counts.txt')
chunk_sizes_file_path = current_directory + '/counts.txt'

# Function to get the chunk sizes from the file
chunk_sizes = get_chunk_sizes(chunk_sizes_file_path)

# Storing the initial search results.
initial_search = search_EGGNOG_matches(file_path, search_value, window_size, column_index, context_column_index, chunk_sizes)

# Writing the initial results to a tsv file for subsequent analysis.
clean_input_name_1 = re.sub(r'[^\w\-_\. ]|@1', '_', search_value) # Removing any characters that could interfere with the filenaming.
clean_input_name_2 = clean_input_name_1[:-2] # Removing the last two characters for ease of reading multiple files.
window_int = int(window_size) # Incorporating the window size into the file name.
output_file_path = f"OUTPUT_INITIAL_SEARCH_{clean_input_name_2}_WINDOW_{window_int}.tsv" # Final overall file size

# Function to save the results as a tsv file
def saving_results(initial_search, output_file_path):
    with open(output_file_path, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerows(initial_search)

saving_results(initial_search, output_file_path) # Calling the function to write results to a tsv file.

print(f"Data saved to {output_file_path}") # Message to let the user know this portion of the script has been completed and where the results have been saved to.


# In[6]:


from collections import Counter

#This directs the output from out initial search into a variable to be used in the following script.
filename = output_file_path

#This creates an object to store the array from which the previous tsv file will be read in to and analysed.
split_entries = []

#Reading the tsv file
with open(filename, newline='', encoding='utf-8') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    
    for row in reader:
        #The column of interest is column 12 as it contains the EGGNOG annotations
        column_to_split = row[12]
        
        #In order to differentiate the degrees of EGGNOG ID annotation the column is split by the | character.
        split_values = column_to_split.split('|')
        
        #The entries are added to the array
        if split_values and split_values[0].endswith('@1'): #This ensures the entries denoting the EGGNOG LUCA ancestor are the only ones considered.
            split_entries.append(split_values[0])

#The ten most common genes are counted from the dataset and printed showing their order from highest to lowest and the number of instances that they occur.
filtered_entries = [entry[:-2] for entry in split_entries] #Removing the last @1 from each entry for ease of reading
most_common_genes = Counter(filtered_entries).most_common(8)
for result in most_common_genes:
    print(result)

#In order to plot the data the array is saved as a counter object.
counter = Counter(filtered_entries)

with open('output_file_path', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['Gene', 'Count'])  #This denotes the EGGNOG annotation and the number of times it occurs.
    for entry, count in counter.items():
        writer.writerow([entry, count])

#matplotlib is imported to plot the resulting data.
import collections
from matplotlib import pyplot as plt

entries, counts = zip(*most_common_genes)

#Aesthetics of the plotted data.
plt.figure(figsize=(10, 6))
plt.bar(entries, counts, color='skyblue')
plt.xlabel('Genes')
plt.ylabel('Number')
plt.title('8 Most Common Genes in the initial search window')
plt.show()


# In[10]:


#The following section of code is designed to take the list of EGGNOG annotations from the first round or results, count them and match their descriptions with the EGGNOG list of annotations.
#The results are then saved to a tsv file. Having the EGGNOG annotations as a separate file from our master dataset allows for easier incorporation of EGGNOG database
# updates into use of the script.
import pandas as pd
EGGNOG_annotation_file = current_directory + '/EGGNOG_IDS.tsv'
output_file = f"OUTPUT_INITIAL_SEARCH_{clean_input_name_2}_WINDOW_{window_int}_EGGNOG_data.tsv" 

#The number of occurrances in the array are counted.
counter = collections.Counter(filtered_entries)

entry_column_index = 0 #This denotes the first column of our EGGNOG file, the EGGNOG ID.
description_column_index = 1 #This denotes the second column, the EGGNOG descriptions.

#The tsv file of EGGNOG descriptions is read in to an array.
EGGNOG_description_data = {}
with open(EGGNOG_annotation_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)  # Skip the header row
    for row in reader:
        entry = row[entry_column_index]
        info = row[description_column_index]
        EGGNOG_description_data[entry] = info

#The data is sorted by the number of entries in the middle column.
data_to_sort = []
for entry, count in counter.items():
    info = EGGNOG_description_data.get(entry, 'N/A')  #N/A is used where no description is present in the EGGNOG annotation file.
    data_to_sort.append([entry, count, info])

#Sorting the data by descending order.
data_to_sort.sort(key=lambda x: x[1], reverse=True)

#Results are saved to a tsv file with appropriate headers.
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['EGGNOG ID', 'Count', 'Protein name'])
    writer.writerows(data_to_sort)
df = pd.read_csv(output_file, delimiter='\t')
print(df.head(8))
#Printing our results to a tsv file appropriately named.
print(f"Results saved to {output_file}")


# In[12]:


#Querying the user for gene to use in the refined search and how large of a search window to define
second_value = input('What gene would you like to consider as a second search?') + "@1"
windows_size = int(input('How large of a search window would you like to consider around gene 1?'))

#Produce a function that searches through our test dataset, obtains the matches for the EGGNOG annotations and appends them to an array.
def second_search_EGGNOG_matches(file_path, search_value, second_value, window_size, column_index, context_column_index):
    results = []
    all_rows = []
    included_indices = set()
    
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
#Specific flag to ensure our second gene search is being considered in the search window.
            second_value_found = False
            
            #Searching within the window for the second gene of interest.
            for j in range(start, end):
                context_row = all_rows[j]
                if context_row[context_column_index] == search_entry_context_value and second_value in context_row[column_index]:
                    second_value_found = True
                    break
            
            #If the second gene is found in the window the window is then appended to our results.
            if second_value_found:
                for j in range(start, end):
                    if j not in included_indices:
                        context_row = all_rows[j]
                        if context_row[context_column_index] == search_entry_context_value:
                            results.append(context_row)
                            included_indices.add(j)

    return results

#Setting the variables for the second search function.
file_path = current_directory + "/" + output_file_path

#Calling the function to produce a second search and saving the results to an object.
refined_search = second_search_EGGNOG_matches(file_path, search_value, second_value, window_size, column_index, context_column_index)

#Incorporating the second search into the tsv file name.
clean_second_input_name_1 = re.sub(r'[^\w\-_\. ]|@1', '_', second_value) #Removing any characters that could interfere with the filenaming.
clean_second_input_name_2 = clean_second_input_name_1[:-2] #Removing the last two characters for ease of reading multiple files.

window_int = int(window_size) #Incorporating the window size into the file name.

refined_output_file_path = f"OUTPUT_REFINED_SEARCH_{clean_input_name_2}_{clean_second_input_name_2}_WINDOW_{window_int}.tsv" #Final overall file namwe.

#Function to be called to save the results as a tsv file.
def saving_results(refined_search, refined_output_file_path):
    with open(refined_output_file_path, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerows(refined_search)

saving_results(refined_search, refined_output_file_path) #Calling the function to write our results to a tsv file.

print(f"Data saved to {refined_output_file_path}")


# In[14]:


#Ensuring the file from the second search is considered in the subsequent function.
filename = refined_output_file_path

#Again an empty list is created to hold values.
split_entries = []

#Reading in th tsv file.
with open(filename, newline='', encoding='utf-8') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    
    for row in reader:
        #Column 12 holds the EGGNOG annotation information.
        column_to_split = row[12]
        
        #The column is split on the | character.
        split_values = column_to_split.split('|')
        
        #Split values ending in @1 are added to the array.
        if split_values and split_values[0].endswith('@1'): #This ensures only the entries ending in @1 are counted. This is important as our search considers only the LUCA EGGNOG annotation.
            split_entries.append(split_values[0])


#The resulting array is ran through a counter and the 5 most common results are printed.
filtered_entries_2 = [entry[:-2] for entry in split_entries] #Removing the last @1 from each entry for ease of reading
most_common_genes = Counter(filtered_entries_2).most_common(8)
for result in most_common_genes:
    print(result)
import collections
from matplotlib import pyplot as plt

entries, counts = zip(*most_common_genes)

#Plotting the data as a visual bar chart as done before.
plt.figure(figsize=(10, 6))
plt.bar(entries, counts, color='skyblue')
plt.xlabel('Genes')
plt.ylabel('Number')
plt.title('8 Most Common Genes in each window')
plt.show()


# In[64]:


#The following code is a repeat of the earlier EGGNOG tsv file except tailerod to the results of the second, refined search. The results are again saved to a separate, appropriately named file.
EGGNOG_annotation_file = current_directory + '/EGGNOG_IDS.tsv'
output_file = f"OUTPUT_REFINED_SEARCH_{clean_input_name_2}_{clean_second_input_name_2}_WINDOW_{window_int}_EGGNOG_data.tsv" 

#The number of occurrances in the array are counted.
counter = collections.Counter(filtered_entries_2)

entry_column_index = 0 #This denotes the first column of our EGGNOG file, the EGGNOG ID.
description_column_index = 1 #This denotes the second column, the EGGNOG descriptions.

#The tsv file of EGGNOG descriptions is read in to an array.
EGGNOG_description_data = {}
with open(EGGNOG_annotation_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)  # Skip the header row
    for row in reader:
        entry = row[entry_column_index]
        info = row[description_column_index]
        EGGNOG_description_data[entry] = info

#The data is sorted by the number of entries in the middle column.
data_to_sort = []
for entry, count in counter.items():
    info = EGGNOG_description_data.get(entry, 'N/A')  #N/A is used where no description is present in the EGGNOG annotation file.
    data_to_sort.append([entry, count, info])

#Sorting the data by descending order.
data_to_sort.sort(key=lambda x: x[1], reverse=True)

#Results are saved to a tsv file with appropriate headers.
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['EGGNOG ID', 'Count', 'Protein name'])
    writer.writerows(data_to_sort)
df = pd.read_csv(output_file, delimiter='\t')
print(df.head(8))
#Printing our results to a tsv file appropriately named.
print(f"Results saved to {output_file}")


# In[22]:


#The follwing script produces a taxonomic tree from the results of the initial search conducted above.

#The GTDB taxonomy file and refined search file are assigned for later use.
import pandas as pd 
search_file = refined_output_file_path
GTDB_taxonomy_file = current_directory + '/' + 'GTDB_taxonomy.tsv'
#Both files are read into data frames.
df1 = pd.read_csv(search_file, sep='\t')
df2 = pd.read_csv(GTDB_taxonomy_file, sep='\t')

last_column_file1 = df1.iloc[:, -1]

first_column_file2 = df2.iloc[:, 0]

matches = df2[first_column_file2.isin(last_column_file1)] #The two files are organised based on their column of interest and matched then to create a new file.

output = refined_output_file_path + 'GTDB_tree.tsv'
matches.to_csv(output, sep='\t', index=False) #Here the output file is created.



# In[24]:


from anytree import Node, RenderTree, ContRoundStyle
from collections import defaultdict

#This function reads through the output file and counts thhe number of columns(taxonomic classifications). It then prompts the user to select from those available.
def get_max_column(columns):
    print("Available taxonomic classifications")
    for i, col in enumerate(columns):
        print(f"{i + 1}: {col}")
    while True:
        try:
            column_index = int(input("Enter the number corresponding to the taxonomic ranking up to which you'd like to construct the tree: "))
            if 1 <= column_index <= len(columns):
                return column_index
            else:
                print(f"Please enter a number between 1 and {len(columns)}.")
        except ValueError:
            print("Invalid input. Please enter an integer.")

#The tsv file is read in with the first column being ignored as it only contains the matching Biosample ID.
filename = output #Output corresponds to the second search output.

def read_tsv(filename):
    taxonomy_dict = defaultdict(lambda: defaultdict(dict))
    counters = defaultdict(lambda: defaultdict(int)) #This specifically counts the number of entries at each node.

    with open(filename, newline='', encoding='utf-8') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        fieldnames = reader.fieldnames[1:]  #Skipping the first column.
        max_column_index = get_max_column(fieldnames)
        selected_fieldnames = fieldnames[:max_column_index]

        for row in reader:
            current_level = taxonomy_dict
            row_items = list(row.items())[1:max_column_index + 1]  #The first column is skipped and the input specified by the user is called to decide the number of columns used.
            for rank, name in row_items:
                if name not in current_level:
                    current_level[name] = defaultdict(dict)
                counters[rank][name] += 1  
                current_level = current_level[name]
                
    return taxonomy_dict, counters, selected_fieldnames

#The following function creates the tree with matching counters for each entry.
def create_tree(parent, tree_dict, counters, rank_keys, rank_index):
    for key, value in tree_dict.items():
        rank = rank_keys[rank_index] if rank_index < len(rank_keys) else None
        counter = counters[rank][key] if rank else 0
        child = Node(f"{key} ({counter})", parent=parent)
        create_tree(child, value, counters, rank_keys, rank_index + 1)

#The tree roots nodes are created by this function. 
def create_trees(taxonomy_dict, counters, rank_keys):
    roots = []
    for root_name in taxonomy_dict:
        root_rank = rank_keys[0] if rank_keys else None
        counter = counters[root_rank][root_name] if root_rank else 0
        root = Node(f"{root_name} ({counter})")
        create_tree(root, taxonomy_dict[root_name], counters, rank_keys, 1) 
        roots.append(root)
    return roots

#This functions creates the final tree, specifying the style to be used and calling upon the roots created by the create_trees function.
def render_trees(roots):
    for root in roots:
        print(f"Tree for root: {root.name}")
        for pre, _, node in RenderTree(root, style=ContRoundStyle()):
            print("%s%s" % (pre, node.name))
        print("\n")

#Upon creating the tree this function creates a for loop which prompts the user to reproduce the tree with a new taxonomic level if they desire based on readability etc.
def main():
    filename = output
    final_roots = None
    
    while True:
        taxonomy_dict, counters, selected_fieldnames = read_tsv(filename)
        roots = create_trees(taxonomy_dict, counters, selected_fieldnames)
        final_roots = roots  # Save the final tree to an object

        render_trees(roots)
        
        again = input("Would you like to remake the tree? (yes/no): ").strip().lower()
        if again != 'yes':
            break
    
    return final_roots

if __name__ == "__main__":
    final_trees = main()

