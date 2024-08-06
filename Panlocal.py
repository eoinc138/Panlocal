#!/usr/bin/env python
# coding: utf-8

# In[44]:


#Importing packages for later use.
import csv
import re
import os
import pandas as pd
import collections
from collections import defaultdict
import io

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


#Set the working directory to our sample dataset
current_directory = os.getcwd()
file_path = current_directory + '/Finished_test_dataset.tsv'
count_first_column_entries(current_directory + '/Finished_test_dataset.tsv', 'counts.txt')
chunk_sizes_file_path = current_directory + '/counts.txt'

#Function to get the chunk sizes from the file
chunk_sizes = get_chunk_sizes(chunk_sizes_file_path)

#Storing the initial search results.
initial_search = search_EGGNOG_matches(file_path, search_value, window_size, column_index, context_column_index, chunk_sizes)

#Writing the initial results to a tsv file for subsequent analysis.
clean_input_name_1 = re.sub(r'[^\w\-_\. ]|@1', '_', search_value) # Removing any characters that could interfere with the filenaming.
clean_input_name_2 = clean_input_name_1[:-2] # Removing the last two characters for ease of reading multiple files.
window_int = int(window_size) # Incorporating the window size into the file name.
output_file_path = f"OUTPUT_INITIAL_SEARCH_{clean_input_name_2}_WINDOW_{window_int}.tsv" # Final overall file size

#Function to save the results as a tsv file
def saving_results(initial_search, output_file_path):
    with open(output_file_path, 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerows(initial_search)

saving_results(initial_search, output_file_path) #Calling the function to write results to a tsv file.

print(f"Data saved to {output_file_path}") #Message to let the user know this portion of the script has been completed and where the results have been saved to.


# In[46]:


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
from matplotlib import pyplot as plt

entries, counts = zip(*most_common_genes)
PDF_figure_initial_search_output = 'Initial search figure.png'
fig = plt.figure(figsize=(10, 6))
plt.bar(entries, counts, color='skyblue')
plt.xlabel('Genes')
plt.ylabel('Number')
plt.title('8 Most Common Genes in the initial search window')
plt.savefig(PDF_figure_initial_search_output)
plt.show()


# In[48]:


#EGGNOG and COG description files are read in and assigned to statements, additionally the output file is assigned with a unique name corresponding to the search performed
EGGNOG_annotation_file = current_directory + '/EGGNOG_IDS.tsv'
COG_annotation_file = current_directory + '/COGnarCOG_descriptors.tsv'
output_file = f"OUTPUT_INITIAL_SEARCH_{clean_input_name_2}_WINDOW_{window_int}_EGGNOG_data.tsv"

#The entries from the first search are counted and stored as an object
counter = collections.Counter(filtered_entries)
#Column indices are assigned for later calling
ID_column_index = 0  #EGGNOG/COG ID column
GENE_column_index = 1  #EGGNOG description column/COG gene name column
DESCRIPTION_column_index = 2  #Cog description column. Note that what EGGNOG considers the name column corresponds to the COG description column.

#A dictionary is created to store the gene description results.
description_data = {}

#Read the COG annotation file assigning the ID to entry, Gene name to description and additional info to the description column.
with open(COG_annotation_file, 'r',encoding='ISO 8859-1') as f: #Characters in the COG file required a different form of encoding.
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)  #Skip the header row as it is not required.
    for row in reader:
        entry = row[ID_column_index]
        gene = row[GENE_column_index]
        additional_info = row[DESCRIPTION_column_index]
        description_data[entry] = (gene, additional_info)
        
#Read in the EGGNOG annotation file assigning ID to the entry column and Gene name to the description column.
with open(EGGNOG_annotation_file, 'r') as eggnog_file:
    eggnog_reader = csv.reader(eggnog_file, delimiter='\t') #No encoding required
    eggnog_header = next(eggnog_reader)  #Skipping the header row again.
    for row in eggnog_reader:
        entry = row[ID_column_index]
        description = row[GENE_column_index]
        if entry not in description_data:
            description_data[entry] = ('N/A', description) #Leave blank the gene name column.

#To a dataframe the data from the description data dictionary is parsed for matches in our counter for results. Matches have their EGGNOG ID, number of occurances 
# Gene name and description added where possible
data_to_sort = []
for entry, count in counter.items():
    description, additional_info = description_data.get(entry, ('N/A', 'N/A'))  #N/A designates where no information was present for the specific ID.
    data_to_sort.append([entry, count, description, additional_info])
    
#The data frame is sorted based on the counts
data_to_sort.sort(key=lambda x: x[1], reverse=True)

#The dataframe is then written to a tsv file for future reference.
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['EGGNOG ID', 'Count', 'Gene name', 'Description'])
    writer.writerows(data_to_sort)

#The first 8 rows are shown in the terminal for the user.
initial_results = pd.read_csv(output_file, delimiter='\t')
PDF_content_initial_results = initial_results.head(8) #Saving the top 8 results for later writing to a pdf
print(initial_results.head(8))

#Message informing the user of where the file has been saved.
print(f"Results saved to {output_file}")

#The following section of code queries the user on specific gene they would like to see in the results.
results = pd.DataFrame()#Dataframe to hold the users results.
with open(output_file, 'r'): #Using the results file.
    while True: 
        entry = input("Enter the EGGNOG ID or Gene name to see specific details:")
        result = initial_results[initial_results.apply(lambda row: row.astype(str).str.contains(entry).any(), axis=1)] #Considers the entire row, allowing the user to Query by EGGNOG ID, gene name or descriptive key words.
        if not result.empty:
            print("Found entries:")
            print(result)
            results = pd.concat([results,result], ignore_index=True) #Appending results to the dataframe.
        else:
            print("Gene was not present in results") #In the event of the gene not being found.
        again = input("Would you like to search again? (yes/no): ").strip().lower()
        if again != 'yes':
            print("Search concluded, summary of queries below. Proceeding with script.") #Message to inform the user of the scripts' progress.
            break

print(results)#Displaying results in terminal once more for the user

PDF_content_queried_initial_results = results#Saving results for later writing to pdf.


# In[54]:


#Produce a function that searches through our test dataset, obtains the matches for the EGGNOG annotations and appends them to an array.
#Note that for this search, reading in chunks is no longer necessary due to the small size of the tsv file produced by the initial search.
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
repeat_process = True
#The following code is indented to facilitate a prompt for the user to repeat the refined search using a different gene. This is due to the length of time
# it would take to run the analysis from the very beginning
while repeat_process:
    #Querying the user for gene to use in the refined search and how large of a search window to define
    second_value = input('What gene would you like to consider as a second search?') + "@1"
    windows_size = int(input('How large of a search window would you like to consider around gene 1?'))

    #Setting the variables for the second search function.
    file_path = current_directory + "/" + output_file_path
    
    #Calling the function to produce a second search and saving the results to an object.
    refined_search = second_search_EGGNOG_matches(file_path, search_value, second_value, window_size, column_index, context_column_index)
    
    #Incorporating the second search into the tsv file name.
    clean_second_input_name_1 = re.sub(r'[^\w\-_\. ]|@1', '_', second_value) #Removing any characters that could interfere with the filenaming.
    clean_second_input_name_2 = clean_second_input_name_1[:-2] #Removing the last two characters for ease of reading multiple files.
    window_int = int(window_size) #Incorporating the window size into the file name.
    
    #Final overall file name.
    refined_output_file_path = f"OUTPUT_REFINED_SEARCH_{clean_input_name_2}_{clean_second_input_name_2}_WINDOW_{window_int}.tsv" 
    
    #Calling the function to write our results to a tsv file.
    saving_results(refined_search, refined_output_file_path) 
    #Message to let the user know the file location.
    print(f"Data saved to {refined_output_file_path}")
    
    #Ensuring the file from the second search is considered in the subsequent code.
    filename = refined_output_file_path
    
    #Again an empty list is created to hold values.
    split_entries = []
    
    #Reading in the tsv file.
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
    
    #Saving the entries and their corresponding counts for plotting by matplotlib.
    entries, counts = zip(*most_common_genes)
    
    PDF_figure_refined_search_output = 'Refined search figure.png' #Saving our figure so it can later be written to a tsv file.
    fig = plt.figure(figsize=(10, 6))
    plt.bar(entries, counts, color='skyblue')
    plt.xlabel('Genes')
    plt.ylabel('Number')
    plt.title('8 Most Common Genes in the refined search window')
    plt.savefig(PDF_figure_refined_search_output)
    plt.show()
    
    #The following code is a repeat of the earlier EGGNOG tsv file except tailerod to the results of the second, refined search. The results are again saved to a separate, appropriately named file.
    EGGNOG_annotation_file = current_directory + '/EGGNOG_IDS.tsv'
    COG_annotation_file = current_directory + '/COGnarCOG_descriptors.tsv'
    output_file = f"OUTPUT_REFINED_SEARCH_{clean_input_name_2}_{clean_second_input_name_2}_WINDOW_{window_int}_EGGNOG_data.tsv" 
    
    #The number of occurrances in the array are counted.
    counter = collections.Counter(filtered_entries_2)
    
    #Note that the indexes for the specific columns in our description file have already been set and no longer need to be called.
    #A dictionary is created to store the gene description results.
    description_data = {}

    #Read the COG annotation file assigning the ID to entry, Gene name to description and additional info to the description column.
    with open(COG_annotation_file, 'r',encoding='ISO 8859-1') as f: #Characters in the COG file required a different form of encoding.
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)  #Skip the header row as it is not required.
        for row in reader:
            entry = row[ID_column_index]
            gene = row[GENE_column_index]
            additional_info = row[DESCRIPTION_column_index]
            description_data[entry] = (gene, additional_info)
            
    #Read in the EGGNOG annotation file assigning ID to the entry column and Gene name to the description column.
    with open(EGGNOG_annotation_file, 'r') as eggnog_file:
        eggnog_reader = csv.reader(eggnog_file, delimiter='\t') #No encoding required
        eggnog_header = next(eggnog_reader)  #Skipping the header row again.
        for row in eggnog_reader:
            entry = row[ID_column_index]
            description = row[GENE_column_index]
            if entry not in description_data:
                description_data[entry] = ('N/A', description) #Leave blank the gene name column.
        
        #To a dataframe the data from the description data dictionary is parsed for matches in our counter for results. Matches have their EGGNOG ID, number of occurances 
        # Gene name and description added where possible
        data_to_sort = []
        for entry, count in counter.items():
            description, additional_info = description_data.get(entry, ('N/A', 'N/A'))  #N/A designates where no information was present for the specific ID.
            data_to_sort.append([entry, count, description, additional_info])
    
    #The data frame is sorted based on the counts
    data_to_sort.sort(key=lambda x: x[1], reverse=True)
    
    #The dataframe is then written to a tsv file for future reference.
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['EGGNOG ID', 'Count', 'Gene name', 'Description'])
        writer.writerows(data_to_sort)
    
    #The top 8 results are written shown in the terminal to the user
    refined_results = pd.read_csv(output_file, delimiter='\t')
    PDF_content_refined_results = refined_results.head(8) #Saved for writing to a pdf
    print(refined_results.head(8))
    
    #Printing our results to a tsv file appropriately named.
    print(f"Results saved to {output_file}")
    
    #The code below is a repeat of the previous code for the purpose of assess genes of interest for the user.
    results = pd.DataFrame()
    with open(output_file, 'r'):
        while True: 
            entry = input("Enter the EGGNOG ID or Gene name to see specific details:")
            result = refined_results[refined_results.apply(lambda row: row.astype(str).str.contains(entry).any(), axis=1)]
            if not result.empty:
                print("Found entries:")
                print(result)
                results = pd.concat([results,result], ignore_index=True)
            else:
                print("Gene was not present in results")
            again = input("Would you like to search again? (yes/no): ").strip().lower()
            if again != 'yes':
                print("Search concluded, summary of queries below. Proceeding with script.")
                break
    
    print(results)
    
    PDF_content_queried_refined_results = results #Saving results for later writing to pdf.
    
    
    #The follwing script produces a taxonomic tree from the results of the initial search conducted above.
    #The GTDB taxonomy file and refined search file are assigned for later use.
    search_file = refined_output_file_path
    GTDB_taxonomy_file = current_directory + '/GTDB_taxonomy.tsv'
    
    #Both files are read into data frames.
    df1 = pd.read_csv(search_file, sep='\t')
    df2 = pd.read_csv(GTDB_taxonomy_file, sep='\t')
    #The refined search output file is sorted by the last column.
    last_column_file1 = df1.iloc[:, -1]
    #The GTDB taxonomy file is sorted by the first column.
    first_column_file2 = df2.iloc[:, 0]
    #The two files now organised on the biosample column are matched based on this column, duplicates are removed to account for multiple loci.
    matches = df2[first_column_file2.isin(last_column_file1)].drop_duplicates()
    #Here the output files are named for the resulting matches.
    output1 = f"OUTPUT_REFINED_SEARCH_{clean_input_name_2}_{clean_second_input_name_2}_WINDOW_{window_int}_GTDB_tree.tsv"
    output2 = f"OUTPUT_REFINED_SEARCH_{clean_input_name_2}_{clean_second_input_name_2}_WINDOW_{window_int}_GTDB_results.tsv"
    
    matches.to_csv(output1, sep='\t', index=False) #Here the output file is created and written as defined by the output.
    
    #Additionally a permanent file is created which removes the biosample column leaving it possible to later analyse the number of species represented in the search. 
    matches.iloc[:, 1:].to_csv(output2, sep='\t', index=False) 

    #Anytree is a python package designed to enable the conversion of tsv files into renderings as graphic trees.
    from anytree import Node, RenderTree, ContRoundStyle
    
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
    filename = output1 #Output corresponds to the second search output.
    
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
    
    #The tree root nodes are created by this function. 
    def create_trees(taxonomy_dict, counters, rank_keys):
        roots = []
        for root_name in taxonomy_dict:
            root_rank = rank_keys[0] if rank_keys else None
            counter = counters[root_rank][root_name] if root_rank else 0
            root = Node(f"{root_name} ({counter})") #Specifies the 
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
    
    #Upon creating the tree this function creates a loop which prompts the user to reproduce the tree with a new taxonomic level if they desire based on readability etc.
    def main():
        filename = output1
        final_roots = None
        
        while True:
            taxonomy_dict, counters, selected_fieldnames = read_tsv(filename)
            roots = create_trees(taxonomy_dict, counters, selected_fieldnames)
            final_roots = roots  
    
            render_trees(roots)
            
            again = input("Would you like to remake the tree? (yes/no): ").strip().lower()
            if again != 'yes':
                break
        
        return final_roots
    
    if __name__ == "__main__":
        final_trees = main()

    #Sys is required here for exiting the script if the user desires.
    import sys
    #This function queries the user on wether they would like to produce a PDF summary of the bar plots and tables created.
    #Note that as creating pdfs in python can be occasionally difficult and depend on packages which may break in future. As such the user can exit 
    #the script at this point, using screenshots from the terminal for reference.
    def main():
        response = input("Would you like to produce a summary pdf of the results found? Warning this may encounter errors arrising from packages. (yes/no): ").strip().lower()
        
        if response == 'yes':
            print("Producing summary pdf")
            pdf_output = 'SUMMARY_' + clean_input_name_2 + '_' + clean_second_input_name_2 + '.pdf'
            print(f"PDF saved as {pdf_output}")
        elif response == 'no':
            print("No pdf produced") #Message to inform user of script proceeding
    
    if __name__ == "__main__":
        main()
    
    #Reportlab is a python package designed primarily to create pdfs from python scripts. It does not require any additional packages/dependancies.
    from reportlab.lib.pagesizes import letter
    from reportlab.lib import colors
    from reportlab.platypus import Image, SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    styles = getSampleStyleSheet()

    #This function when called later will wrap the text to be inserted in our table, in the vast majority of cases all text will be displayed in a clear legible manner.
    def wrap_text(text, style):
        return Paragraph(text, style)
    #This function takes our dataframes as inputs and converts them to a basic table format.
    def create_table_from_df(df):
        #Here the dataframe entries are converted to paragraphs
        data_for_table = [[wrap_text(str(cell), styles['BodyText']) for cell in row] for row in [df.columns.to_list()] + df.values.tolist()]
    
        #The table is created with a dynamic row height related to the len("size") of the paragraphs from the dataframe.
        table = Table(data_for_table, rowHeights=[0.9*inch]*len(data_for_table)) 
        #This specific table style was chosen for aim of the readability, the header section is clearly distinguished by the Grey background while each entry is readable against the white.
        #This was chosen with the aim of allowing these tables to be easily transfered/inserted into articles, documents etc.
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('ALIGN',(0,0),(-1,-1),'CENTER'),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),  #Text in all sections is aligned to the top to allow for readability.
        ])
    
        table.setStyle(style)
        return table
    
    #The following function appends our dataframe(now table) to the list which will be made into a PDF. It allows for a changing title for each table.
    def add_dataframe_to_elements(elements, df, title):
        elements.append(Paragraph(title, styles['Heading2']))
        if df.empty:
            elements.append(Paragraph("No queried genes/descriptions were found in the search window.", styles['BodyText'])) #In the event of none of the user's queries being found a message is added to the pdf informing them of this.
        else:
            elements.append(create_table_from_df(df))
            elements.append(Spacer(1, 20))
            elements.append(PageBreak()) #Inserting a page break ensures that when possible each table will not be divided by the placement between pages.
    
    #The following function takes the saved Matplotlib figures and appends them to our list.
    def add_plot_to_elements(elements, plot_path, title, width=8 * inch, height=6 * inch):
        img = Image(plot_path)
        img.drawHeight = height #Both height and width were chosen with readability in mind when converted to pdf format.
        img.drawWidth = width  
        elements.append(img)
        elements.append(Spacer(1, 20))
        elements.append(Paragraph(title, styles['Heading2']))
        elements.append(PageBreak()) #Again a page break is inserted to prevent disruption of following figures/tables.
    
    #PDF is named dynamically based on previous inputs given. Although pdf_output is already defined earlier, deleting it here prevents script from completing, will remove if time allows.
    pdf_output = 'SUMMARY_' + clean_input_name_2 + '_' + clean_second_input_name_2 + '.pdf'
    #This set the template and output path of our pdf
    pdf = SimpleDocTemplate(pdf_output, pagesize=letter)
    
    #The elements list is empty pending the appending of the contents of this most recent search.
    elements = []
    
    #Step by step the tables/figures are added by the add_... functions. Initially a dynamic title is added based on the search carried out.
    #Each figure/table is then added by calling their add.. function which call the table/figure creating functions within them.
    elements.append(Paragraph("Summary for " + clean_input_name_2 + " and " + clean_second_input_name_2 + " Search", styles['Heading2']))
    add_dataframe_to_elements(elements, PDF_content_initial_results, "Table 1: Initial Search results")
    add_plot_to_elements(elements, PDF_figure_initial_search_output, "Figure 1: Initial Search results plot")
    add_dataframe_to_elements(elements, PDF_content_queried_initial_results, "Table 2: Queried Initial Search results")
    add_dataframe_to_elements(elements, PDF_content_refined_results, "Table 3: Refined Search results")
    add_plot_to_elements(elements, PDF_figure_refined_search_output, "Figure 2: Refined Search results plot")
    add_dataframe_to_elements(elements, PDF_content_queried_refined_results, "Table 4: Queried Refined Search results")
    
    #This last line builds the pdf from the elements list once it contains our data.
    pdf.build(elements)

    #Here the user is queried on wether they would like to repeat their search, stating yes reverts to the beginning of the refined search. Stating no exits the script.
    repeat = input("Would you like to repeat the refined search using a different gene? (yes/no): ").strip().lower()
    if repeat != 'yes':
        #These last few lines remove any unneeded files from the directory, leaving only the OUTPUT and SUMMARY files for later use.
        os.remove('output_file_path')
        os.remove('Initial search figure.png')
        os.remove('Refined search figure.png')
        os.remove('counts.txt')
        print('Panlocal script has concluded')
        sys.exit()


# In[ ]:




