# -*- coding: utf-8 -*-
"""
Created on Sun May  7 20:33:49 2023

@author: Allison Walker
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

def readActivities(activities, activities_filebase):
    probabilities = {}
    paths = {}
    genome_bgcs = {}
    bgc_to_genome_dict = {}
    for activity in activities:
        probabilities[activity] = {}
        infile = open(activities_filebase + activity + ".txt")
        for line in infile:
            split_line = line.split(",")
            cluster_path = split_line[0]
            cluster = split_line[1]
            cluster = cluster[0:cluster.rfind(".")]
            genome = cluster[0:cluster.rfind("/")]
            if ".gbff" in genome:
                genome = genome[0:genome.rfind(".")]
                cluster_part1 = cluster[0:cluster.find("/")]
                cluster_part2 = cluster[cluster.find("/"):len(cluster)]
                cluster = genome + cluster_part2
            if "_copy" in genome:
                genome = genome[0:genome.rfind("_")]
                cluster_part1 = cluster[0:cluster.find("/")]
                cluster_part2 = cluster[cluster.find("/"):len(cluster)]
                cluster = genome + cluster_part2
            if "genomic_2" in genome:
                genome = genome[0:genome.rfind("_")]
                cluster_part1 = cluster[0:cluster.find("/")]
                cluster_part2 = cluster[cluster.find("/"):len(cluster)]
                cluster = genome + cluster_part2
            if genome not in genome_bgcs:
                genome_bgcs[genome] = []
            if activity == "antibacterial" and cluster not in genome_bgcs[genome]:
                genome_bgcs[genome].append(cluster)
            bgc_to_genome_dict[cluster] = genome
            probability = float(split_line[2])
            if cluster not in paths:
                paths[cluster] = cluster_path
            if cluster not in probabilities[activity]:
                probabilities[activity][cluster] = probability
    return probabilities, paths, genome_bgcs, bgc_to_genome_dict

def readTaxonomyFile(taxonomy_file):
    genus_members = {}
    genus_list = []
    full_taxonomy = {}
    for line in taxonomy_file:
        if "# Genome" in line:
            continue
        split_line = line.split("\t")
        kingdom = split_line[1]
        phylum = split_line[2]
        tax_class = split_line[3]
        order = split_line[4]
        family = split_line[5]
        genus = split_line[6]
        species = split_line[7]
        genome = split_line[0].replace("/","")
        
        if ".gbff" in genome:
            genome = genome[0:genome.rfind(".")] 
        if "_copy" in genome:
            genome = genome[0:genome.rfind("_")]
        if "genomic_2" in genome:
            genome = genome[0:genome.rfind("_")]
        if genus not in genus_list:
            genus_list.append(genus)
            genus_members[genus] = []
        if genome not in genus_members[genus]:
            genus_members[genus].append(genome)
        full_taxonomy[genome] = [kingdom, phylum, tax_class, order, family, genus, species]
    return genus_members, genus_list, full_taxonomy

def getClusterCounts(bgcs, probabilities, activity, characterized_gcfs, l2_gcfs, activity_threshold):
    total_counts = 0
    uncharacterized_counts = 0
    characterized_counts = 0
    for bgc in bgcs:
        if bgc not in probabilities[activity]:
            print("HERE!")
            return -1, -1, -1
        prob = probabilities[activity][bgc]
        if prob < activity_threshold:
            continue
        total_counts += 1
        if bgc not in l2_gcfs:
            return -1, -1, -1
        gcf = l2_gcfs[bgc]
        characterized = characterized_gcfs[gcf]
        if characterized:
            characterized_counts += 1
        else:
            uncharacterized_counts += 1
    return total_counts, characterized_counts, uncharacterized_counts



#make this the input directory with the output from the l2 norm clustering script
#all outputs will also be written to this directory
input_dir = ""
#which activity type to analyze options are "antibacterial", "antigrampos","antigramneg", "antieuk","antifungal","antitumor"
activity_to_plot = "antibacterial"

l2_gcf_file = open(input_dir+"membership")
original_gcf_file = open(input_dir+"bgc_names.txt")
taxonomy_file = open("bigslice_input/uncharacterized_clusters/bigscape_taxonomy.tsv")
no_bgcs_list_file = open("genomes_without_bgcs.txt")
activities_filebase = "predictions_"
activities = ["antibacterial", "antigrampos","antigramneg", "antieuk","antifungal","antitumor"]
print("reading activity predictions")
probabilities, paths, genome_bgcs, bgc_to_genome_dict = readActivities(activities, activities_filebase)

print("reading taxonomy")
genus_members, genus_list, full_taxonomy = readTaxonomyFile(taxonomy_file)
original_gcfs = {}
l2_gcfs = {}
bgc_names = {}
characterized_gcfs = {}
gcf_activities= {}
list_of_bgcs_in_l2gcf = {}
genomes_without_bgcs = []
genomes_in_bigslice_results = []
for line in no_bgcs_list_file:
    genome = line.replace("\n","")
    if genome not in genomes_without_bgcs:
        genomes_without_bgcs.append(genome)
        genomes_in_bigslice_results.append(line.replace("\n",""))
    
print("reading gcfs")
for line in original_gcf_file:
    split_line = line.split("\t")
    bgc_id = int(split_line[0])
    bgc_name = split_line[1]
    bgc_name = bgc_name[0:bgc_name.rfind(".")]
    if "genomic_2" in bgc_name or "_copy" in bgc_name:
        bgc_name1 = bgc_name[0:bgc_name.find("/")]
        bgc_name1 = bgc_name1[0:bgc_name1.rfind("_")]
        bgc_name2 = bgc_name[bgc_name.find("/"):len(bgc_name)]
        bgc_name = bgc_name1 + bgc_name2
    if ".gbff" in bgc_name:
        bgc_name1 = bgc_name[0:bgc_name.find("/")]
        bgc_name1 = bgc_name1[0:bgc_name1.rfind(".")]
        bgc_name2 = bgc_name[bgc_name.find("/"):len(bgc_name)]
        bgc_name = bgc_name1 + bgc_name2
    original_gcf = split_line[4]
    original_gcfs[bgc_name] = original_gcf
    bgc_names[bgc_id] = bgc_name

for line in l2_gcf_file:
    if "gcf" in line:
        continue
    split_line = line.split("\t")
    bgc_id = int(split_line[0])
    bgc_name = bgc_names[bgc_id]
    l2_gcfs[bgc_name] = int(split_line[1])
    if int(split_line[1]) not in list_of_bgcs_in_l2gcf:
        list_of_bgcs_in_l2gcf[int(split_line[1])] = []
    list_of_bgcs_in_l2gcf[int(split_line[1])].append(bgc_name)
print("determine gcf activities")    
for bgc_name in l2_gcfs:
    gcf = l2_gcfs[bgc_name]
    if "converted" in bgc_name:
        characterized_gcfs[gcf] = True
        continue
    if bgc_name not in probabilities["antibacterial"]:
            print("missing activity " + bgc_name)
            continue
    if bgc_name not in bgc_to_genome_dict:
            print("not in genome dictionary")
            print(bgc_name)
            continue
    if gcf not in gcf_activities:
        gcf_activities[gcf] = {}
        for activity in activities:
            gcf_activities[gcf][activity] = []
        
    if gcf not in characterized_gcfs:
        characterized_gcfs[gcf] = False
        
    genome = bgc_to_genome_dict[bgc_name]
    if genome not in genomes_in_bigslice_results:
        genomes_in_bigslice_results.append(genome)
    for activity in activities:
         gcf_activities[gcf][activity].append(probabilities[activity][bgc_name])

print("count active BGCs and GCFs")
total_active_counts = {}
uncharacterized_active_counts = {}
characterized_active_counts = {}

genus_in_dataset_counts = {}
for genus in genus_members:
    total_active_counts[genus] = {}
    uncharacterized_active_counts[genus] = {}
    characterized_active_counts[genus] = {}
    genus_in_dataset_counts[genus] = 0
    for activity in activities:
         total_active_counts[genus][activity] = []
         uncharacterized_active_counts[genus][activity] = []
         characterized_active_counts[genus][activity] = []
         for genome in genus_members[genus]:
             total_counts = -1
             if genome not in genome_bgcs:
                 if genome in genomes_without_bgcs:
                     total_counts = 0
                     characterized_counts = 0
                     uncharacterized_counts = 0
             else:                  
                total_counts, characterized_counts, uncharacterized_counts = getClusterCounts(genome_bgcs[genome], probabilities, activity, characterized_gcfs, l2_gcfs, 0.5)
             if total_counts == -1:
                 continue
             genus_in_dataset_counts[genus] += 1 
             total_active_counts[genus][activity].append(total_counts)
             uncharacterized_active_counts[genus][activity].append(uncharacterized_counts)
             characterized_active_counts[genus][activity].append(characterized_counts)
             
#sort genera by uncharacterized active antibiotic counts
average_uncharacterized_counts = {}
average_characterized_counts = {}
average_total_counts = {}
total_unique_unchar_gcf_by_genus = {}
list_of_unique_unchar_gcf_by_genus = {}
total_unique_char_gcf_by_genus = {}
list_of_unique_char_gcf_by_genus = {}

gcf_activities_avg = {}
active_gcfs = {}
gcf_presence_absence = {}
#classify GCFs based on activity
for activity in activities:
    gcf_activities_avg[activity] = {}
    active_gcfs[activity] = []
    gcf_presence_absence[activity] = {}
    for gcf in gcf_activities:
        avg = np.average(gcf_activities[gcf][activity])
        gcf_activities_avg[activity][gcf] = avg
        if avg > 0.5:
            active_gcfs[activity].append(gcf)
            
print(gcf_activities_avg)

for activity in activities:
    total_unique_unchar_gcf_by_genus[activity] = {}
    total_unique_char_gcf_by_genus[activity] = {}
    list_of_unique_unchar_gcf_by_genus[activity] = {}
    list_of_unique_char_gcf_by_genus[activity] = {}
    for genus in genus_members:
        total_unique_unchar_gcf_by_genus[activity][genus] = 0
        list_of_unique_unchar_gcf_by_genus[activity][genus] = []
        list_of_unique_char_gcf_by_genus[activity][genus] = []
        total_unique_char_gcf_by_genus[activity][genus] = 0
        for genome in genus_members[genus]:
            if genome not in genomes_in_bigslice_results or genome in genomes_without_bgcs:
                continue
            for bgc in genome_bgcs[genome]:
                if bgc not in l2_gcfs:
                    print("missing bgc from bigslice: " + bgc)
                    continue
                gcf = l2_gcfs[bgc]
                if bgc not in probabilities[activity]:
                    continue
                if not characterized_gcfs[gcf] and gcf not in list_of_unique_unchar_gcf_by_genus[activity][genus] and  gcf in active_gcfs[activity]:
                    list_of_unique_unchar_gcf_by_genus[activity][genus].append(gcf)
                if  characterized_gcfs[gcf] and gcf not in list_of_unique_char_gcf_by_genus[activity][genus] and  gcf in active_gcfs[activity]:
                    list_of_unique_char_gcf_by_genus[activity][genus].append(gcf)
        total_unique_unchar_gcf_by_genus[activity][genus] = len(list_of_unique_unchar_gcf_by_genus[activity][genus])
        total_unique_char_gcf_by_genus[activity][genus] = len(list_of_unique_char_gcf_by_genus[activity][genus])

for activity in activities:
    average_uncharacterized_counts[activity] = {}
    average_characterized_counts[activity] = {}
    average_total_counts[activity] = {}
    for genus in genus_members:
        unchar_counts = uncharacterized_active_counts[genus][activity]
        char_counts = characterized_active_counts[genus][activity]
        total_counts = total_active_counts[genus][activity]
        average_uncharacterized_counts[activity][genus] = 0
        average_characterized_counts[activity][genus] = 0
        average_total_counts[activity][genus] = 0
        if len(unchar_counts) != 0:
            avg_count = np.average(np.array(unchar_counts))
            average_uncharacterized_counts[activity][genus] = avg_count
        if len(char_counts) != 0:
            avg_count = np.average(np.array(char_counts))
            average_characterized_counts[activity][genus] = avg_count
        if len(total_counts)  != 0:
            avg_count = np.average(np.array(total_counts))
            average_total_counts[activity][genus] = avg_count
            
        
        
print("number of genomes: " + str(len(genomes_in_bigslice_results)))

sorted_avg_counts = sorted(average_uncharacterized_counts[activity_to_plot].items(), key=lambda x:x[1],reverse=True)
sorted_total_counts = sorted(total_unique_unchar_gcf_by_genus[activity_to_plot].items(), key=lambda x:x[1],reverse=True)

labels = []
labels_total = []
uncharacterized_data  = {}
characterized_data = {}
uncharacterized_total_data = {}
characterized_total_data = {}
for i in range(0,25):
    labels.append(sorted_avg_counts[i][0])
    labels_total.append(sorted_total_counts[i][0])
    for activity in activities:
        if activity not in uncharacterized_data:
            uncharacterized_data[activity] = []
            characterized_data[activity] = []
            uncharacterized_total_data[activity] = []
            characterized_total_data[activity] = []
        uncharacterized_data[activity].append(average_uncharacterized_counts[activity][sorted_avg_counts[i][0]])
        characterized_data[activity].append(average_characterized_counts[activity][sorted_avg_counts[i][0]])
        uncharacterized_total_data[activity].append(total_unique_unchar_gcf_by_genus[activity][sorted_total_counts[i][0]])
        characterized_total_data[activity].append(total_unique_char_gcf_by_genus[activity][sorted_total_counts[i][0]])

#make bar graphs
# Define the colors for the bars
colors = ['r', 'g', 'b', 'c', 'm', 'y']

fig, axs = plt.subplots(2, figsize=(15,10))

width = 0.25  # Set the width of each bar
set_width = width * 6 + 0.2  # Set the width of each set of bars, including spacing
x_positions = np.arange(len(labels)) * set_width  # Compute the x positions of the sets of bars
for i in range(6):
    x = x_positions + i * width  # Compute the x positions of the bars
    axs[0].bar(x, uncharacterized_data[activities[i]], width=width, color=colors[i], label=activities[i])
    axs[1].bar(x, characterized_data[activities[i]], width=width, color=colors[i], label=activities[i])
axs[0].legend() 
xtick_locs = x_positions + 2.5 * width
axs[1].set_xticks(xtick_locs)
axs[0].margins(x=0.01, tight=True)
axs[1].margins(x=0.01, tight=True)
axs[0].set_xticks(xtick_locs)
axs[0].xaxis.set_tick_params(labelbottom=False)
axs[1].set_xticklabels(labels, rotation=65, ha="right")
axs[0].set_title('Clusters in GCF with no MIBiG Clusters')
axs[1].set_title('Clusters in GCF with MIBiG Cluster(s)')
axs[1].set_xlabel('genus')
axs[0].set_ylabel('# of BGCs per Genome')   
axs[1].set_ylabel('# of BGCs per Genome')   
fig.suptitle("Number of BGCs per Genome")
fig.savefig(input_dir + activity_to_plot + '_per_genome_active_cluster_counts.pdf',bbox_inches = "tight")



            
#make plot based on total active gcfs per genome
#find top ten genera for each activity, make active gcf presence absence tables
fig3, axs3 = plt.subplots(2, figsize=(15,10))
for i in range(6):
    x = x_positions + i * width  # Compute the x positions of the bars
    axs3[0].bar(x, uncharacterized_total_data[activities[i]], width=width, color=colors[i], label=activities[i])
    axs3[1].bar(x, characterized_total_data[activities[i]], width=width, color=colors[i], label=activities[i])
axs3[0].legend() 
xtick_locs = x_positions + 2.5 * width
axs3[1].set_xticks(xtick_locs)
axs3[0].margins(x=0.01, tight=True)
axs3[1].margins(x=0.01, tight=True)
axs3[0].set_xticks(xtick_locs)
axs3[0].xaxis.set_tick_params(labelbottom=False)
axs3[1].set_xticklabels(labels_total, rotation=65, ha="right")
axs3[0].set_title('GCFs with no MIBiG Clusters')
axs3[1].set_title('GCFs with MIBiG Cluster(s)')
axs3[1].set_xlabel('genus')
axs3[0].set_ylabel('# of GCFs per Genera')   
axs3[1].set_ylabel('# of GCFs per Genera')   
fig3.suptitle("Number of GCF per Genera")
fig3.savefig(input_dir + activity_to_plot + '_per_genera_total_active_gcf_counts.pdf',bbox_inches = "tight")


#output sorted lists to text file
output_file = open(input_dir + activity_to_plot + '_per_genome_active_cluster_counts.txt','w')
output_file.write("genus,number of BGCs in cluster without MIBiG,number of BGCs in cluster with MIBiG,number of genomes analyzed\n")
for i in range(0,len(sorted_avg_counts)):
    output_file.write(sorted_avg_counts[i][0] + "," + str(average_uncharacterized_counts[activity_to_plot][sorted_avg_counts[i][0]]) + "," + str(average_characterized_counts[activity_to_plot][sorted_avg_counts[i][0]]) + "," + str(genus_in_dataset_counts[sorted_avg_counts[i][0]]) + "\n")
output_file.close()

output_file = open(input_dir + activity_to_plot + '_per_genera_active_gcf_counts.txt','w')
output_file.write("genus,number of GCFs without MIBiG,number of GCFs with MIBiG,number of genomes analyzed\n")
for i in range(0, len(sorted_total_counts)):
    output_file.write(sorted_total_counts[i][0] + "," + str(total_unique_unchar_gcf_by_genus[activity_to_plot][sorted_total_counts[i][0]]) + "," + str(total_unique_char_gcf_by_genus[activity_to_plot][sorted_total_counts[i][0]]) + "," + str(genus_in_dataset_counts[sorted_total_counts[i][0]]) + "\n")
output_file.close()
s



    