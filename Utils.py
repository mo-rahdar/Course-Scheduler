#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import pandas as pd
import math
import re, os
from collections import defaultdict, OrderedDict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import colors
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
from IPython.display import FileLink, Image, display, HTML
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)


# In[ ]:

def what_semester(input_file_path):

    sheetnames = dict()
    
    df = pd.read_excel(input_file_path, sheet_name='Info')
    semester = df.loc[:, 'Semester'].dropna().item()

    if semester =='Fall' or semester == 'fall':
        sheetnames['OurCourses'] = 'Fall Courses'
        sheetnames['MathChem'] = 'Fall Math&Chem'
        sheetnames['Conflicting'] = 'Fall Conflicting'
        sheetnames['Labs'] = 'Fall Labs'
        sheetnames['NonTrad'] = 'Fall Non-Trad'
        sheetnames['DesiredTimes'] = 'Fall Desired Times'
        sheetnames['UndesiredTimes'] = 'Fall Undesired Times'
    elif semester == 'Spring' or semester == 'spring':
        sheetnames['OurCourses'] = 'Spring Courses'
        sheetnames['MathChem'] = 'Spring Math&Chem'
        sheetnames['Conflicting'] = 'Spring Conflicting'
        sheetnames['Labs'] = 'Spring Labs'
        sheetnames['NonTrad'] = 'Spring Non-Trad'
        sheetnames['DesiredTimes'] = 'Spring Desired Times'
        sheetnames['UndesiredTimes'] = 'Spring Undesired Times'
    else:
        return print('Check the data you entered for the semester in the "Info" worksheet')

    return sheetnames


# In[1]:


def read_data(input_file_path, sheetnames):

    for key, value in sheetnames.items():
        globals()[key] = value

    data = dict()

    # Dictionary of time slot info    
    time_slot_info = {}
    
    # Read the Excel file into a DataFrame
    df = pd.read_excel(input_file_path, sheet_name='Info')
    
    # Extract data from the first column and remove NaN values
    times = df.loc[:, 'Time slots'].dropna().tolist()
    
    for slot in times:
        
        # Extract day, start time, and duration from the slot
        d, s, u = re.search(r'([A-Za-z]+)(\d+)-(\d+)', slot).groups()
        
        # Convert start time to float
        s = float(s)
        
        # Adjust start and end time if the start time is 15 or 18 (it is 15:30 or 18:30)
        if s == 15 or s == 18:
            s += 0.5
            
        # Calculate end time
        e = s + int(u) / 60
        
        # Create an inner dictionary for the time slot
        info = {'M': 1 if 'M' in d else 0,
                'T': 1 if 'T' in d else 0,
                'W': 1 if 'W' in d else 0,
                'R': 1 if 'R' in d else 0,
                'F': 1 if 'F' in d else 0,
                'duration': u,
                'start': s,
                'end': e
               }
        
        # Add the inner dictionary to the outer dictionary
        time_slot_info[slot] = info
    data['time_slot_info'] = time_slot_info

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of courses info (engr, phys, math, chem)
    df1 = pd.read_excel(input_file_path, sheet_name=OurCourses)
    df2 = pd.read_excel(input_file_path, sheet_name=MathChem)
    
    # Convert the DataFrame to a nested dictionary
    courses_info_dict = {}
    
    for index, row in df1.iterrows():
        key = row['Item']
        inner_dict = {}
        for column in df1.columns[1:8]:
            inner_dict[column] = row[column]
        courses_info_dict[key] = inner_dict
        
    for index, row in df2.iterrows():
        key = row['Item']
        inner_dict = {}
        for column in df2.columns[1:5]:
            inner_dict[column] = row[column]
        inner_dict['Duration'] = f"{row['Course']}-u"
        courses_info_dict[key] = inner_dict
            
    data['courses_info_dict'] = courses_info_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Department meetings durations
    df = pd.read_excel(input_file_path, sheet_name='Info')
    
    # Extract data and remove NaN values
    meetings_duration = df.loc[:, 'Meetings'].dropna().item()
    data['meetings_duration'] = meetings_duration

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of meeting times -> Set M
    M_list = []
    
    # Days of the week with their respective end times
    days_end_hours = {'M': 17, 'T': 17, 'W': 17, 'R': 17, 'F': 14}
    
    # Generate the times
    for d, e in days_end_hours.items():
        h = 8.0  # Start time
        while h + meetings_duration <= e:
            M_list.append(f"{d}{h}")
            h += meetings_duration
    data['M_list'] = M_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of meeting times info
    meeting_info = {}
    
    for m in M_list:
        # Extract day and start time
        match = re.match(r'([A-Z]+)(\d+(\.\d+)?)$', m)
        if match:
            d, start_str, _ = match.groups()
    
            # Convert start time to float
            s = float(start_str)
    
            # Calculate end time (start time + meeting duration - 1 minute)
            e = s + meetings_duration - 1 / 60
    
            # Create an inner dictionary for the time slot
            info = {
                'M': 1 if 'M' in d else 0,
                'T': 1 if 'T' in d else 0,
                'W': 1 if 'W' in d else 0,
                'R': 1 if 'R' in d else 0,
                'F': 1 if 'F' in d else 0,
                'duration': meetings_duration * 60,  # Duration in minutes
                'start': s,
                'end': e
            }
    
            # Add to the dictionary
            meeting_info[m] = info
    data['meeting_info'] = meeting_info

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of all time slots -> set T
    
    # We have already read time slots from the Excel file and stored it in times
    T_list = times
    data['T_list'] = T_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of 1-hr time blocks -> set B
    B_list = []
    
    # Days of the week with their respective end times
    days_end_hours = {'M': 18, 'T': 18, 'W': 18, 'R': 18, 'F': 14}
    
    # Generate the times
    for d, e in days_end_hours.items():
        # Ensure the last block does not exceed the end time
        for h in range(8, e):
            if h + 1 > e:
                break
            B_list.append(f'{d}{h}')
    data['B_list'] = B_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of 1-hr block times info
    block_info = {}
    
    for m in B_list:
        
        # Extract the day and start time from the block times
        d, s = re.search(r'([A-Z]+)(\d+)', m).groups()
        
        # Convert start time to float
        s = float(s)
    
        # Calculate end time (start time + 59 minutes)
        e = s + 59/60
        
        # Create an inner dictionary for the time slot
        info = {'M': 1 if 'M' in d else 0,
                'T': 1 if 'T' in d else 0,
                'W': 1 if 'W' in d else 0,
                'R': 1 if 'R' in d else 0,
                'F': 1 if 'F' in d else 0,
                'duration': meetings_duration*60,
                'start': s,
                'end': e
               }
        
        # Add the inner dictionary to the outer dictionary
        block_info[m] = info
    data['block_info'] = block_info

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of courses -> set J
    J_list = []
    
    # Iterate over the values of the courses_info_dict dictionary
    for c in courses_info_dict.values():
        
        # Extract the desired value (e.g., 'Course' or 'Professor') and add it to the list
        J_list.append(c['Course'])
    
    # Remove duplicate values by converting the list to a set
    J_list = list(set(J_list))
    data['J_list'] = J_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of sections -> set S
    df = pd.read_excel(input_file_path, sheet_name='Info')
    
    # Extract data from the first column and remove NaN values
    S_list = df.loc[:, 'Sections'].dropna().tolist()
    data['S_list'] = S_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of weekdays -> set D
    df = pd.read_excel(input_file_path, sheet_name='Info')
    
    # Extract data from the first column and remove NaN values
    D_list = df.loc[:, 'Days'].dropna().tolist()
    data['D_list'] = D_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of class durations (set U) in this format: 'duration - number of days in a week'
    df1 = pd.read_excel(input_file_path, sheet_name='Info')
    df2 = pd.read_excel(input_file_path, sheet_name=MathChem)
    
    # Extract data from the first column and remove NaN values
    U_ourcourses = df1['Durations'].dropna().tolist()
    U_mathchem = df2['Course'].unique().tolist()
    U_mathchem = [f"{item}-u" for item in U_mathchem]
    U_list = U_ourcourses + U_mathchem
    
    data['U_list'] = U_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of professors -> set P
    P_list = []
    
    # Iterate over the values of the courses_info_dict dictionary
    for p in courses_info_dict.values():
        
        # Check if the 'Professor' key exists in the inner dictionary
        if 'Professor1' in p:
            # Extract the professor's name and add it to the list
            P_list.append(p['Professor1'])
            
    # Remove duplicate values by converting the list to a set
    P_list = sorted(list(set(P_list)))
    data['P_list'] = P_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of preferences for what time a day to teach -> set R
    df = pd.read_excel(input_file_path, sheet_name='Info')
    
    # Extract data from the first column and remove NaN values
    R_list = df.loc[:, 'Preferences'].dropna().tolist()
    data['R_list'] = R_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of courses and sections -> set JS
    JS_list = []
    
    # Iterate over the items of courses_info_dict
    for i, c_info in courses_info_dict.items():
        c = c_info['Course']
        s = c_info['Section']
        
        # Append the course and section as a tuple to the list
        JS_list.append((c, s))
    data['JS_list'] = JS_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create a dictionary where course names are keys and sections are values -> set JS
    JS_dict = {}
    
    for k, d in courses_info_dict.items():
        c = d['Course']
        s = d['Section']
        
        # If the course is not already in the dictionary, initialize it as an empty list
        if c not in JS_dict:
            JS_dict[c] = []
        
        # Append the section to the course
        JS_dict[c].append(s)
    data['JS_dict'] = JS_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of courses, sections, and durations for our courses -> set JSU
    JSU_dict = {}
    
    # Iterate over the items of courses_info_dict
    for u in courses_info_dict.values():
        
        # Check if the 'Duration' key exists in the inner dictionary
        if 'Duration' in u:
            c = u['Course']
            s = u['Section']
            d = u['Duration']
    
            # Append the course and section as a tuple to the list
            JSU_dict[(c, s)] = [d]
    data['JSU_dict'] = JSU_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary to store union of time slots for each course j -> set JU(j)
    JU_dict = defaultdict(set)
    
    # Iterate over the original dictionary
    for (c, s), t in JSU_dict.items():
        JU_dict[c].update(t)
    
    # Convert sets to lists
    JU_dict = {c: list(t) for c, t in JU_dict.items()}
    data['JU_dict'] = JU_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of professors who don't want to teach in a specific day -> set PD
    df = pd.read_excel(input_file_path, sheet_name='Undesired Days')
    
    # Filter the DataFrame to include only rows where 'Not desirable day' is not empty
    filtered_df = df.dropna(subset=['Not desirable day'])
    
    # Create a professor day list from the DataFrame
    PD_list = list(filtered_df[['Professors',
                                'Not desirable day']].itertuples(index=False,
                                                                 name=None))
    data['PD_list'] = PD_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of courses that desired to be scheduled at a more desirable time -> set JST1
    df = pd.read_excel(input_file_path, sheet_name=DesiredTimes)
    
    # Create a preferred time list from the DataFrame
    JST1_list = list(df[['Course',
                         'Section',
                         'Time']].itertuples(index=False, name=None))
    data['JST1_list'] = JST1_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of courses that are not desired to be scheduled at a certain time -> set JST2
    df = pd.read_excel(input_file_path, sheet_name=UndesiredTimes)
    
    # Create an avoid time list from the DataFrame
    JST2_list = list(df[['Course',
                         'Section',
                         'Time']].itertuples(index=False, name=None))
    data['JST2_list'] = JST2_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of Math and Chem courses with fixed schedule -> set JST3
    
    # Read the Excel file
    df = pd.read_excel(input_file_path, sheet_name=MathChem)
    
    # Create a fixed time list from the DataFrame
    JST3_list = list(df[['Course',
                         'Section',
                         'Time']].itertuples(index=False, name=None))
    data['JST3_list'] = JST3_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of professors and their courses -> set PJ(p)
    PJ_dict = {}
    
    # Iterate over the items of courses_info_dict
    for item in courses_info_dict.values():
        # Check if the 'Professor' key exists in the inner dictionary
        if 'Professor1' in item:
            p1 = item['Professor1']
            p2 = item['Professor2']
            c = item['Course']
            s = item['Section']
            
            # Add Professor1 and their courses to PJ_dict
            if p1:  # Ensure p1 is not None or empty
                if p1 in PJ_dict:
                    PJ_dict[p1].append((c, s))
                else:
                    PJ_dict[p1] = [(c, s)]
    
            # Add Professor2 and their courses to PJ_dict, skipping if it's nan
            if p2 and not (isinstance(p2, float) and math.isnan(p2)):  # Skip if p2 is nan
                if p2 in PJ_dict:
                    PJ_dict[p2].append((c, s))
                else:
                    PJ_dict[p2] = [(c, s)]
    data['PJ_dict'] = PJ_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of time slots for duration u -> set UT(u)    
    UT_dict = {'50min-1':    [i for i in times if '50' in i and '150' not in i and i.find(next(filter(str.isdigit, i))) == 1],
               '50min-2':    [i for i in times if '50' in i and i.find(next(filter(str.isdigit, i))) == 2],
               '50min-3':    [i for i in times if '50' in i and i.find(next(filter(str.isdigit, i))) == 3],
               '75min-1':    [i for i in times if '75' in i and i.find(next(filter(str.isdigit, i))) == 1],
               '75min-2':    [i for i in times if '75' in i and i.find(next(filter(str.isdigit, i))) == 2],
               '75min-3':    [i for i in times if '75' in i and i.find(next(filter(str.isdigit, i))) == 3],
               '100min-1':   [i for i in times if '100' in i and i.find(next(filter(str.isdigit, i))) == 1],
               '110min-1':   [i for i in times if '110' in i and i.find(next(filter(str.isdigit, i))) == 1],
               '150min-1':   [i for i in times if '150' in i and i.find(next(filter(str.isdigit, i))) == 1],
               '170min-1':   [i for i in times if '170' in i and i.find(next(filter(str.isdigit, i))) == 1],
              }
    
    # Add math and chem durations to the dictionary
    df = pd.read_excel(input_file_path, sheet_name=MathChem)
    # Group by course name and collect unique time slots
    MC_times = df.groupby('Course')['Time'].unique().apply(list).to_dict()
    # Update keys to include '-u'
    MC_times = {f"{key}-u": value for key, value in MC_times.items()}

    UT_dict.update(MC_times)

    # MC_times = list(df[['Duration',
    #                     'Time']].itertuples(index=False, name=None))
    
    # for k, v in MC_times:
    #     UT_dict[k] = [v]
    data['UT_dict'] = UT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of all available time slots for course section (j,s) -> set JST(j,s)    
    JST_dict = {k: UT_dict[v[0]] for k, v in JSU_dict.items()}
    
    # # Add math and chem durations to the dictionary
    # df = pd.read_excel(input_file_path, sheet_name=MathChem)
    # # Group by course name and collect unique time slots
    # MC_times = df.groupby('Course')['Time'].unique().apply(list).to_dict()
    # # Update keys to include '-u'
    # MC_times = {f"{key}-u": value for key, value in MC_times.items()}
    
    # Adding math and chem courses to JST_dict
    # for (j, s, t) in JST3_list:
    #     JST_dict[(j, s)] = [t]
    data['JST_dict'] = JST_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary to store union of time slots for each course j -> set JT(j)
    JT_dict = defaultdict(set)
    
    # Iterate over the original dictionary
    for (c, s), t in JST_dict.items():
        JT_dict[c].update(t)
    
    # Convert sets to lists
    JT_dict = {c: list(t) for c, t in JT_dict.items()}
    data['JT_dict'] = JT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of professors and available time slots for their courses -> set PT(p)
    PT_dict = {}
    
    # Loop through each professor and their associated (j, s) pairs
    for p, js_pairs in PJ_dict.items():
        # Initialize an empty list to store unique time slots for professor p
        ts = []
    
        # Loop through each (j, s) pair for this professor
        for js in js_pairs:
            # Add all time slots from JST corresponding to this (j, s) pair
            if js in JST_dict:
                for t in JST_dict[js]:
                    if t not in ts:  # Avoid duplicates while preserving order
                        ts.append(t)
    
        # Assign the list of time slots to PT[p] without sorting
        PT_dict[p] = ts
    data['PT_dict'] = PT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of time slots for day d -> set DT(d)
    DT_dict = {}
    
    # Define a mapping of day initials to full day names
    day_mapping = {'M': 'Mon', 'T': 'Tue', 'W': 'Wed', 'R': 'Thu', 'F': 'Fri'}
    
    # Iterate over the list of times
    for t in T_list:
        
        # Extract the day initials using a regular expression
        day_initials = re.search(r'[A-Za-z]+', t).group()
        
        # Iterate over the extracted day initials
        for i in day_initials:
            
            # Map the day initial to the full day name
            d = day_mapping.get(i)
            
            # Check if the day already exists in the dictionary
            if d in DT_dict:
                # Append the time to the list of times for the corresponding day
                DT_dict[d].append(t)
            else:
                # Create a new list with the time for the corresponding day
                DT_dict[d] = [t]
    data['DT_dict'] = DT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of preference r and their time slots belong to it -> set RT(r)
    
    # Extract the start time of classes from time slots
    start_times = [int(re.search(r'[A-Za-z]+(\d+)-', i).group(1))
                   for i in times if re.search(r'[A-Za-z]+(\d+)-', i)]
    
    RT_dict = {'Morning':         [i for (i, j) in zip(times, start_times) if             j <  12],
               'ExtendedMorning': [i for (i, j) in zip(times, start_times) if             j <  14],
               'WholeDay':        [i for (i, j) in zip(times, start_times) if             j <= 16],
               'Midday':          [i for (i, j) in zip(times, start_times) if j >= 10 and j <= 14],
               'ExtendedMidday':  [i for (i, j) in zip(times, start_times) if j >= 10 and j <= 16],
               'Afternoon':       [i for (i, j) in zip(times, start_times) if j >= 12 and j <= 16],
               'Evening':         [i for (i, j) in zip(times, start_times) if j >= 16            ]
              }
    data['RT_dict'] = RT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of professors and their preferred time to teach -> set PR(p)
    
    # Read the Excel file into a Pandas DataFrame
    df = pd.read_excel(input_file_path, sheet_name='Teaching Times')
    
    # Create a dictionary with Professors as keys and Preferred time to teach as values
    PR_dict = pd.Series(df['Preferred time to teach'].values,
                        index=df['Professors']).to_dict()
    PR_dict = {p: [PR_dict[p]] for p in PR_dict}
    data['PR_dict'] = PR_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of levels with incompatible courses -> set F
    df = pd.read_excel(input_file_path, sheet_name=Conflicting)
    
    # Extract the 'Group' column, drop any NaN values, and find unique values
    F_list = df['Group'].dropna().unique()
    data['F_list'] = F_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of incompatible courses -> set FJ(f)
    df = pd.read_excel(input_file_path, sheet_name=Conflicting)
    
    # Initialize a dictionary to store groups and their courses
    FJ_dict = defaultdict(list)
    
    # Track the current group while iterating through the rows
    current_group = None
    
    for index, row in df.iterrows():
        group = row['Group']
        course = row['Courses']
        
        # If we encounter a new group, update the current group
        if pd.notna(group):
            current_group = group
        
        # If there's a valid group and course, add the course to the group
        if current_group and pd.notna(course):
            FJ_dict[current_group].append(course)
    
    # Convert defaultdict to a regular dictionary
    FJ_dict = dict(FJ_dict)
    data['FJ_dict'] = FJ_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of the number of courses in each group of incompatible courses -> set FN(f)
    FN_dict = {g: len(cs) for g, cs in FJ_dict.items()}
    data['FN_dict'] = FN_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of conflicted time slots -> set GT(g)
    GT_dict = {}
    
    for slot1, info1 in time_slot_info.items():
        
        for slot2, info2 in time_slot_info.items():
            
            # Skip if comparing the same slot or if they are on different days
            if ((info1['M'] == 1 and info2['M'] == 1) or
                (info1['T'] == 1 and info2['T'] == 1) or
                (info1['W'] == 1 and info2['W'] == 1) or
                (info1['R'] == 1 and info2['R'] == 1) or
                (info1['F'] == 1 and info2['F'] == 1)):
             
                # Check for conflicts based on start and end times
                if (
                    (info1['start'] <= info2['start'] < info1['end']) or
                    (info1['start'] < info2['end'] <= info1['end']) or
                    (info2['start'] <= info1['start'] < info2['end']) or
                    (info2['start'] < info1['end'] <= info2['end'])       ):
                    
                    # Add conflicted slots to the conflicts dictionary
                    if slot1 in GT_dict:
                        GT_dict[slot1].append(slot2)
                    else:
                        GT_dict[slot1] = [slot2]
    data['GT_dict'] = GT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of times conflicted with other time slots -> set G
    G_list = []
    
    # Iterate over the items of GT_dict
    for item in GT_dict.keys():
        
        G_list.append(item)
    data['G_list'] = G_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load Excel data into a Pandas DataFrame
    df = pd.read_excel(input_file_path, sheet_name='Busy Times')
    
    # Mapping full day names to shorthand format
    day_mapping = {'Mon': 'M', 'Tue': 'T', 'Wed': 'W', 'Thu': 'R', 'Fri': 'F'}
    
    # Parse restrictions into a dictionary
    restrictions = {}
    for _, row in df.iterrows():
        p = row['Professor']
        d = row['Busy day']
        s = row['From']
        e = row['To']
        
        if pd.notna(p) and pd.notna(d) and pd.notna(s) and pd.notna(e):        
            day_short = day_mapping.get(d, None) # Convert full day name to shorthand
            if day_short:
                if p not in restrictions:
                    restrictions[p] = {}
                restrictions[p][day_short] = (s, e)
    data['restrictions'] = restrictions

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of courses that cannot be at certain time slots (busy times for a professor) -> set JST4
    JST4_list = []
    
    for p, ds in restrictions.items():
        if p in PJ_dict:  # Skip professors without courses
            for j, s in PJ_dict[p]:
                if (j, s) in JST_dict:
                    for t in JST_dict[j, s]:  # Times available for (j, s)
                        time_data = time_slot_info[t]
                        start_time = time_slot_info[t]['start']
                        end_time = time_slot_info[t]['end']
    
                        # Check if time slot t conflicts with restrictions
                        for d, is_active in time_data.items():
                            if is_active == 1 and d in ds:  # Only consider active days and restricted days
                                restricted_start, restricted_end = ds[d]
                                
                                # Check for overlap
                                if (restricted_start <= start_time < restricted_end or
                                    restricted_start <= end_time < restricted_end):
                                    JST4_list.append((j, s, t))
    data['JST4_list'] = JST4_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # List of non-traditional courses that we may not want to consider them in the working span -> set JST5
    df = pd.read_excel(input_file_path, sheet_name=NonTrad)
    
    # Create a preferred time list from the DataFrame
    JST5_list = list(df[['Course', 'Section']].itertuples(index=False, name=None))
    data['JST5_list'] = JST5_list

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of conflicted time slots with block b -> set BT(b)
    BT_dict = {}
    
    for slot1, info1 in block_info.items():
        
        for slot2, info2 in time_slot_info.items():
            
            # Skip if comparing the same slot or if they are on different days
            if ((info1['M'] == 1 and info2['M'] == 1) or
                (info1['T'] == 1 and info2['T'] == 1) or
                (info1['W'] == 1 and info2['W'] == 1) or
                (info1['R'] == 1 and info2['R'] == 1) or
                (info1['F'] == 1 and info2['F'] == 1)):
             
                # Check for conflicts based on start and end times
                if (
                    (info1['start'] <= info2['start'] < info1['end']) or
                    (info1['start'] < info2['end'] <= info1['end']) or
                    (info2['start'] <= info1['start'] < info2['end']) or
                    (info2['start'] < info1['end'] <= info2['end'])       ):
                    
                    # Add conflicted slots to the conflicts dictionary
                    if slot1 in BT_dict:
                        BT_dict[slot1].append(slot2)
                    else:
                        BT_dict[slot1] = [slot2]
    data['BT_dict'] = BT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of conflicted time slots with block b -> set MT(m)
    MT_dict = {}
    
    for slot1, info1 in meeting_info.items():
        
        for slot2, info2 in time_slot_info.items():
            
            # Skip if comparing the same slot or if they are on different days
            if ((info1['M'] == 1 and info2['M'] == 1) or
                (info1['T'] == 1 and info2['T'] == 1) or
                (info1['W'] == 1 and info2['W'] == 1) or
                (info1['R'] == 1 and info2['R'] == 1) or
                (info1['F'] == 1 and info2['F'] == 1)):
             
                # Check for conflicts based on start and end times
                if (
                    (info1['start'] <= info2['start'] < info1['end']) or
                    (info1['start'] < info2['end'] <= info1['end']) or
                    (info2['start'] <= info1['start'] < info2['end']) or
                    (info2['start'] < info1['end'] <= info2['end'])       ):
                    
                    # Add conflicted slots to the conflicts dictionary
                    if slot1 in MT_dict:
                        MT_dict[slot1].append(slot2)
                    else:
                        MT_dict[slot1] = [slot2]
    data['MT_dict'] = MT_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionary of the courses or labs that need to be scheduled back to back on a same day
    df = pd.read_excel(input_file_path, sheet_name=Labs)
    
    # Initialize course_data dictionary
    LS_dict = {}
    
    # Iterate over rows of the DataFrame
    for _, row in df.iterrows():
        c = row['Course']
        # Drop any NaN or blank values in the section's columns
        s = [sec for sec in row[1:].values if pd.notna(sec)]
        # Add course and sections to the dictionary
        LS_dict[c] = {'sections': s, 'offsets': list(range(len(s)))}
    data['LS_dict'] = LS_dict

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Dictionaries of max hours a day and max days a week
    df = pd.read_excel(input_file_path, sheet_name='Teaching Load')
    
    hrs_a_day = {}
    
    # Build the dictionary from the DataFrame
    max_hrs_dict = df.set_index('Professors')['Max hrs a day'].to_dict()
    max_days_dict = df.set_index('Professors')['Max days a week'].to_dict()
    
    # If there is no data in the dictionary, replace it with a default value
    max_hrs_dict = {k: 10. if np.isnan(v) else v for k, v in max_hrs_dict.items()}
    max_days_dict = {k: 5. if np.isnan(v) else v for k, v in max_days_dict.items()}

    data['max_hrs_dict'] = max_hrs_dict
    data['max_days_dict'] = max_days_dict

    return (data)


# In[ ]:

def build_model(data):

    for key, value in data.items():
        globals()[key] = value
    
    model = pyo.ConcreteModel() # create a concrete model using Pyomo
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The indices for constraints and variables (Sets)
    model.T = pyo.Set(initialize=T_list)  # set of time slots (T)
    model.B = pyo.Set(initialize=B_list)  # set of block times (B)
    model.J = pyo.Set(initialize=J_list)  # set of courses (J)
    model.S = pyo.Set(initialize=S_list)  # set of sections (S)
    model.D = pyo.Set(initialize=D_list)  # set of weekdays (D)
    model.U = pyo.Set(initialize=U_list)  # set of class durations (U)
    model.P = pyo.Set(initialize=P_list)  # set of professors (P)
    model.R = pyo.Set(initialize=R_list)  # set of preference (R)
    model.G = pyo.Set(initialize=G_list)  # set of conflicted time slots (G)
    model.F = pyo.Set(initialize=F_list)  # set of incompatible course groups (F)
    model.L = pyo.Set(initialize=LS_dict.keys()) # set of labs whose sections need to be back to back (L)
    
    model.JST1 = pyo.Set(initialize=JST1_list) # set of courses wanted to be scheduled at certain times
    model.JST2 = pyo.Set(initialize=JST2_list) # set of courses wanted not to be scheduled at certain times
    model.JST3 = pyo.Set(initialize=JST3_list) # set of courses have to be scheduled at certain times
    model.JST4 = pyo.Set(initialize=JST4_list) # set of courses cannot be scheduled at certain times (busy times)
    model.JST5 = pyo.Set(initialize=JST5_list) # set of non traditional courses
    model.PD   = pyo.Set(initialize=PD_list)   # set of professors and days they don't want to teach
    model.JS1  = pyo.Set(initialize=JS_list)   # set of courses and sections (tuple)

    model.JS = pyo.Set(model.J, within=model.S, initialize=JS_dict)   # set of courses and sections (dic)
    model.LS = pyo.Set(model.L, within=model.S, initialize={l: data['sections'] for l, data in LS_dict.items()}) # set of sections for back to back courses
    model.JSU= pyo.Set(model.JS1, within=model.U, initialize=JSU_dict) # set of our courses, sections, and durations
    model.JU = pyo.Set(model.J, within=model.U, initialize=JU_dict)    # set of all courses and possible durations
    model.JST= pyo.Set(model.JS1, within=model.T, initialize=JST_dict) # set of all courses, sections, and possible times
    model.JT = pyo.Set(model.J, within=model.T, initialize=JT_dict)    # set of all courses and possible times
    model.PJ = pyo.Set(model.P, within=model.JS1, initialize=PJ_dict)   # set of courses for professor p
    model.PT = pyo.Set(model.P, within=model.T, initialize=PT_dict)    # set of time slots for professor p
    model.PR = pyo.Set(model.P, within=model.R, initialize=PR_dict)    # set of professors and their preferred time to teach
    model.UT = pyo.Set(model.U, within=model.T, initialize=UT_dict)    # set of time slots for duration u
    model.DT = pyo.Set(model.D, within=model.T, initialize=DT_dict)    # set of time slots in day d
    model.RT = pyo.Set(model.R, within=model.T, initialize=RT_dict)    # set of time slots for preference r
    model.GT = pyo.Set(model.G, within=model.T, initialize=GT_dict)    # set of times conflicted with a time slot
    model.FJ = pyo.Set(model.F, within=model.J, initialize=FJ_dict)    # set of conflicting courses in a semester
    model.BT = pyo.Set(model.B, within=model.T, initialize=BT_dict)    # set of blocks and time slots conflicts

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create the Pyomo parameter back to back courses offsets
    
    # Flatten the LS_dict to create the offset dictionary
    offset_data = {(c, s): data['offsets'][idx]
                   for c, data in LS_dict.items()
                   for idx, s in enumerate(data['sections'])}
    
    model.q = pyo.Param(offset_data.keys(),  # The keys are course-section pairs
                        initialize=offset_data,
                        within=pyo.NonNegativeIntegers)  # Assuming offsets are non-negative

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Define the parameters with [p, d] domain: h (max hrs a day), k (max days a week)
    def max_hrs_init(model, p, d):
        return max_hrs_dict[p]  # Use professor-specific max hours
    
    model.h = pyo.Param(model.P, model.D, initialize=max_hrs_init)
    model.k = pyo.Param(model.P, initialize=max_days_dict)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create an empty list to hold the variable x indices
    x_ind_list = []
    
    # Loop through each (j, s) pair in JS and corresponding timeslots in JST
    for (j, s) in JS_list:
        for t in model.JST[(j, s)]:
            # Initialize the variable for each valid (j, s, t) combination
            x_ind_list.append((j, s, t))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # A function to create variable y indices
    def valid_index_y(model):
        return ((f, j, t) for f in model.F for j in model.FJ[f] for t in model.JT[j])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Nested dictionary for time slots {l: {s: {d: [t]}}}
    time_slots_dict = defaultdict(lambda: defaultdict(dict))
    
    # Populate the dictionary
    for l in model.L:  # Courses
        for s in model.LS[l]:  # Sections for the course
            for d in model.D:  # Days
                # Filter time slots for the specific day and course-section pair
                day_slots = [t for t in model.DT[d] if t in model.JST[l, s]]
                time_slots_dict[l][s][d] = day_slots  # Assign to the dictionary
        
    time_slots_dict = dict(time_slots_dict)

    # A function to create variable z indices
    def valid_index_z(model):
        valid_indices = []
        for l in model.L:  # Iterate over courses
            for t in model.JT[l]:  # Iterate over time slots for course l
                day = next(d for d in model.D if t in model.DT[d])  # Get the day corresponding to t
                day_slots = time_slots_dict[l][list(time_slots_dict[l].keys())[0]][day]  # Use section A as a representative to check slots
                
                t_index = day_slots.index(t)  # Get index of time slot t
                max_offset = max(model.q[l, s] for s in model.LS[l])  # Maximum offset for course l
                
                # Ensure there is room for the offset
                if t_index + max_offset < len(day_slots):
                    valid_indices.append((l, t))
        
        return valid_indices

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # A function to create variable u indices
    valid_u_indices = {(p, d, t)
                       for p in model.P
                       for d in model.D
                       for t in model.DT[d] if t in model.PT[p]}

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Variables of the model
    
    # Binary variable x(j,s,t) is 1 if course j, section s, scheduled in time t
    model.x = pyo.Var(x_ind_list, domain=pyo.Binary)
    
    # Binary variable y(f,j,t) is 1 if course j has an available section at time t
    model.y = pyo.Var(valid_index_y, domain=pyo.Binary)
    
    # Binary variable z(l,t) is 1 if all sections of course l are scheduled back-to-back starting at t
    model.z = pyo.Var(valid_index_z, domain=pyo.Binary)
    
    # Binary variable tau(j,s,t) is 1 if course j, section s, not scheduled at desirable time t
    model.tau = pyo.Var(x_ind_list, domain=pyo.NonNegativeReals)
    
    # Binary variable pi(j,s,t) is 1 if course j, section s, scheduled at not desirable time t
    model.pi = pyo.Var(x_ind_list, domain=pyo.NonNegativeReals)
    
    # Binary variable gamma(j,s,t) is 1 if course j, section s, scheduled at time t,
    # which is not a desirable day for the course professor
    model.gamma = pyo.Var(x_ind_list, domain=pyo.NonNegativeReals)
    
    # Integer variable eta(f,b) indicates the number of conflicted courses in group f with block b
    model.eta = pyo.Var(model.F*model.B, domain=pyo.NonNegativeReals)
    
    # Integer variable mu(p) indicates the number of courses of Professor P
    # that do not satisfy their preferred time to teach
    model.mu = pyo.Var(model.P, domain=pyo.NonNegativeReals)
    
    # Binary variable alpha(p,d,t) indicates if t is the first class of professor p on day d
    model.alpha = pyo.Var(valid_u_indices, within=pyo.Binary)
    
    # Binary variable beta(p,d,t) indicates if t is the last class of professor p on day d
    model.beta = pyo.Var(valid_u_indices, within=pyo.Binary)
    
    # Continuous variable w(p,d) indicates the working span of professor p on day d
    model.w = pyo.Var(model.P, model.D, within=pyo.NonNegativeReals)
    
    # Continuous variable theta(p,d) indicates the extra working hours of professor p on day d
    model.theta = pyo.Var(model.P, model.D, within=pyo.NonNegativeReals)
    
    # Binary variable m(p,d) indicates if professor p teaches on day d
    model.m = pyo.Var(model.P, model.D, within=pyo.Binary)
    
    # Continuous variable delta(p,d) indicates the extra working days of professor p each week
    model.delta = pyo.Var(model.P, within=pyo.NonNegativeReals)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Define the objective function expression
    def obj_expression(model):
        return  (pyo.summation(model.tau)
                 + pyo.summation(model.pi)
                 + pyo.summation(model.gamma)
                 + pyo.summation(model.mu)
                 + pyo.summation(model.eta)
                 + pyo.summation(model.theta)
                 + pyo.summation(model.delta)
                )
    
    model.Obj = pyo.Objective(rule=obj_expression, sense=pyo.minimize)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule1(model, f, j, t):
        if f in model.F and j in model.FJ[f] and t in model.JT[j]:
            return model.y[f,j,t] <= sum(model.x[j,s,t]
                                         for s in model.JS[j]
                                         if (j, s) in model.JST and t in model.JST[j, s])
        else:
            return pyo.Constraint.Skip
    
    model.Const1 = pyo.Constraint(model.F, model.J, model.T, rule=constraint_rule1)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule2(model, f, j):
        if j in model.FJ[f]:
            return (sum(sum(model.y[f,j,t]
                            for t in model.UT[u])
                        for u in model.JU[j])
                    == 1)
        else:
            return pyo.Constraint.Skip
    
    model.Const2 = pyo.Constraint(model.F, model.J, rule=constraint_rule2)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule3(model, f, b):
        if f in model.F and b in model.B:
            # Collect valid terms
            valid_terms = [model.y[f, j, t] for j in model.FJ[f] for t in model.BT[b] if t in model.JT[j]]
            
            if valid_terms:  # Only create constraint if there are valid terms
                expr = sum(valid_terms)
                return expr <= 1 + model.eta[f,b]
            else:  # If no valid terms, return a feasible constraint
                return pyo.Constraint.Feasible
        else:
            return pyo.Constraint.Skip
    
    model.Const3 = pyo.Constraint(model.F, model.B, rule=constraint_rule3)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Courses taught by one professor cannot be assigned to conflicted time slots
    def constraint_rule4(model, p, g):
        if p in model.P:
            # Check if g is a valid key in model.GT
            if g not in model.GT:
                return pyo.Constraint.Skip
            
            if g in model.PT[p]:
                # Filter valid time slots (Use a set to avoid duplicate time slots)
                valid_t2s = set([t for (j, s) in model.PJ[p] if (j, s) in model.JST
                                 for t in model.GT[g] if t in model.JST[(j, s)]])
                
                # Skip constraint if no valid time slots
                if not valid_t2s:
                    return pyo.Constraint.Skip
                
                # Sum over valid (j, s, t) combinations
                return sum(sum(model.x[j, s, t]
                               for t in valid_t2s if (j, s, t) in model.x)
                           for (j, s) in model.PJ[p]) <= 1
            else:
                return pyo.Constraint.Skip
        
    # Create one constraint for each member of sets G and P
    model.Const4 = pyo.Constraint(model.P, model.T, rule=constraint_rule4)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Course j section s should be assigned to exactly one time slot
    def constraint_rule5(model, j, s, u):
        if j in model.J and s in model.JS[j] and u in model.JSU[j,s]:
            return sum(model.x[j,s,t] for t in model.UT[u]) == 1
        else:
            return pyo.Constraint.Skip
    
    # Create one constraint for each member of the set JS
    model.Const5 = pyo.Constraint(model.J, model.S, model.U, rule=constraint_rule5)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The schedule of math and chem courses is fixed

    def constraint_rule6(model, j, s, t):
        return model.x[j,s,t] == 1

    # Create one constraint for each member of set JST3
    model.Const6 = pyo.Constraint(model.JST3, rule=constraint_rule6)


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The schedule of courses preferred to be at certain times
    def constraint_rule7(model, j, s, t):
        if (j, s, t) in model.JST1:
            return model.x[j,s,t] == 1 - model.tau[j,s,t]
        else:
            return pyo.Constraint.Skip
    
    # Create one constraint for each member of set JST1
    model.Const7 = pyo.Constraint(model.JST1, rule=constraint_rule7)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The schedule of courses preferred not to be at certain times
    def constraint_rule8(model, j, s, t):
        if (j, s, t) in model.JST2:
            return model.x[j,s,t] == model.pi[j,s,t]
        else:
            return pyo.Constraint.Skip
    
    # Create one constraint for each member of set JST2
    model.Const8 = pyo.Constraint(model.JST2, rule=constraint_rule8)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The schedule of courses that cannot be at busy times of their professors
    def constraint_rule9(model, j, s, t):
        return model.x[j,s,t] == 0
    
    # Create one constraint for each member of set JST4
    model.Const9 = pyo.Constraint(model.JST4, rule=constraint_rule9)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # If a professor prefers not to teach on a specific day
    def constraint_rule10(model, j, s, t, p, d):
        if (p, d) in model.PD and (j, s) in model.PJ[p] and t in model.DT[d] and t in model.JST[j,s]:
            # Ensure x[j, s, t] is a valid index
            if (j, s, t) in model.x:
                # Check if (p, d) is in model.PD
                if (p, d) in model.PD:
                    return model.x[j, s, t] == model.gamma[j, s, t]
        else:
            return pyo.Constraint.Skip
    
    # Create one constraint for each member of sets JS, T, P, D
    model.Const10 = pyo.Constraint(model.J, model.S, model.T, model.P, model.D, rule=constraint_rule10)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # If professor p teaches on day d
    def constraint_rule11(model, p, d):
        return (sum(sum(model.x[j,s,t]
                        for t in model.DT[d] if t in model.JST[j,s])
                    for (j,s) in model.PJ[p])
               <= 10 * model.m[p, d])
        
    # Create one constraint for each member of sets D and P
    model.Const11 = pyo.Constraint(model.P, model.D, rule=constraint_rule11)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # the max number of days to teach in a week for professor p
    def constraint_rule12(model, p):
        return sum(model.m[p, d] for d in model.D) <= model.k[p] + model.delta[p]
    
    model.Const12 = pyo.Constraint(model.P, rule=constraint_rule12)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # A constraint to satisfy preferred teaching time in a day
    def constraint_rule13(model, p):
        if p in model.P:
            
            # skip the constraint if 'WholeDay' is the preference
            if 'WholeDay' in model.PR[p]:
                return pyo.Constraint.Skip
            
            return model.mu[p] == (
                sum(sum(model.x[j,s,t]
                        for t in model.JST[j,s])
                    for (j, s) in model.PJ[p]
                    if (j, s) not in model.JST5) -
                sum(sum(model.x[j,s,t]
                        for t in model.RT[r] if t in model.JST[j, s])
                    for (j, s) in model.PJ[p] if (j, s) not in model.JST5
                    for r in model.PR[p]) )
        
    # Create one constraint for each member of sets P and R
    model.Const13 = pyo.Constraint(model.P, rule=constraint_rule13)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule14(model, d, l, s, t):
        # Ensure z[l, t] is valid
        if (l, t) in model.z:
            # Get the time slots for the course l, section s, and day d
            day_slots = time_slots_dict[l][s].get(d, [])  # Safely retrieve slots for the day
            if t in day_slots:  # Ensure t exists in the time slots for the day
                t_index = day_slots.index(t)  # Find the index of t
                offset = model.q[l, s]  # Offset for course l, section s
                new_index = t_index + offset  # Calculate the new index after offset
    
                if new_index < len(day_slots):  # Ensure the new index is within bounds
                    t_new = day_slots[new_index]  # New time slot after offset
                    return model.z[l, t] <= model.x[l, s, t_new]
    
        return pyo.Constraint.Skip  # Skip invalid constraints
    
    model.Const14 = pyo.Constraint(model.D, model.L, model.S, model.T, rule=constraint_rule14)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule15(model, l, t):
        # Only include valid z indices
        if (l, t) in model.z:
            sections = model.LS[l]  # Sections for course l
            time_slots = time_slots_dict[l]  # Access the time slots dictionary for course l
            total_sections = len(sections)  # Number of sections for course l
    
            # Sum over sections
            expr = 0
            for s in sections:
                # Get the day from time slot 't' and the list of slots for that day
                day_slots = None
                for d in time_slots[s]:
                    if t in time_slots[s][d]:  # Check if 't' belongs to this day's slots
                        day_slots = time_slots[s][d]
                        break
    
                if day_slots is not None and t in day_slots:
                    t_index = day_slots.index(t)
                    offset = model.q[l, s]
                    new_index = t_index + offset
    
                    # Ensure new_index is within bounds
                    if new_index < len(day_slots):
                        t_new = day_slots[new_index]
                        expr += model.x[l, s, t_new]
    
            return model.z[l, t] >= expr - (total_sections - 1)
    
        return pyo.Constraint.Skip
    
    model.Const15 = pyo.Constraint(model.L, model.T, rule=constraint_rule15)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule16(model, l):
        # Use a set to collect unique time slots t for course l
        unique_times = set()
        for d in model.D:
            for t in model.DT[d]:
                if (l, t) in model.z:
                    unique_times.add((l, t))  # Add (l, t) to the set to ensure uniqueness
        
        # Sum z[l, t] for all unique time slots t
        return sum(model.z[l, t] for (l, t) in unique_times) == 1
    
    # Add the constraint to the model
    model.Const16 = pyo.Constraint(model.L, rule=constraint_rule16)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule17(model, p, d, t):
        if p in model.P and d in model.D and t in model.DT[d] and t in model.PT[p]:
            return (model.alpha[p,d,t] <= sum(model.x[j,s,t]
                                          for (j, s) in model.PJ[p]
                                          if (j, s) not in model.JST5
                                          and t in model.JST[j, s]))
        else:
            return pyo.Constraint.Skip
    
    model.Const17 = pyo.Constraint(model.P, model.D, model.T, rule=constraint_rule17)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule18(model, p, d, t):
        if p in model.P and d in model.D and t in model.DT[d] and t in model.PT[p]:
            
            # Get the numerical start time of t
            start_time_t = time_slot_info[t]['start']
            
            # Filter time slots t1 that occur earlier than or equal to t
            relevant_time_slots = [t1 for t1 in model.PT[p]
                                   if t1 in model.DT[d] and time_slot_info[t1]['start'] <= start_time_t]
            
            return (sum(model.alpha[p, d, t1]
                        for t1 in relevant_time_slots)
                    >=
                    sum(model.x[j, s, t]
                        for (j, s) in model.PJ[p]
                        if (j, s) not in model.JST5
                        and t in model.JST[j, s]))
        else:
            return pyo.Constraint.Skip
    
    model.Const18 = pyo.Constraint(model.P, model.D, model.T, rule=constraint_rule18)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule19(model, p, d):
    
        # Collect valid terms
        valid_terms = [model.alpha[p,d,t] for t in model.DT[d] if t in model.PT[p]]
    
        if valid_terms:  # Only create constraint if there are valid terms
            expr = sum(valid_terms)
            return expr <= 1
        else:  # If no valid terms, return a feasible constraint
            return pyo.Constraint.Feasible
    
    model.Const19 = pyo.Constraint(model.P, model.D, rule=constraint_rule19)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule20(model, p, d, t):
        if p in model.P and d in model.D and t in model.DT[d] and t in model.PT[p]:
            return (model.beta[p,d,t] <= sum(model.x[j,s,t]
                                          for (j, s) in model.PJ[p]
                                          if (j, s) not in model.JST5
                                          and t in model.JST[j, s]))
        else:
            return pyo.Constraint.Skip
    
    model.Const20 = pyo.Constraint(model.P, model.D, model.T, rule=constraint_rule20)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule21(model, p, d, t):
        if p in model.P and d in model.D and t in model.DT[d] and t in model.PT[p]:
            
            # Get the numerical start time of t
            start_time_t = time_slot_info[t]['start']
            
            # Filter time slots t1 that occur earlier than or equal to t
            relevant_time_slots = [t1 for t1 in model.PT[p]
                                   if t1 in model.DT[d] and time_slot_info[t1]['start'] >= start_time_t]
            
            return (sum(model.beta[p, d, t1]
                        for t1 in relevant_time_slots)
                    >=
                    sum(model.x[j, s, t]
                        for (j, s) in model.PJ[p]
                        if (j, s) not in model.JST5
                        and t in model.JST[j, s]))
        else:
            return pyo.Constraint.Skip
    
    model.Const21 = pyo.Constraint(model.P, model.D, model.T, rule=constraint_rule21)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule22(model, p, d):
        
        # Collect valid terms
        valid_terms = [model.beta[p,d,t] for t in model.DT[d] if t in model.PT[p]]
    
        if valid_terms:  # Only create constraint if there are valid terms
            expr = sum(valid_terms)
            return expr <= 1
        else:  # If no valid terms, return a feasible constraint
            return pyo.Constraint.Feasible
    
    model.Const22 = pyo.Constraint(model.P, model.D, rule=constraint_rule22)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule23(model, p, d):
        
        # Ensure only valid indices are considered
        valid_time_slots = [t for t in model.DT[d] if (p, d, t) in model.alpha]
        
        # Skip constraint creation if there are no valid time slots
        if not valid_time_slots:
            return pyo.Constraint.Skip
        
        return (model.w[p,d] == (sum(time_slot_info[t]['end'] * model.beta[p,d,t]
                                    for t in valid_time_slots) -
                                 sum(time_slot_info[t]['start'] * model.alpha[p,d,t]
                                     for t in valid_time_slots))  )
    
    model.Const23 = pyo.Constraint(model.P, model.D, rule=constraint_rule23)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule24(model, p, d):
        return model.w[p, d] <= model.h[p,d] + model.theta[p, d]
    
    model.Const24 = pyo.Constraint(model.P, model.D, rule=constraint_rule24)
    
    return model

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # print('The model has been built, and the solver is solving the model ...')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
# In[ ]:
    
def solve_model(model):
    
    # Create a solver instance using the GLPK solver
    # solver = SolverFactory('glpk', executable=r'C:\Users\rahda\anaconda3\Library\bin\glpsol.exe')
    solver = SolverFactory('glpk')
    
    # Set the termination limit in seconds
    solver.options['tmlim'] = 120        # Time limit
    # solver.options['mipgap'] = 0.05      # Allow a 5% optimality gap

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Call the solve method to solve the optimization problem defined by the Pyomo model
    solution = solver.solve(model, tee=False)

    # print('\nThe model has been solved -->', solution.Solver._list[0].Termination_condition)

    return model, solution


# In[ ]:


def write_result(x, y, output_file_path_excel):
    df1 = pd.DataFrame(x, columns=['Course', 'Section', 'Time'])
    df2 = pd.DataFrame(y, columns=['Group', 'Course', 'Time'])
    
    with pd.ExcelWriter(output_file_path_excel, engine='xlsxwriter') as writer:
        df1.to_excel(writer, sheet_name='x', index=False)
        df2.to_excel(writer, sheet_name='y', index=False)


# In[ ]:


def plot_schedule(x, output_file_path_schedule, lgnd=False, save=False):
    # Generate a color map for courses
#     colors = colormaps['jet'](np.linspace(0, 1, len(J_list)))
    jet = cm = plt.get_cmap('jet') 
    colors = jet(np.linspace(0, 1, len(J_list))) 
    course_colors = {c: colors[i] for i, c in enumerate(sorted(J_list))}

    # Define the preferred time ranges for visualization
    preferred_times = {
        'Morning': (8, 12),
        'ExtendedMorning': (8, 14),
        'WholeDay': (8, 16),
        'Midday': (10, 14),
        'ExtendedMidday': (10, 16),
        'Afternoon': (12, 16),
        'Evening': (16, 22)
    }
    
    # Create a figure with one subplot per professor
    num_groups = len(PJ_dict)
    fig, axes = plt.subplots((num_groups + 1) // 2, 2, figsize=(16, 6 * ((num_groups + 1) // 2)))
    axes = axes.flatten()  # Flatten to easily iterate over axes
    
    Days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
    days = ['M', 'T', 'W', 'R', 'F']
    day_index = {day: i for i, day in enumerate(days)}
    
    # Iterate over the professors and corresponding axes
    for idx, (p, ax) in enumerate(zip(P_list, axes)):
        schedule = []
        # Collect the schedule for the professor
        for [j, s, t] in x:
            if (j, s) in PJ_dict[p]:
                schedule.append((j, s, t))
                
        # Shade the preferred time range
        for preference in PR_dict[p]:
            start_pref, end_pref = preferred_times[preference]
            for day in days:
                ax.broken_barh([(day_index[day] - 0.5, 1)],
                               (start_pref, end_pref - start_pref),
                               facecolors='lightgray', alpha=0.7)
            
        if p in restrictions:
            for day in days:
                if day in restrictions[p]:
                    start_pref, end_pref = restrictions[p][day]
                    rect = Rectangle(
                        (day_index[day] - 0.4, start_pref), # Bottom-left corner
                        0.8,                                # Width
                        end_pref - start_pref,              # Height
                        facecolor='peachpuff',              # fill color
                        edgecolor='lightsalmon',            # Border color
                        hatch='xxxx',                       # Hatching pattern
                        linewidth=1                         # border line width
                    )
                    ax.add_patch(rect)
                    ax.text(day_index[day], start_pref + (end_pref - start_pref) / 2,
                            'Busy time', ha='center', va='center', color='k', fontsize=11)
    
        # Plot each course in the professor's schedule
        for course, section, time_slot in schedule:
            days_of_week = [day for day in days if time_slot_info[time_slot][day]]
            start_hour = time_slot_info[time_slot]['start']
            end_hour = time_slot_info[time_slot]['end']
            course_color = course_colors[course]
    
            for day in days_of_week:
                ax.broken_barh([(day_index[day] - 0.4, 0.8)],
                               (start_hour, end_hour - start_hour),
                               facecolors=(course_color,), ec='k', alpha=.4)
                ax.text(day_index[day], start_hour + (end_hour - start_hour) / 2,
                        f'{course}-{section}', ha='center', va='center', color='k')
    
        # Formatting the plot
        ax.set_xticks(range(len(Days)))
        ax.set_xticklabels(Days)
        ax.set_yticks(range(8, 22))
        ax.set_yticklabels([f'{hour}:00' for hour in range(8, 22)])
        ax.yaxis.grid(True)
        ax.set_ylim(8, 22)  # Show time from 8:00 AM to 10:00 PM
        ax.set_xlim(-0.5, 4.5)
        ax.invert_yaxis()
        ax.set_title(f'{p}\'s Teaching Schedule (' + PR_dict[p][0] + ')')

        # fig.suptitle('Course Schedules', fontsize=14, fontweight='bold')
    
    # Hide unused subplots
    for ax in axes[num_groups:]:
        ax.axis('off')
    
    if lgnd:
        # Adding a legend to the plot
        handles = [plt.Line2D([0], [0], color=color, lw=4) for color in course_colors.values()]
        labels = course_colors.keys()
        # Add the legend above the subplots
        fig.legend(handles, labels, loc='upper center', ncol=8)

    if save:
        # Save the plot to a file
        plt.savefig(output_file_path_schedule, dpi=300, bbox_inches='tight')
        plt.savefig(output_file_path_schedule+'.pdf', dpi=300, bbox_inches='tight')
    
    plt.show()


# In[ ]:


def plot_conflicts(x, y, output_file_path_conflicts, lgnd=False, save=False):

    # Find all conflicts for a single day
    def detect_conflicts(courses):
        num_courses = len(courses)
        is_conflict = [False] * num_courses  # Initialize conflict flags for all courses
    
        # Check each pair of courses for conflicts
        for i in range(num_courses):
            start1, end1, _, _ = courses[i]
            for j in range(num_courses):
                if i != j:  # Don't compare a course with itself
                    start2, end2, _, _ = courses[j]
                    if start1 < end2 and end1 > start2:  # Overlap condition
                        is_conflict[i] = True
                        is_conflict[j] = True  # Propagate the conflict flag
    
        return is_conflict
    
    # Generate a color map for courses
#     colors = colormaps['jet'](np.linspace(0, 1, len(J_list))) 
    jet = cm = plt.get_cmap('jet') 
    colors = jet(np.linspace(0, 1, len(J_list))) 
    course_colors = {c: colors[i] for i, c in enumerate(sorted(J_list))}
    
    # Create a figure with one subplot per group
    num_groups = len(FJ_dict)
    fig, axes = plt.subplots((num_groups + 1) // 2, 2, figsize=(16, 6 * ((num_groups + 1) // 2)))
    axes = axes.flatten()  # Flatten to easily iterate over axes
    
    Days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
    days = ['M', 'T', 'W', 'R', 'F']
    day_index = {day: i for i, day in enumerate(days)}
    
    full_width = 0.8  # Original width for non-conflicting boxes
    
    # Iterate over each group and corresponding subplot axis
    for idx, (group, ax) in enumerate(zip(FJ_dict, axes)):
        schedule = {}  # Store schedules as {day: [(start, end, course, section)]}
        
        # Collect schedules based on `y` and `x`
        for j in FJ_dict[group]:
            for t in JT_dict[j]:
                if [group, j, t] in y:  # If y[j,t] = 1, find corresponding section
                    selected_section = next((s for s in JS_dict[j] if [j, s, t] in x), None)
                    if selected_section:
                        days_of_week = [day for day in days if time_slot_info[t][day]]
                        start_hour = time_slot_info[t]['start']
                        end_hour = time_slot_info[t]['end']
                        for day in days_of_week:
                            if day not in schedule:
                                schedule[day] = []
                            schedule[day].append((start_hour, end_hour, j, selected_section))
    
        # Detect conflicts for each day
        for day, courses in schedule.items():
            # Determine conflict status for all courses on the same day
            is_conflict = detect_conflicts(courses)
            count = 0
            for i, (start, end, course, section) in enumerate(courses):
                course_color = course_colors[course]
                
                # Set box width and offset based on conflict status
                if is_conflict[i]:
                    width = full_width/sum(is_conflict)
                    offset = -0.4 + count * full_width/sum(is_conflict)
                    count += 1
                    face_color ='pink'
                    line_width = 3
                    edge_color = 'red'  # Red border for conflicts
                    font_color = 'black'
                    txt = course
                    font_size = 9
                else:
                    width = full_width
                    offset = -0.4
                    face_color = course_color
                    line_width = 1
                    edge_color = 'black'  # Default black border
                    font_color = 'black'
                    txt = f'{course}-{section}'
                    font_size = 10
    
                # Plot the course
                ax.broken_barh([(day_index[day] + offset, width)],
                               (start, end - start),
                               facecolors=face_color,
                               edgecolors=edge_color,
                               linewidth = line_width,
                               alpha=0.4)
                ax.text(day_index[day] + offset + width / 2, 
                        start + (end - start) / 2,
                        txt, ha='center', va='center',
                        color=font_color, fontsize = font_size)
        
        # Formatting the plot
        ax.set_xticks(range(len(Days)))
        ax.set_xticklabels(Days)
        ax.set_yticks(range(8, 22))
        ax.set_yticklabels([f'{hour}:00' for hour in range(8, 22)])
        ax.yaxis.grid(True)
        ax.set_ylim(8, 22)  # Show time from 8:00 AM to 10:00 PM
        ax.set_xlim(-0.5, 4.5)
        ax.invert_yaxis()
        ax.set_title(f'{FN_dict[group]} courses in {group} that cannot be conflicted')
    
    # fig.suptitle('Conflicting Courses', fontsize=14, fontweight='bold')
        
    # Hide unused subplots
    for ax in axes[num_groups:]:
        ax.axis('off')
        
    
    if lgnd:
        # Adding a legend to the plot
        handles = [plt.Line2D([0], [0], color=color, lw=4) for color in course_colors.values()]
        labels = course_colors.keys()
        # Add the legend above the subplots
        fig.legend(handles, labels, loc='upper center', ncol=8)
    
    if save:
        # Save the plot to a file
        plt.savefig(output_file_path_conflicts+'.png', dpi=300, bbox_inches='tight')
        plt.savefig(output_file_path_conflicts+'.pdf', dpi=300, bbox_inches='tight')
        
    plt.show()


# In[ ]:


def plot_meetings(x, output_file_path_meetings, save=False):

    possible_meeting_times = []
    aa = 0
    for m in M_list:
        for p in P_list:
            a = sum(1 for j, s, t in x if t in MT_dict[m] and t in JST_dict[j, s] and (j, s) in PJ_dict[p])
            aa += a
        if aa == 0:
            possible_meeting_times.append(m)
        aa = 0

    # Combine all professor schedules into a single plot
    Days = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday']
    days = ['M', 'T', 'W', 'R', 'F']
    day_index = {day: i for i, day in enumerate(days)}
    day_map = {'M': 'Mon', 'T': 'Tue', 'W': 'Wed',  'R': 'Thu', 'F': 'Fri'}
    
    # Initialize the figure and axis
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title("Course Schedule and Possible Meeting Times")
    
    # Collect all schedules across professors
    all_schedules = []  # (course, section, time_slot, day, start, end)
    for p in P_list:
        for [j, s, t] in x:
            if (j, s) in PJ_dict[p]:
                days_of_week = [day for day in days if time_slot_info[t][day]]
                start_hour = time_slot_info[t]['start']
                end_hour = time_slot_info[t]['end']
                for day in days_of_week:
                    all_schedules.append((j, s, t, day, start_hour, end_hour))
    
    # Plot all courses (no color maps, uniform color)
    for course, section, time_slot, day, start, end in all_schedules:
        ax.broken_barh([(day_index[day] - 0.4, 0.8)],
                       (start, end - start),
                       facecolors="lightblue", edgecolor='gray', alpha=0.5)
    
    # Highlight possible meeting times
    for time_slot in possible_meeting_times:
        days_of_week = [day for day in days if meeting_info[time_slot][day]]
        start_hour = meeting_info[time_slot]['start']
        end_hour = meeting_info[time_slot]['end']

        # Extract the day and time (in decimal format)
        d = time_slot[0]
        f = float(time_slot[1:])
    
        # Convert decimal time to hours and minutes
        h = int(f)
        m = int((f - h) * 60)  # Get the minutes from the decimal part
    
        # Format the output string
        txt = f"{day_map[d]} {h:02}:{m:02}"

        for day in days_of_week:
            ax.broken_barh([(day_index[day] - 0.4, 0.8)],
                           (start_hour, end_hour - start_hour),
                           facecolors='lightgreen', edgecolor='k', alpha=0.5)
            
            ax.text(day_index[day], start_hour + (end_hour - start_hour) / 2,
                    txt, ha='center', va='center', fontsize=10, color='k')
    
    # Formatting the plot
    ax.set_xticks(range(len(Days)))
    ax.set_xticklabels(Days)
    ax.set_yticks(range(8, 22))
    ax.set_yticklabels([f'{hour}:00' for hour in range(8, 22)])
    ax.yaxis.grid(True)
    ax.set_ylim(8, 22)  # Show time from 8:00 AM to 10:00 PM
    ax.set_xlim(-0.5, 4.5)
    ax.invert_yaxis()
    ax.set_ylabel('Time of Day')
    ax.set_xlabel('Day of the Week')
    
    # Add a legend
    handles = [
        plt.Line2D([0], [0], color="lightblue", lw=6, alpha=0.7, label="Courses"),
        plt.Line2D([0], [0], color="lightgreen", lw=6, alpha=0.5, label="Possible Meeting Times")
    ]
    ax.legend(handles=handles, loc='lower center', ncol=2)

    if save:
        plt.savefig(output_file_path_meetings+'.png', dpi=300)
        plt.savefig(output_file_path_meetings+'.pdf', dpi=300)

    plt.show()


# In[ ]:


def read_schedule(output_file_path_excel):
    df1 = pd.read_excel(output_file_path_excel, sheet_name='x')
    x = df1.values.tolist()
    return x


# In[ ]:


def build_check_model(data, x):
    
    for key, value in data.items():
        globals()[key] = value
    
    model = pyo.ConcreteModel() # create a concrete model using Pyomo
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The indices for constraints and variables (Sets)
    model.T = pyo.Set(initialize=T_list)  # set of time slots (T)
    model.B = pyo.Set(initialize=B_list)  # set of block times (B)
    model.J = pyo.Set(initialize=J_list)  # set of courses (J)
    model.S = pyo.Set(initialize=S_list)  # set of sections (S)
    model.D = pyo.Set(initialize=D_list)  # set of weekdays (D)
    model.U = pyo.Set(initialize=U_list)  # set of class durations (U)
    model.P = pyo.Set(initialize=P_list)  # set of professors (P)
    model.F = pyo.Set(initialize=F_list)  # set of incompatible course groups (F)
    model.PD   = pyo.Set(initialize=PD_list)   # set of professors and days they don't want to teach
    model.JS1  = pyo.Set(initialize=JS_list)   # set of courses and sections (tuple)

    model.JS = pyo.Set(model.J, within=model.S, initialize=JS_dict)   # set of courses and sections (dic)
    model.JSU= pyo.Set(model.JS1, within=model.U, initialize=JSU_dict) # set of our courses, sections, and durations
    model.JU = pyo.Set(model.J, within=model.U, initialize=JU_dict)    # set of all courses and possible durations
    model.JST= pyo.Set(model.JS1, within=model.T, initialize=JST_dict) # set of all courses, sections, and possible times
    model.JT = pyo.Set(model.J, within=model.T, initialize=JT_dict)    # set of all courses and possible times
    model.PT = pyo.Set(model.P, within=model.T, initialize=PT_dict)    # set of time slots for professor p
    model.UT = pyo.Set(model.U, within=model.T, initialize=UT_dict)    # set of time slots for duration u
    model.DT = pyo.Set(model.D, within=model.T, initialize=DT_dict)    # set of time slots in day d
    model.FJ = pyo.Set(model.F, within=model.J, initialize=FJ_dict)    # set of incompatible courses
    model.BT = pyo.Set(model.B, within=model.T, initialize=BT_dict)    # set of blocks and time slots conflicts

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create an empty list to hold the variable x indices
    x_ind_list = []
    
    # Loop through each (j, s) pair in JS and corresponding timeslots in JST
    for (j, s) in JS_list:
        for t in model.JST[(j, s)]:
            # Initialize the variable for each valid (j, s, t) combination
            x_ind_list.append((j, s, t))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # A function to create variable y indices
    def valid_index_y(model):
#         return ((j, t) for j in model.J for t in model.JT[j])
        return ((f, j, t) for f in model.F for j in model.FJ[f] for t in model.JT[j])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # A function to create variable u indices
    valid_u_indices = {(p, d, t)
                       for p in model.P
                       for d in model.D
                       for t in model.DT[d] if t in model.PT[p]}

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Variables of the model
    
    # Binary variable x(j,s,t): it is 1 if course j, section s, scheduled in time t
    model.x = pyo.Var(x_ind_list, domain=pyo.Binary)
    
    # Binary variable y(j,t) is 1 if course j has an available section at time t
    model.y = pyo.Var(valid_index_y, domain=pyo.Binary)
    
    # Integer variable eta(f,b) indicates the number of conflicted courses in group f block b
    model.eta = pyo.Var(model.F*model.B, domain=pyo.NonNegativeReals)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Define the objective function expression
    def obj_expression(model):
        return  pyo.summation(model.eta)
    
    model.Obj = pyo.Objective(rule=obj_expression, sense=pyo.minimize)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule1(model, f, j, t):
#         if j in model.J and t in model.JT[j]:
        if f in model.F and j in model.FJ[f] and t in model.JT[j]:
            return model.y[f,j,t] <= sum(model.x[j,s,t]
                                         for s in model.JS[j]
                                         if (j, s) in model.JST and t in model.JST[j, s])
        else:
            return pyo.Constraint.Skip
    
    model.Const1 = pyo.Constraint(model.F, model.J, model.T, rule=constraint_rule1)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule2(model, f, j):
        if j in model.FJ[f]:
            return (sum(sum(model.y[f,j,t]
                            for t in model.UT[u])
                        for u in model.JU[j])
                    == 1)
        else:
            return pyo.Constraint.Skip
    
    model.Const2 = pyo.Constraint(model.F, model.J, rule=constraint_rule2)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule3(model, f, b):
        if f in model.F and b in model.B:
            # Collect valid terms
            valid_terms = [model.y[f, j, t] for j in model.FJ[f] for t in model.BT[b] if t in model.JT[j]]
            
            if valid_terms:  # Only create constraint if there are valid terms
                expr = sum(valid_terms)
                return expr <= 1 + model.eta[f,b]
            else:  # If no valid terms, return a feasible constraint
                return pyo.Constraint.Feasible
        else:
            return pyo.Constraint.Skip
    
    model.Const3 = pyo.Constraint(model.F, model.B, rule=constraint_rule3)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def constraint_rule4(model, j, s, t):
        if [j, s, t] in x:
            return model.x[j,s,t] == 1
        else:
            return model.x[j,s,t] == 0
    
    model.Const4 = pyo.Constraint(x_ind_list, rule=constraint_rule4)
 

    
    # Create a solver instance using the GLPK solver
    solver = SolverFactory('glpk')
    
    # Set the termination limit in seconds
    solver.options['tmlim'] = 120

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Call the solve method to solve the optimization problem defined by the Pyomo model
    solution = solver.solve(model, tee=False)

    return model, solution


# In[ ]:

def save_uploaded_file(upload_widget, saved_name=None, mode='scheduler'):
    """
    Save uploaded FileUpload content to disk and return the filename.
    
    Parameters
    ----------
    upload_widget : ipywidgets.FileUpload
        The widget holding the uploaded file.
    saved_name : str, optional
        Name to save file as. If None, a default will be used based on mode.
    mode : str, {'scheduler', 'check'}
        Determines default saved filename if saved_name not provided.
    """
    if not upload_widget.value:
        return None
    
    # Determine filename
    if saved_name is None:
        if mode == 'scheduler':
            saved_name = 'Courses_info_run.xlsx'
        elif mode == 'check':
            saved_name = 'Schedule_results_run.xlsx'
        else:
            raise ValueError("mode must be either 'scheduler' or 'check'")

    # Handle both new (tuple) and old (dict) behaviors
    if isinstance(upload_widget.value, dict):
        file_info = list(upload_widget.value.values())[0]
    else:
        file_info = upload_widget.value[0]

    content = file_info['content']
    
    with open(saved_name, 'wb') as f:
        f.write(content)
    return saved_name



# In[ ]:


def run_pipeline(input_filename):
    """Run your scheduler pipeline and save Excel + images."""
    
    # 1) Read input
    print("\n📥 Reading data...")
    sheetnames = what_semester(input_filename)
    data = read_data(input_filename, sheetnames)

    # 2) Build & solve
    print("🛠️ Building the model... (this may take a moment)")
    model = build_model(data)
    print("🧮 Solving the model...\n")
    model, solution = solve_model(model)
    
    # Show the status of the solution and if it is optimal or not
    result = solution.Solver._list
    print('     Status:   ', result[0].Status)
    print('     Condition:', result[0].Termination_condition)
    print('     Time:      %.2f seconds\n' %result[0].Time)

    # 3) Extract results
    print("📊 Extracting results...")
    x, y = [], []
    for j, s, t in model.x:
        if model.x[j,s,t].value == 1:
            x.append([j,s,t])
    for f, j, t in model.y:
        if model.y[f,j,t].value == 1:
            y.append([f,j,t])

    # 4) Write Excel + plots
    print("💾 Writing results to Excel and generating plots...\n")
    excel_out = 'Schedule_results.xlsx'
    write_result(x, y, excel_out)

    # Show Excel download link
    print("🎉 Scheduling complete.")
    print("⬇️ Download the Excel results here:\n")
    display(FileLink(excel_out))

    print('\n🖼️ Find the scheduling plots here:\n')    
    plot_schedule(x, 'schedule.png')
    print('\n\n\n')
    plot_conflicts(x, y, 'conflicts.png')
    print('\n\n\n')
    plot_meetings(x, 'meetings.png')
    
    return excel_out


# In[ ]:

def run_check_pipeline(schedule_filename, input_filename="Courses_info_run.xlsx"):
    """Check a user-modified schedule against constraints."""
    # Read data
    sheetnames = what_semester(input_filename)
    data = read_data(input_filename, sheetnames)

    # Read modified schedule
    x = read_schedule(schedule_filename)

    # Build & solve check model
    print("🔍 Checking modified schedule...")
    check_model, check_solution = build_check_model(data, x)

    y = []
    for f, j, t in check_model.y:
        if check_model.y[f,j,t].value == 1:
            y.append([f, j, t])

    # Plots
    plot_schedule(x, "schedule_checked")
    plot_conflicts(x, y, "conflicts_checked")
    plot_meetings(x, "meetings_checked")

    print("✅ Check complete. See plots above.")

# In[ ]:

def create_download_link(filepath, label=None, mime=None):
    """
    Create and display a download link for `filepath` by embedding it as a
    base64 data URI. Works reliably in Binder/Voila.
    """
    if not os.path.exists(filepath):
        print(f"❌ File not found: {filepath}")
        return

    if label is None:
        label = f"⬇️ Download {os.path.basename(filepath)}"

    # Guess mime type if not provided
    if mime is None:
        ext = os.path.splitext(filepath)[1].lower()
        if ext in ('.xlsx', '.xlsm', '.xls'):
            mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        elif ext in ('.csv',):
            mime = "text/csv"
        elif ext in ('.png', '.jpg', '.jpeg'):
            mime = f"image/{ext[1:]}"
        else:
            mime = "application/octet-stream"

    with open(filepath, "rb") as f:
        b = f.read()
    b64 = base64.b64encode(b).decode()
    href = f"data:{mime};base64,{b64}"
    html = f'<a download="{os.path.basename(filepath)}" href="{href}">{label}</a>'
    display(HTML(html))


# In[ ]:

def make_download_html(path, label=None, mime=None):
    """Return an IPython HTML object with an <a download> link that embeds the file in base64.
    Works in Voila/Binder because it does not rely on server file serving."""
    if not os.path.exists(path):
        return HTML(f"<b>File not found: {path}</b>")
    label = label or f"⬇️ Download {os.path.basename(path)}"
    with open(path, "rb") as f:
        data = f.read()
    if mime is None:
        # common Excel mime
        if path.endswith(".xlsx"):
            mime = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        elif path.endswith(".xls"):
            mime = "application/vnd.ms-excel"
        elif path.endswith(".png"):
            mime = "image/png"
        else:
            mime = "application/octet-stream"
    b64 = base64.b64encode(data).decode()
    href = f"data:{mime};base64,{b64}"
    # target=_blank recommended for some browsers
    html = f'<a download="{os.path.basename(path)}" href="{href}" target="_blank">{label}</a>'
    return HTML(html)


# In[ ]:

