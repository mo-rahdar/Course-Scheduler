# 🎓 Course Scheduler Web App

This project provides an interactive **Course Scheduling** and **Schedule Checker** tool built with Pyomo and Voila.  
It allows you to generate feasible schedules based on course, professor, and their preferences — and later validate or adjust them.

---

## 🚀 Launch the Apps on Binder

Click the badges below to launch the apps directly in your browser. No installation required (powered by [Binder](https://mybinder.org)):

- **Course Scheduler** → Upload course info and generate a schedule.  
- **Schedule Checker** → Upload a modified schedule to check for conflicts.  
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mo-rahdar/Course-scheduler/HEAD?urlpath=voila/render/Course_Scheduling7.ipynb)

---

## 📂 Files

- `Course_Scheduling7.ipynb` → Main scheduling app  
- `Utils.py` → Helper functions (data reading, solving, plotting, etc.)  
- `requirements.txt` → Dependencies for Binder
- `apt.txt` → to install solver

---

## 🛠 How to Use

### 1. Run the Scheduler
1. Click the **Course Scheduler** badge above.  
2. Upload your `Courses_info.xlsx` file (with course, professor, and other data).  
3. Click **Run Scheduler**.  
4. The system will:  
   - Generate a feasible schedule.  
   - Save results to `Schedule_results.xlsx`.  
   - Show plots for course allocations, conflicts, and meeting times.  

---

### 2. Check a Modified Schedule
1. If you want to make manual changes:  
   - Download `Schedule_results.xlsx`.  
   - Edit it locally (e.g., swap times).  
2. Click the **Schedule Checker** badge above.  
3. Upload your modified `Schedule_results.xlsx`.  
4. The system will:  
   - Validate the new schedule.  
   - Highlight any conflicts (e.g., overlapping courses).  
   - Re-plot the updated schedule and conflicts.

