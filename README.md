# MN_benchmakring
The GitHub repository for the paper: Network Topology Evaluation and Construction for Molecular Networking

## Project Structure 📂

Below is the hierarchical outline of the key directories and files within this project:

```plaintext
📦 project_directory
 ┣ 📂 data
 ┃ ┣ 📂 converted
 ┃ ┣ 📂 merged_pairs
 ┃ ┣ 📂 MS2DeepScoreModel
 ┃ ┣ 📂 raw_ms2deepscore
 ┃ ┣ 📂 Network_barebone
 ┃ ┗ 📂 summary
 ┣ 📂 src
 ┃ ┣ 📂 Benchmarking
 ┃ ┣ 📂 Transirive_Alignment
 ┃ ┣ 📂 MS2DeepScore
 ┃ ┗ 📂 Plot
 ┗ 📂 results
   ┣ 📂 alignment_results
   ┣ 📂 MS2DeepScore_Network
   ┣ 📂 results-baseline
   ┣ 📂 results-cast
   ┣ 📂 results-MS2DeepScore
   ┣ 📂 results-re
   ┣ 📂 results-re-cast
   ┗ 📂 results-ClassyFire
📜 requirement.txt
📜 README.md
```

## Installation Guide 🛠️

Before running the project, ensure you have Python installed on your system. This project requires Python 3.8 or newer. You can download Python from [python.org](https://www.python.org/downloads/).

### Setting Up a Virtual Environment

It's recommended to use a virtual environment for Python projects. This keeps dependencies required by different projects separate by creating isolated environments for them. Here's how you can set it up:

1. Open a terminal or command prompt.
2. Navigate to the project's root directory.
3. Run the following command to create a virtual environment named `MN_benchmarking`:

```bash
python3 -m venv MN_benchmarking
source MN_benchmarking/bin/activate
pip install -r requirements.txt
```

## Prepare the Testing Dataset 📚

Before you can start the benchmarking process, you need to prepare the dataset that will be used for testing. This involves either downloading a pre-existing dataset or generating your own raw data in the required format. Follow the steps below to prepare your testing dataset:

### Downloading the Pre-existing Dataset

1. Visit [Zenodo](https://zenodo.org/records/10724765) to download the pre-prepared benchmarking dataset.
2. After downloading, extract the dataset to your local machine.

### Generating Your Own Dataset

If you prefer to use your own raw data for benchmarking, follow the GNPS2 networking_barebone_workflow to generate data in the correct format. Here's a brief overview of the steps involved:

1. Collect your raw mass spectrometry data.
2. Process the data using the GNPS2 networking_barebone_workflow. Detailed instructions for this workflow can be found on the GNPS documentation website.
3. Ensure the data is in the appropriate format for the benchmarking process.

### Organizing the Dataset

Once you have your dataset ready, organize it into the correct directories within the project structure. Here's how to organize the files:

- **Merged Pairs**: Place all merged pairs files in the `data/merged_pairs` directory.
- **Converted Data**: Put all converted data files in the `data/converted` directory.
- **Summary Data**: Summary data files should go into the `data/summary` directory.

This organization is crucial for the benchmarking scripts to correctly locate and process the data. Ensure that each file type is placed in its respective directory as outlined above.

### Finalizing the Dataset Preparation

After organizing the dataset into the appropriate directories, you're ready to proceed with the benchmarking guide outlined in the previous sections. Ensure that your `input_library.txt` file reflects the datasets you've prepared for testing.

By following these steps, you'll have a well-organized and ready-to-use testing dataset for benchmarking the performance of different methods on your specific datasets.

## Transitive Alignment 🔄

Transitive Alignment is a novel approach proposed in our paper, "Network Topology Evaluation and Construction for Molecular Networking". This method focuses on re-aligning two spectra of molecules that may have undergone multiple modifications, leveraging the molecular network topology to enhance the accuracy of these alignments.

The concept behind Transitive Alignment is to utilize the inherent relationships within a molecular network to facilitate the alignment of spectra. This approach helps in accounting for the variations and modifications that molecules can undergo, providing a more robust framework for molecular networking analysis.

### Running Transitive Alignment

To perform Transitive Alignment on your dataset, follow these steps:

1. Prepare your dataset and list all the items you wish to align in a text file, such as `input_library.txt`. This file should be placed in the designated input library folder.

2. Execute the following command from the terminal, making sure to replace `PATH_TO_PROJECT_FOLDER` and `PATH_TO_INPUT_LIBRARY_FOLDER` with the actual paths to your project folder and the folder containing your `input_library.txt` file, respectively:

```bash
python3 PATH_TO_PROJECT_FOLDER/src/Transitive_Alignment/alignment.py --input PATH_TO_INPUT_LIBRARY_FOLDER/input_library.txt
```
After running the command, the alignment results will be stored in `results/alignment_results`. 

## Generate the MS2DeepScore Network 🌐

 For this project, we utilize a retrained MS2DeepScore model that is specifically tailored for the dataset mentioned in our paper, excluding the benchmarking dataset to ensure accurate evaluation.

### Downloading the Retrained MS2DeepScore Model

First, download the retrained MS2DeepScore model file from [Zenodo](https://zenodo.org/records/10724765). It's important to note that this retrained model is designed for use with the dataset mentioned in the paper and may not yield accurate results for other datasets without retraining.

### Using Your Own Data

If you wish to analyze your own data, you have two options:

1. Retrain the MS2DeepScore model with your dataset.
2. Use the original MS2DeepScore model published by the developers.

Retraining the model would be necessary to adapt to the specifics of your data to avoid overfitting.

### Running the MS2DeepScore Network Generation

To generate the MS2DeepScore Network, execute the following command in your terminal. Make sure to replace `PATH_TO_PROJECT_FOLDER`, `PATH_TO_INPUT_MGF_FILE_FOLDER/NAME_FOR_MGF_FILE.mgf`, and `PATH_TO_INPUT_MS2DEEPSCORE_MODEL_FOLDER/NAME_FOR_MODEL.hdf5` with the actual paths and filenames relevant to your project setup:

```bash
python3 PATH_TO_PROJECT_FOLDER/src/MS2DeepScore/Network_ms2deepscore.py -m PATH_TO_INPUT_MGF_FILE_FOLDER/NAME_FOR_MGF_FILE.mgf -w PATH_TO_INPUT_MS2DEEPSCORE_MODEL_FOLDER/NAME_FOR_MODEL.hdf5
```

The results from this operation will be stored in the `results/MS2DeepScore_Network` folder. 

## Benchmarking Guide 🚀

This section outlines the steps to benchmark the performance of different methods on the datasets specified in `input_library.txt`.

### Preparing the Datasets

First, list all the datasets you want to benchmark in a text file named `input_library.txt`. Place this file in the `PATH_TO_INPUT_LIBRARY_FOLDER` directory.

### Running Benchmarks

To benchmark the various methods, navigate to your project directory and execute the corresponding commands in the terminal.

#### Classic Method

```bash
python3 PATH_TO_PROJECT_FOLDER/scr/benchmarking/benchmark_classic.py --input PATH_TO_INPUT_LIBRARY_FOLDER/input_library.txt
```

#### CAST Method

```bash
python3 PATH_TO_PROJECT_FOLDER/scr/benchmarking/benchmark_cast.py --input PATH_TO_INPUT_LIBRARY_FOLDER/input_library.txt
```

#### CAST + Alignment Method

```bash
python3 PATH_TO_PROJECT_FOLDER/scr/benchmarking/benchmark_re_cast.py --input PATH_TO_INPUT_LIBRARY_FOLDER/input_library.txt
```

#### MS2DeepScore or Baseline Method

```bash
python3 PATH_TO_PROJECT_FOLDER/scr/benchmarking/benchmark_baseline.py --input NAME_FOR_THE_LIBRARY --mehtod "baseline" or "MS2DeepScore"
```

Make sure to replace `PATH_TO_PROJECT_FOLDER` and `PATH_TO_INPUT_LIBRARY_FOLDER` with the actual paths to your project and input library folder, respectively.

### Analyzing the Results

After running the benchmarks, analyze the output to compare the performance of the Classic, CAST, and CAST + Alignment methods. This comparison will help identify which method performs best for your specific datasets.

## Visualize the Benchmarking Results 📈

After successfully running the benchmarks, you'll find the results for each method stored in dedicated folders named `results-METHOD_NAME`. To proceed with visualizing these results and comparing the performance across different methods, follow these steps:

### Step 1: Prepare Your Results

Ensure that the results are correctly placed in their respective folders, such as `results-classic`, `results-cast`, and `results-re-cast`.

### Step 2: Adjust the Visualization Script

Navigate to the `benchmarking_plot.py` script located in `/src/plot/`. Adjust the dataset names within this script to align with those present in your benchmark results. This ensures that the script accurately retrieves and processes the data for visualization.

### Step 3: Generate the Plots

Run the `benchmarking_plot.py` script to visualize the results. Use the following command, adjusting the path as necessary:

```bash
python3 /path/to/MN_benchmarking/src/plot/benchmarking_plot.py
```

### Example Output

An example output figure will be saved to `results/plot/demo.png`, showcasing the comparative performance of the benchmarked methods. 

![Benchmarking Results Example](results/plot/demo.png)