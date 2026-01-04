# SeqNado Web Interface

SeqNado includes a web-based graphical interface to simplify managing genomes, generating experiment designs, and running workflows.

## Starting the Web App

To start the web server, run the following command in your terminal:

```bash
seqnado-web
```

This will launch the application on port `5001` by default.

### Connecting via SSH Tunnel

If you are running SeqNado on a remote server (e.g., an HPC cluster), you will need to set up an SSH tunnel to access the interface from your local browser.

Run this command on your **local machine**:

```bash
ssh -L 5001:localhost:5001 <your-username>@<server-address>
```

Replace `<your-username>` and `<server-address>` with your actual login details (e.g., `user@hpc.edu`).

Once the tunnel is established, open your web browser and navigate to:

[http://localhost:5001](http://localhost:5001)

## Features

### Dashboard
The dashboard provides a central hub for navigating the application. You can quickly access the Experiment Runner, Genome Manager, Configs, and Designs.

### Genome Management
Manage your reference genomes easily.
- **Add Genomes**: Configure paths to FASTA, GTF, and index files (STAR, Bowtie2).
- **Auto-Prediction**: If you leave the "Organism" field blank, SeqNado will attempt to guess it from the genome name (e.g., "hg38" -> "Homo sapiens").
- **Browse Files**: Use the built-in file browser to select files directly from the server.

### Design Editor
Create and modify experiment design CSV files without leaving the browser.
- **Generate from FASTQs**: Point the app to a folder of FASTQ files. It will automatically discover samples, pair reads (R1/R2), and populate metadata columns based on SeqNado's heuristics.
- **Interactive Editing**: Edit the design matrix in a spreadsheet-like view. Add columns (Condition, Replicate, etc.) and save changes directly to the server.
- **Scan for Designs**: Find existing design files in your project directories.

### Experiment Runner
(Work in Progress)
Configure and submit SeqNado workflows directly from the UI. Select your assay, project name, genome, and design file to generate the necessary `config.yaml` and execution commands.
