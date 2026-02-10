# micb475_project2
## 2026WT: group 9 MICB475 project 2 repository

## Feb 10, 2026
## Team Meeting 3 Agenda
- Review the demux.qzv generated for the Parkinson's dataset
- Review the table.qzv generated for the Parkinson's dataset
- Assignment 6 discussion
  - Option I: Matthew, Mikaela
  - Option II: Jaclyn, Ella
  - Option III: Narjis
- Discuss proposal plan
  - We had a group call and set up a document for our proposal. We started to discuss roles and timeline and hope to discuss this further in the meeting
  - We aim to do an overview of the literature to set up our research question and aims based on the pre-existing data
  - We're hoping to lay out our proposal aims and discuss the details of each analysis discussed in the last team meeting


### Attendance: Narjis, Matthew, Jaclyn, Mikaela, Ella, Evelyn, Ritu

## Team Meeting 3 Notes
## Proposal notes
- For the proposal, people writing different parts, set a deadline for group and go through each other's parts to make sure the tone/flow of language is coherent.
- RQ section:
  - Mention in bold, then write paragraph about hypothesis and background (could be different diseases; can refer to intro again).
- Title of proposal should be 'investigating' or 'exploring' rather than a conclusion.
- Background:
  - Mention Parkinson's, what this looks like pathologically and what the severity scores indicate. 
- Dataset discussion:
  - Overview of individuals in the data (age, etc.)
- Aims:
  - Mention that we will be looking at x metric, x statistical test, ensure to write everything explicitly.
  - Initial high-level analysis, aim is to look at different levels of severity and gut microbiome diversity.
  - When we look into the data further (Core Microbiome, Species Indicator).
    - Core Microbiome = differences between groups.
    - Species Indicator = specific species associated/rare species of bacteria present in conditions.
    - DESeq = microbiome influences species present (upregulation/downregulation).
    - Functional analysis? Evelyn mentioned we will determine this later if we have certain findings. 
  - For example, looking at Core Microbiome across different severity levels (comparative).
- Analysis: 
    - Run both alpha and beta diversity on all 3 variables selected.
    - May all be significant, can decide what to move forward with based on this.
    - Can find studies that look at the variables that bin the scores.
    - Coding should be done up to rarefaction.
    - Comparing between groups with continuity - may need to categorize afterwards, but first one is correlative analysis with alpha diversity.
    - Be explicit when we write this - the reader needs to know we are categorizing afterwards. 

## Feb 3, 2026
## Team Meeting 2 Agenda
- Reviewing Datasets & Research Questions
- Planning out project 2

### Attendance: Narjis, Matthew, Jaclyn, Mikaela, Ella, Evelyn, Ritu

## Team Meeting 2 Notes
- Reading over the options for research questions
- Using Parkinson's dataset, all continuous variables:
  - MOCA: cognitive
  - Duration of disease
  - UPDRS: motor function

**Research Question: Is there a correlation between the severity of Parkinson's disease and the gut microbiome?
**
- Can also make the question more specific based on our variables
  
### Project Experimental Plan
- Basic alpha diversity continuous plot
- Not time-based so can't really do a longitudinal
- Take each of the variables and do the alpha diversity correlation (3 separate plots for each metric)
1. Run an alpha diversity first as a correlation
- this step will define the rest + the binning depending on what we see in the preliminary analysis
2. Bin your categories (based on what we find in 1, will be the most challenging part since the 3 categories are confouding with eachother, can compare to the healthy patient profile here as a baseline but would not include healthy in the linear correlation)
3. Beta diversity
4. Core microbiome (Venn diagram of how many are shared between the most and least severe)
5. Indicator taxa
6. Deseq (look at how the micobes fluctuate in terms of abundance)

### Post-Meeting Action Items
1. Team Project Proposal, flesh it out by next meeting: Feb 10th
2. Do the QIIME2 processing before (how many samples we have etc.): Get it to the table.qzv point
   - Be sure to email Ritu with how we're going to trim
   - Tips: use a detached screen, let each other know who's working on each file, keep a copy of the original metadata file and make changes to a different copy,     upload everything to GitHub repo (any, .qza, .qzv)
3. Server credentials for group project from Evelyn


## Jan 27, 2026
## Team Meeting 1 Agenda
- Introductions
- Course logistics
- Datasets

### Attendance: Narjis, Matthew, Jaclyn, Mikaela, Ella, Evelyn, Ritu

## Team Meeting 1 Notes
- Add Evelyn and Ritu to the repo
- Broadly interested in health topics, mental health, virology, immunology etc over environmental microbiomes
- Depression dataset: we are the first ones analyzing it, could be tricky if there are any hiccups
  - Will probably need to do some extra psych research
  - Current depression dataset has HIV/HCV coinfection but patient group not specified in the metadata so we shouldn't use it 
- Have a research question decided by the end of next meeting (**Tues Feb 3**)
- Can find the paper associated with the dataset and send it to Evelyn by the end of **Friday Jan 30**, she'll look at the data and metadata
- Have a backup of another dataset to use in case we can't find one

### Post-Meeting Action Items
1. Send papers with datasets to Evelyn: by EOD **Friday Jan 30**
2. Set research question: by end of team meeting **Tuesday Feb 3**
