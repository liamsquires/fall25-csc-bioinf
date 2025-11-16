Assignment Description:

Your job is to implement the following steps:

    1.fetch some single-cell data and create a cell-gene expression matrix and the associated AnnData object;
    2.cluster the data to understand the cell structure; and
    3.annotate the data with biologically relevant information.

To do this, you should:

    1.Get the sample data here: https://app.box.com/s/lx2xownlrhz3us8496tyu9c4dgade814. This data contains the single-cell FASTQ files, as well as the reference genome (chr5) and the GTF file with transcript information.

    The list of whitelist barcodes is available here (https://github.com/f0t1h/3M-february-2018/raw/refs/heads/master/3M-february-2018.txt.gz).

    2. Use Alevin-fry to align and quantify this data to the reference genome. See here (https://www.sc-best-practices.org/introduction/raw_data_processing.html#a-real-world-example) for the installation steps.

    3. Perform cell clustering (via Leiden modularity algorithm). Output the clustering plot.

    4. Perform automatic cell annotation via CellTypist. Annotate the plot with the cell types.

    Yes, all steps are based on Single-cell Best Practices Book (https://www.sc-best-practices.org/). Feel free to use it! You will still need to automate the process and ensure that it can be run.

Included files:
- tempfiles/toy_read_ref_set.tar.gz (from https://app.box.com/s/lx2xownlrhz3us8496tyu9c4dgade814)

Goal:
Make a jupyter notebook (week6.ipynb) that automates the steps described above.
In the end, the notebook will be run in github actions with:
    run: |
          cd week6
          jupyter nbconvert --execute week6.ipynb --to markdown --output week6
          cat week6.md

I don't know if I'll need to install conda to get this working, and I don't know if it's possible to get conda working in github actions, so I may need a workaround.