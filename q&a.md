1. **PIPEX only works in linux?**
  - PIPEX has been developed using Linux due to its flexibility (e.g. possibility to allocate Swap memory if it is needed). It also can be installed as a docker image with some knowledge. UPDATED: not anymore! Windows compatibility added.
2. **But I have to write python commands in that TXT file and I don't know python...**
  - Actually there is no need to add additional commands: as it is stated in the documentation, page 4, PIPEX runs sequentally a series of commands written in the file `pipex_batch_list.txt`. What is needed is to add the data folder ́s path, the cell segmentation parameters and the markers to be used. There is a template example named `pipex_batch_list_EXAMPLE.txt`.
3. **I don't have a membrane marker, so it means I can't use PIPEX?**
  - PIPEX is highly versatile meaning that it is possible to cell segment cells using only the DAPI marker through Stardist.
4. **Well, PIPEX seems to be just Stardist through command line...**
  - There are significant differences between running Stardist cell segmentation (via python or in QuPath, for example) vs PIPEX:
    - Stardist provides cell segmentation using only the DAPI marker (and then QuPath has integrated an algorithm to expand the nuclei into the cell). This approach may be useful for rounded cells where the cell has similar shape than the nuclei itself; however this scenario is not likely to be the case in a large number of tissues with different cell types.
    - PIPEX allows to cell segment cells using the DAPI marker and optionally to improve the cell segmentation with a membrane marker. Steps to follow:
      1. Cell segmentation using the DAPI marker (Stardist)
      2. Cell expansion (without overlap)
      3. Another cell segmentation step using the membrane marker (e.g. CDH1, CTNNB1) through a watershed algorithm
      4. Join both cell segmentations: the nuclei segmentation has priority meaning that if one nucleus has been divided in the watershed but not in the nuclei segmentation, it will be kept intact. Afterwards, the priority is given to the membrane marker that will refine and cut the expanded nucleus to the membrane limit.
5. **Can you make it run automatically when CODEX has finished to process the images?**
  - That is an additional integration work that may be implemented if it is really needed.
6. **I do my own analysis in R, not a solution for me...**
  - PIPEX generates a CSV file after the cell segmentation with all the necessary information to perform downstream data analysis which may be done using any programming language.
7. **I don't like the cluster color selection shown in QuPath, how can I change it?**
  - When PIPEX generates the clustering, a column is added in the CSV file with the cluster_color. It is accessible and it may be modified.
8. **Not useful to me, I don't usually work with CODEX images...**
  - Actually, PIPEX has been developed with the vision to be used in any scenario where there are fluorescence images that need to be cell segmented and analysed. It really just requires the name of the marker in the filename, nothing else.
9. **But will it work with Cytomap?**
  - PIPEX generates a CSV file after the cell segmentation with all the necessary information to perform downstream data analysis which may be done using any analysis tool.
10. **I don't have a powerful laptop, I'm sure PIPEX will not work in it...**
  - At the moment, PIPEX has been tested in not so powerful linux laptops with and without Swap memory added, with very satisfactory results. UPDATED: you can even modify the max_resolution through a command to instruct PIPEX to work downsampling your images for segmentation calculations.
11. **What if we want PIPEX to generate another/more plots in the analysis results?**
  - PIPEX has a python script named `analysis.py` which may be modified to add additional code for customized analysis.
12. **I can see the resulting images are in JPG. Aren't those bad? They should be TIFF...**
  - The resulting images are generated for QC purposes providing good enough information in JPG without the need of using the extra space required by TIFF format.
13. **Plain Stardist cell segmentation is better for my project, PIPEX is trying to be too smart. But my Stardist crashes with large images, what can I do?**
  - PIPEX allows to run basic Stardist alrgorithm with the advantage of smart memory usage to handle large images.
14. **Who will be responsible for PIPEX if we run into problems?**
  - Frederic Ballllosera et al.
15. **I would do it differently and, honestly, better by using [INSERT BUZZWORD HERE] and the latest [INSERT TRENDY NEW FRAMEWORK HERE]...**
  - That sounds very good! PIPEX was delivered in the requested time window to be able to analyse and present the ESPACE project data without delay: it was developed targeting the pressing need to have a better cell segmentation than the one provided by Akoya and offer a generalistic solution that covers other people's scenarios (including an integrated workflow with the data analysis). As always, the more tools the lab has, the merrier everyone!

